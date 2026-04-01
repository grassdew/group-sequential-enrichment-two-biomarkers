# =============================================================================
# Adaptive enrichment trial design functions (TWO continuous biomarkers)
# Biomarkers: x1, x2 ~ Uniform(0,1)
# Enrichment pattern output: JOINT 20x20 bin counts for ENROLLED subjects only
# (No screening counts are returned; screening is still tracked internally to
#  enforce max_screen and for GSE-F projection futility.)
# =============================================================================

# ---- Enrollment helper: enroll subjects based on PBI -------------------------
enroll_subjects <- function(model, subjects_to_screen, threshold, enrolled_indices,
                            enrolled_indices_current_stage, max_enrollments, epsilon) {
  pbis <- calculate_pbi(model, subjects_to_screen, epsilon)
  eligible_indices <- which(pbis > threshold)

  n_take <- min(length(eligible_indices), max_enrollments)
  if (n_take > 0) {
    take <- if (n_take < length(eligible_indices)) sample(eligible_indices, n_take) else eligible_indices
    enrolled_subjects <- subjects_to_screen[take, ]
    enrolled_indices  <- c(enrolled_indices, enrolled_subjects$subject_id)
    enrolled_indices_current_stage <- c(enrolled_indices_current_stage, enrolled_subjects$subject_id)
  } else {
    enrolled_subjects <- NULL
  }

  list(
    enrolled_subjects = enrolled_subjects,
    enrolled_indices = enrolled_indices,
    enrolled_indices_current_stage = enrolled_indices_current_stage
  )
}

# ---- Update regression model -------------------------------------------------
# Uses a tensor-product smooth surface for (x1, x2), plus a treatment-specific
# surface via a numeric 'trt' indicator (0=Control, 1=Treatment).
update_model <- function(cumulative_data) {

  if (!("trt" %in% names(cumulative_data))) {
    cumulative_data$trt <- as.numeric(cumulative_data$group == "Treatment")
  }

  updated_model <- stan_gamm4(
    formula = y ~ group + s(x1, x2) + s(x1, x2, by = trt),
    data    = cumulative_data,
    prior   = normal(0, 10),
    family  = gaussian(),
    chains  = 4,
    iter    = 2000,
    cores   = 1
  )

  updated_model
}

# ---- Personalized benefit index (PBI) ---------------------------------------
# epsilon: minimal clinically significant difference
calculate_pbi <- function(model, new_subject_data, epsilon) {
  if (is.null(model) || is.null(new_subject_data) || nrow(new_subject_data) == 0) {
    return(numeric(0))
  }

  # Ensure factor levels are consistent for posterior_predict
  levs <- c("Control", "Treatment")

  nd_treat <- new_subject_data
  nd_treat$group <- factor("Treatment", levels = levs)
  nd_treat$trt   <- 1

  nd_ctrl <- new_subject_data
  nd_ctrl$group <- factor("Control", levels = levs)
  nd_ctrl$trt   <- 0

  posterior_treatment <- posterior_predict(model, newdata = nd_treat)
  posterior_control   <- posterior_predict(model, newdata = nd_ctrl)

  # PBI = P( Y_T - Y_C > epsilon )
  colMeans((posterior_treatment - posterior_control) > epsilon)
}

# ---- Adaptive threshold ------------------------------------------------------
adaptive_threshold <- function(n_k, N, d, g) {
  d * (n_k / N)^g
}

# ---- t-statistic at each interim stage --------------------------------------
calculate_t_statistic <- function(interim_data) {
  treatment_y <- interim_data$y[interim_data$group == "Treatment"]
  control_y   <- interim_data$y[interim_data$group == "Control"]
  n_t <- length(treatment_y)
  n_c <- length(control_y)
  if (n_t < 2 || n_c < 2) return(0)

  mean_t <- mean(treatment_y)
  mean_c <- mean(control_y)
  var_t  <- var(treatment_y)
  var_c  <- var(control_y)

  den <- sqrt(var_t / n_t + var_c / n_c)
  if (!is.finite(den) || den == 0) return(0)
  (mean_t - mean_c) / den
}

# ---- Stratified test statistic (Simon & Simon, 2013) ------------------------
calculate_stratified_statistic <- function(test_statistics, incremental_sample_sizes, k) {
  if (k <= 0) return(0)
  idx <- seq_len(k)
  n_l <- incremental_sample_sizes[idx]
  n_l_sum <- sum(n_l)
  if (n_l_sum <= 0) return(0)
  weights <- sqrt(n_l / n_l_sum)
  sum(weights * test_statistics[idx])
}

# ---- Group sequential stopping boundaries -----------------------------------
stopping_bounds <- function(alpha, beta, type_alpha, type_beta, informationRates) {
  gsd <- getDesignGroupSequential(
    typeOfDesign = type_alpha, alpha = alpha,
    typeBetaSpending = type_beta, bindingFutility = FALSE,
    informationRates = informationRates,
    sided = 1, beta = beta
  )

  list(
    efficacy_bound = gsd$criticalValues,
    futility_bound = c(gsd$futilityBounds, gsd$criticalValues[length(gsd$criticalValues)])
  )
}

# ---- 2D (joint) bin counts for two continuous biomarkers --------------------
get_bin_indices <- function(x_values, bin_breaks) {
  as.numeric(cut(x_values, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE))
}

count_subjects_in_joint_bins <- function(x1_values, x2_values, bin_breaks) {
  b1 <- get_bin_indices(x1_values, bin_breaks)
  b2 <- get_bin_indices(x2_values, bin_breaks)
  n_bins <- length(bin_breaks) - 1L

  # Map (b1,b2) -> 1..(n_bins^2)
  joint_idx <- (b1 - 1L) * n_bins + b2
  tabulate(joint_idx, nbins = n_bins^2)
}

# ---- Projection futility helpers (used in GSE-F) ----------------------------
estimate_acceptance_rate <- function(model, epsilon, threshold, new_subjects) {
  if (is.null(model) || is.null(new_subjects) || nrow(new_subjects) == 0) return(NA_real_)
  pbi <- calculate_pbi(model, new_subjects, epsilon)
  mean(pbi > threshold)
}

projected_total_screens_to_finish <- function(total_enrolled, total_screened,
                                              max_sample_size, acceptance_rate) {
  remaining_enroll <- max(0, max_sample_size - total_enrolled)
  if (is.na(acceptance_rate) || acceptance_rate <= 0) return(Inf)
  expected_additional_screens <- ceiling(remaining_enroll / acceptance_rate)
  total_screened + expected_additional_screens
}

# =============================================================================
# Design 1: GSE (group-seq + enrichment + early stopping)
# =============================================================================
gse <- function(trial_data, max_screen, max_sample_size, info_rates,
                alpha, beta, type_alpha, type_beta, epsilon, d, g) {

  total_screened <- 0
  total_enrolled <- 0
  enrolled_indices <- c()
  enrolled_indices_current_stage <- c()

  model <- NULL
  test_statistics <- c()
  incremental_sample_sizes <- c()
  interim_looks <- c()
  early_stop <- FALSE
  conclusion <- NA
  stop_loop <- FALSE

  # bins (joint 20x20 over [0,1]^2)
  n_bins <- 20
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)
  n_joint <- n_bins^2

  safe_joint_counts <- function(id_vec) {
    if (length(id_vec) == 0L) {
      integer(n_joint)
    } else {
      x1_vals <- trial_data[id_vec, "x1", drop = TRUE]
      x2_vals <- trial_data[id_vec, "x2", drop = TRUE]
      count_subjects_in_joint_bins(x1_vals, x2_vals, bin_breaks)
    }
  }

  interim_enrollments <- c()
  enrollment_joint_counts_list <- list()

  start_idx <- 1
  threshold <- 0
  screen_cohort <- max_sample_size * info_rates[1]

  while (total_enrolled < max_sample_size && total_screened < max_screen && !stop_loop) {

    next_milestone_index <- length(interim_looks) + 1
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    remaining_to_milestone <- next_milestone - total_enrolled

    screen_count <- min(max_screen - total_screened, screen_cohort)
    end_idx <- min(start_idx + screen_count - 1, nrow(trial_data))
    new_subjects <- trial_data[start_idx:end_idx, ]

    # update screened count (internal only)
    total_screened <- total_screened + nrow(new_subjects)

    # enroll
    if (total_enrolled < info_rates[1] * max_sample_size) {
      num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
      subjects_to_enroll <- new_subjects[1:num_to_enroll, ]

      enrolled_indices <- c(enrolled_indices, subjects_to_enroll$subject_id)
      enrolled_indices_current_stage <- c(enrolled_indices_current_stage, subjects_to_enroll$subject_id)
      total_enrolled <- length(enrolled_indices)
    } else {
      res <- enroll_subjects(model, new_subjects, threshold,
                             enrolled_indices, enrolled_indices_current_stage,
                             remaining_to_milestone, epsilon)
      enrolled_indices <- res$enrolled_indices
      enrolled_indices_current_stage <- res$enrolled_indices_current_stage
      total_enrolled <- length(enrolled_indices)
    }

    # interim analysis
    if (total_enrolled >= next_milestone) {

      interim_looks <- c(interim_looks, total_enrolled)

      # update model + threshold
      cumulative_data <- trial_data[enrolled_indices, ]
      model <- update_model(cumulative_data)
      threshold <- adaptive_threshold(total_enrolled, max_sample_size, d, g)

      # incremental t-stat
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))

      # stopping decision
      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, length(test_statistics))
      bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)

      if (total_enrolled < max_sample_size) {
        if (stratified_stat > bounds$efficacy_bound[length(interim_looks)]) {
          early_stop <- TRUE; conclusion <- "Efficacy"
        } else if (stratified_stat < bounds$futility_bound[length(interim_looks)]) {
          early_stop <- TRUE; conclusion <- "Futility"
        }
      } else {
        conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[length(interim_looks)],
                             "Efficacy", "Futility")
      }

      # record interim outputs (enrollment only)
      interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
      enrollment_joint_counts_list[[length(interim_looks)]] <- safe_joint_counts(enrolled_indices_current_stage)

      # reset stage buffer
      enrolled_indices_current_stage <- c()

      if (early_stop) stop_loop <- TRUE
    }

    start_idx <- end_idx + 1
    if (start_idx > nrow(trial_data)) break
  }

  # Final look if needed (no early stop and not all planned looks)
  max_interims <- length(info_rates)
  if (!early_stop && total_enrolled > 0 && length(interim_looks) < max_interims) {
    interim_looks <- c(interim_looks, total_enrolled)

    if (length(enrolled_indices_current_stage) > 0) {
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
    }

    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")

    interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
    enrollment_joint_counts_list[[length(interim_looks)]] <- safe_joint_counts(enrolled_indices_current_stage)

    enrolled_indices_current_stage <- integer(0)
  }

  # status (internal screening still matters)
  if (total_enrolled == max_sample_size) {
    status <- "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    status <- "Screening limit reached"
  } else if (early_stop) {
    status <- "Early stopping"
  } else {
    status <- "Stopped"
  }

  # Guarantee a final decision if missing
  if (is.na(conclusion) && total_enrolled > 0) {
    have_any_look <- length(test_statistics) > 0
    if (!have_any_look) {
      # analyze all enrolled as one look
      interim_data <- trial_data[enrolled_indices, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(t_stat)
      incremental_sample_sizes <- c(nrow(interim_data))
      interim_looks <- c(interim_looks, total_enrolled)

      interim_enrollments <- c(interim_enrollments, length(enrolled_indices))
      enrollment_joint_counts_list[[length(interim_looks)]] <- safe_joint_counts(enrolled_indices)
    }

    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
  }

  # pad interim enrollments to |info_rates|
  interim_enrollments_padded <- c(interim_enrollments, rep(NA, max_interims - length(interim_enrollments)))
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)

  # flatten joint bin counts
  flattened_enroll_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(enrollment_joint_counts_list)) {
      v <- enrollment_joint_counts_list[[k]]
      for (i in seq_len(n_bins)) {
        for (j in seq_len(n_bins)) {
          idx <- (i - 1L) * n_bins + j
          nm  <- paste0("enrollments_x1bin", i, "_x2bin", j, "_interim_", k)
          flattened_enroll_counts[[nm]] <- v[idx]
        }
      }
    } else {
      for (i in seq_len(n_bins)) {
        for (j in seq_len(n_bins)) {
          nm <- paste0("enrollments_x1bin", i, "_x2bin", j, "_interim_", k)
          flattened_enroll_counts[[nm]] <- NA
        }
      }
    }
  }

  c(
    list(
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(flattened_enroll_counts)
  )
}

# =============================================================================
# Design 2: GSE-F (GSE + projection futility check)  [aka gse_plus]
# =============================================================================
gse_plus <- function(trial_data, max_screen, max_sample_size, info_rates,
                     alpha, beta, type_alpha, type_beta, epsilon, d, g, delta) {

  total_screened <- 0
  total_enrolled <- 0
  enrolled_indices <- c()
  enrolled_indices_current_stage <- c()

  model <- NULL
  test_statistics <- c()
  incremental_sample_sizes <- c()
  interim_looks <- c()
  early_stop <- FALSE
  conclusion <- NA

  n_bins <- 20
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)
  n_joint <- n_bins^2

  safe_joint_counts <- function(id_vec) {
    if (length(id_vec) == 0L) {
      integer(n_joint)
    } else {
      x1_vals <- trial_data[id_vec, "x1", drop = TRUE]
      x2_vals <- trial_data[id_vec, "x2", drop = TRUE]
      count_subjects_in_joint_bins(x1_vals, x2_vals, bin_breaks)
    }
  }

  interim_enrollments <- c()
  enrollment_joint_counts_list <- list()

  # commit partial stage as an interim snapshot (used for projection futility)
  commit_current_stage <- function() {
    interim_looks <<- c(interim_looks, total_enrolled)
    interim_enrollments <<- c(interim_enrollments, length(enrolled_indices_current_stage))
    enrollment_joint_counts_list[[length(interim_looks)]] <<- safe_joint_counts(enrolled_indices_current_stage)
  }

  start_idx <- 1
  threshold <- 0
  screen_cohort <- max_sample_size * info_rates[1]

  while (total_enrolled < max_sample_size && total_screened < max_screen) {

    next_milestone_index <- length(interim_looks) + 1
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    remaining_to_milestone <- next_milestone - total_enrolled

    screen_count <- min(max_screen - total_screened, screen_cohort)
    end_idx <- min(start_idx + screen_count - 1, nrow(trial_data))
    new_subjects <- trial_data[start_idx:end_idx, ]

    # ---- projection futility (can trigger mid-stage) -------------------------
    if (!is.null(model) && total_enrolled >= info_rates[1] * max_sample_size) {
      thr_now  <- adaptive_threshold(total_enrolled, max_sample_size, d, g)
      acc_rate <- estimate_acceptance_rate(model, epsilon, thr_now, new_subjects)
      proj_tot <- projected_total_screens_to_finish(total_enrolled, total_screened,
                                                    max_sample_size, acc_rate)
      if (proj_tot > (max_screen + delta)) {
        early_stop <- TRUE
        conclusion <- "Futility"
        if (length(enrolled_indices_current_stage) > 0) {
          commit_current_stage()
        }
        break
      }
    }

    total_screened <- total_screened + nrow(new_subjects)

    # enroll
    if (total_enrolled < info_rates[1] * max_sample_size) {
      num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
      subjects_to_enroll <- new_subjects[1:num_to_enroll, ]
      enrolled_indices <- c(enrolled_indices, subjects_to_enroll$subject_id)
      enrolled_indices_current_stage <- c(enrolled_indices_current_stage, subjects_to_enroll$subject_id)
      total_enrolled <- length(enrolled_indices)
    } else {
      res <- enroll_subjects(model, new_subjects, threshold,
                             enrolled_indices, enrolled_indices_current_stage,
                             remaining_to_milestone, epsilon)
      enrolled_indices <- res$enrolled_indices
      enrolled_indices_current_stage <- res$enrolled_indices_current_stage
      total_enrolled <- length(enrolled_indices)
    }

    # milestone look
    if (total_enrolled >= next_milestone) {

      cumulative_data <- trial_data[enrolled_indices, ]
      model <- update_model(cumulative_data)
      threshold <- adaptive_threshold(total_enrolled, max_sample_size, d, g)

      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))

      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, length(test_statistics))
      bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)

      # record interim outputs BEFORE clearing buffer
      interim_looks <- c(interim_looks, total_enrolled)
      interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
      enrollment_joint_counts_list[[length(interim_looks)]] <- safe_joint_counts(enrolled_indices_current_stage)

      enrolled_indices_current_stage <- c()

      if (total_enrolled < max_sample_size) {
        if (stratified_stat > bounds$efficacy_bound[length(interim_looks)]) {
          early_stop <- TRUE; conclusion <- "Efficacy"; break
        } else if (stratified_stat < bounds$futility_bound[length(interim_looks)]) {
          early_stop <- TRUE; conclusion <- "Futility"; break
        }
      } else {
        conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[length(interim_looks)],
                             "Efficacy", "Futility")
        break
      }
    }

    start_idx <- end_idx + 1
    if (start_idx > nrow(trial_data)) break
  }

  # Final look if needed
  max_interims <- length(info_rates)
  if (!early_stop && total_enrolled > 0 && length(interim_looks) < max_interims) {
    interim_looks <- c(interim_looks, total_enrolled)
    if (length(enrolled_indices_current_stage) > 0) {
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
    }
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")

    interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
    enrollment_joint_counts_list[[length(interim_looks)]] <- safe_joint_counts(enrolled_indices_current_stage)
    enrolled_indices_current_stage <- integer(0)
  }

  if (total_enrolled == max_sample_size) {
    status <- "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    status <- "Screening limit reached"
  } else if (early_stop) {
    status <- "Early stopping"
  } else {
    status <- "Stopped"
  }

  # pad + flatten
  interim_enrollments_padded <- c(interim_enrollments, rep(NA, max_interims - length(interim_enrollments)))
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)

  flattened_enroll_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(enrollment_joint_counts_list)) {
      v <- enrollment_joint_counts_list[[k]]
      for (i in seq_len(n_bins)) for (j in seq_len(n_bins)) {
        idx <- (i - 1L) * n_bins + j
        nm <- paste0("enrollments_x1bin", i, "_x2bin", j, "_interim_", k)
        flattened_enroll_counts[[nm]] <- v[idx]
      }
    } else {
      for (i in seq_len(n_bins)) for (j in seq_len(n_bins)) {
        nm <- paste0("enrollments_x1bin", i, "_x2bin", j, "_interim_", k)
        flattened_enroll_counts[[nm]] <- NA
      }
    }
  }

  c(
    list(
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(flattened_enroll_counts)
  )
}

# =============================================================================
# Design 3: AE (enrichment; final decision only; no early stopping)
# =============================================================================
ae <- function(trial_data, max_screen, max_sample_size, info_rates,
               alpha, beta, type_alpha, type_beta, epsilon, d, g) {

  total_screened <- 0
  total_enrolled <- 0
  enrolled_indices <- c()
  enrolled_indices_current_stage <- c()

  model <- NULL
  test_statistics <- c()
  incremental_sample_sizes <- c()
  interim_looks <- c()
  conclusion <- NA

  n_bins <- 20
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)
  n_joint <- n_bins^2

  safe_joint_counts <- function(id_vec) {
    if (length(id_vec) == 0L) {
      integer(n_joint)
    } else {
      x1_vals <- trial_data[id_vec, "x1", drop = TRUE]
      x2_vals <- trial_data[id_vec, "x2", drop = TRUE]
      count_subjects_in_joint_bins(x1_vals, x2_vals, bin_breaks)
    }
  }

  interim_enrollments <- c()
  enrollment_joint_counts_list <- list()

  start_idx <- 1
  threshold <- 0
  screen_cohort <- max_sample_size * info_rates[1]

  while (total_enrolled < max_sample_size && total_screened < max_screen) {

    next_milestone_index <- length(interim_looks) + 1
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    remaining_to_milestone <- next_milestone - total_enrolled

    screen_count <- min(max_screen - total_screened, screen_cohort)
    end_idx <- min(start_idx + screen_count - 1, nrow(trial_data))
    new_subjects <- trial_data[start_idx:end_idx, ]

    total_screened <- total_screened + nrow(new_subjects)

    if (total_enrolled < info_rates[1] * max_sample_size) {
      num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
      subjects_to_enroll <- new_subjects[1:num_to_enroll, ]
      enrolled_indices <- c(enrolled_indices, subjects_to_enroll$subject_id)
      enrolled_indices_current_stage <- c(enrolled_indices_current_stage, subjects_to_enroll$subject_id)
      total_enrolled <- length(enrolled_indices)
    } else {
      res <- enroll_subjects(model, new_subjects, threshold,
                             enrolled_indices, enrolled_indices_current_stage,
                             remaining_to_milestone, epsilon)
      enrolled_indices <- res$enrolled_indices
      enrolled_indices_current_stage <- res$enrolled_indices_current_stage
      total_enrolled <- length(enrolled_indices)
    }

    if (total_enrolled >= next_milestone) {

      interim_looks <- c(interim_looks, total_enrolled)

      cumulative_data <- trial_data[enrolled_indices, ]
      model <- update_model(cumulative_data)
      threshold <- adaptive_threshold(total_enrolled, max_sample_size, d, g)

      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))

      interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
      enrollment_joint_counts_list[[length(interim_looks)]] <- safe_joint_counts(enrolled_indices_current_stage)

      enrolled_indices_current_stage <- c()
    }

    start_idx <- end_idx + 1
    if (start_idx > nrow(trial_data)) break
  }

  # Final decision (last look if needed)
  max_interims <- length(info_rates)
  if (total_enrolled > 0 && length(interim_looks) < max_interims) {
    interim_looks <- c(interim_looks, total_enrolled)
    if (length(enrolled_indices_current_stage) > 0) {
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
    }
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")

    interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
    enrollment_joint_counts_list[[length(interim_looks)]] <- safe_joint_counts(enrolled_indices_current_stage)
    enrolled_indices_current_stage <- integer(0)
  }

  if (total_enrolled == max_sample_size) {
    status <- "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    status <- "Screening limit reached"
  } else {
    status <- "Stopped"
  }

  if (is.na(conclusion) && total_enrolled > 0) {
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)

    if (length(test_statistics) > 0) {
      k <- length(test_statistics)
      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
      conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
    } else {
      # Edge case: no interim looks were recorded; analyze all enrolled as one look
      interim_data <- trial_data[enrolled_indices, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(t_stat)
      incremental_sample_sizes <- c(nrow(interim_data))
      conclusion <- ifelse(t_stat > bounds$efficacy_bound[1], "Efficacy", "Futility")
    }
  }


  interim_enrollments_padded <- c(interim_enrollments, rep(NA, max_interims - length(interim_enrollments)))
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)

  flattened_enroll_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(enrollment_joint_counts_list)) {
      v <- enrollment_joint_counts_list[[k]]
      for (i in seq_len(n_bins)) for (j in seq_len(n_bins)) {
        idx <- (i - 1L) * n_bins + j
        nm <- paste0("enrollments_x1bin", i, "_x2bin", j, "_interim_", k)
        flattened_enroll_counts[[nm]] <- v[idx]
      }
    } else {
      for (i in seq_len(n_bins)) for (j in seq_len(n_bins)) {
        nm <- paste0("enrollments_x1bin", i, "_x2bin", j, "_interim_", k)
        flattened_enroll_counts[[nm]] <- NA
      }
    }
  }

  c(
    list(
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(flattened_enroll_counts)
  )
}

# =============================================================================
# Design 4: GSD (group sequential; no enrichment)
# =============================================================================
gsd <- function(trial_data, max_screen, max_sample_size, info_rates,
                alpha, beta, type_alpha, type_beta) {

  total_screened <- 0
  total_enrolled <- 0
  enrolled_indices <- c()
  enrolled_indices_current_stage <- c()

  test_statistics <- c()
  incremental_sample_sizes <- c()
  interim_looks <- c()
  early_stop <- FALSE
  conclusion <- NA
  stop_loop <- FALSE

  n_bins <- 20
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)
  n_joint <- n_bins^2

  safe_joint_counts <- function(id_vec) {
    if (length(id_vec) == 0L) {
      integer(n_joint)
    } else {
      x1_vals <- trial_data[id_vec, "x1", drop = TRUE]
      x2_vals <- trial_data[id_vec, "x2", drop = TRUE]
      count_subjects_in_joint_bins(x1_vals, x2_vals, bin_breaks)
    }
  }

  interim_enrollments <- c()
  enrollment_joint_counts_list <- list()

  start_idx <- 1
  screen_cohort <- max_sample_size * info_rates[1]

  while (total_enrolled < max_sample_size && total_screened < max_screen && !stop_loop) {

    next_milestone_index <- length(interim_looks) + 1
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    remaining_to_milestone <- next_milestone - total_enrolled

    screen_count <- min(max_screen - total_screened, screen_cohort)
    end_idx <- min(start_idx + screen_count - 1, nrow(trial_data))
    new_subjects <- trial_data[start_idx:end_idx, ]

    total_screened <- total_screened + nrow(new_subjects)

    # enroll without enrichment
    num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
    subjects_to_enroll <- new_subjects[1:num_to_enroll, ]

    enrolled_indices <- c(enrolled_indices, subjects_to_enroll$subject_id)
    enrolled_indices_current_stage <- c(enrolled_indices_current_stage, subjects_to_enroll$subject_id)
    total_enrolled <- length(enrolled_indices)

    if (total_enrolled >= next_milestone) {

      interim_looks <- c(interim_looks, total_enrolled)

      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))

      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, length(test_statistics))
      bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)

      if (total_enrolled < max_sample_size) {
        if (stratified_stat > bounds$efficacy_bound[length(interim_looks)]) {
          early_stop <- TRUE; conclusion <- "Efficacy"
        } else if (stratified_stat < bounds$futility_bound[length(interim_looks)]) {
          early_stop <- TRUE; conclusion <- "Futility"
        }
      } else {
        conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[length(interim_looks)],
                             "Efficacy", "Futility")
      }

      interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
      enrollment_joint_counts_list[[length(interim_looks)]] <- safe_joint_counts(enrolled_indices_current_stage)

      enrolled_indices_current_stage <- c()

      if (early_stop) stop_loop <- TRUE
    }

    start_idx <- end_idx + 1
    if (start_idx > nrow(trial_data)) break
  }

  # Final look if needed
  max_interims <- length(info_rates)
  if (!early_stop && total_enrolled > 0 && length(interim_looks) < max_interims) {
    interim_looks <- c(interim_looks, total_enrolled)
    if (length(enrolled_indices_current_stage) > 0) {
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
    }
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")

    interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
    enrollment_joint_counts_list[[length(interim_looks)]] <- safe_joint_counts(enrolled_indices_current_stage)
    enrolled_indices_current_stage <- integer(0)
  }

  if (total_enrolled == max_sample_size) {
    status <- "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    status <- "Screening limit reached"
  } else if (early_stop) {
    status <- "Early stopping"
  } else {
    status <- "Stopped"
  }

  if (is.na(conclusion) && total_enrolled > 0) {
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)

    if (length(test_statistics) > 0) {
      k <- length(test_statistics)
      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
      conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
    } else {
      # Edge case: no interim looks were recorded; analyze all enrolled as one look
      interim_data <- trial_data[enrolled_indices, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(t_stat)
      incremental_sample_sizes <- c(nrow(interim_data))
      conclusion <- ifelse(t_stat > bounds$efficacy_bound[1], "Efficacy", "Futility")
    }
  }


  interim_enrollments_padded <- c(interim_enrollments, rep(NA, max_interims - length(interim_enrollments)))
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)

  flattened_enroll_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(enrollment_joint_counts_list)) {
      v <- enrollment_joint_counts_list[[k]]
      for (i in seq_len(n_bins)) for (j in seq_len(n_bins)) {
        idx <- (i - 1L) * n_bins + j
        nm <- paste0("enrollments_x1bin", i, "_x2bin", j, "_interim_", k)
        flattened_enroll_counts[[nm]] <- v[idx]
      }
    } else {
      for (i in seq_len(n_bins)) for (j in seq_len(n_bins)) {
        nm <- paste0("enrollments_x1bin", i, "_x2bin", j, "_interim_", k)
        flattened_enroll_counts[[nm]] <- NA
      }
    }
  }

  c(
    list(
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(flattened_enroll_counts)
  )
}

# =============================================================================
# Trial simulation wrapper (TWO biomarkers)
# =============================================================================
trial_simulation <- function(sim_index, design_type,
                             max_screen,
                             max_sample_size,
                             info_rates,
                             effect_size,
                             alpha, beta,
                             type_alpha, type_beta,
                             epsilon, d, g, delta,
                             scenario_index) {

  set.seed(sim_index)

  data <- data.frame(
    subject_id = 1:max_screen,
    x1 = runif(max_screen, min = 0, max = 1),
    x2 = runif(max_screen, min = 0, max = 1),
    group = rbinom(max_screen, size = 1, prob = 0.5)
  )

  # outcome model params
  beta0 <- 5
  beta1 <- 0.5
  sigma <- 0.1

  baseline  <- 3
  amplitude <- 2
  multiplier <- 1.2

  noise <- rnorm(max_screen, mean = 0, sd = sigma)

  # helpers based on both biomarkers
  x_avg <- (data$x1 + data$x2) / 2
  u_avg <- (4 * (data$x1 - 0.5)^2 + 4 * (data$x2 - 0.5)^2) / 2         # U-shape score in [0,1]
  inv_avg <- (4 * data$x1 * (1 - data$x1) + 4 * data$x2 * (1 - data$x2)) / 2  # inverted-U in [0,1]

  if (scenario_index == 1) {
    # Null: no treatment effect
    effect <- 0
    y_control <- beta0 + beta1 * x_avg + noise
    y_treat   <- y_control + effect
    data$y <- ifelse(data$group == 1, y_treat, y_control)

  } else if (scenario_index == 2) {
    # Constant effect
    effect <- effect_size
    y_control <- beta0 + beta1 * x_avg + noise
    y_treat   <- y_control + effect
    data$y <- ifelse(data$group == 1, y_treat, y_control)

  } else if (scenario_index == 3) {
    # Effect increases with average biomarker level (monotone)
    effect <- effect_size * x_avg * multiplier
    y_control <- beta0 + beta1 * x_avg + noise
    y_treat   <- y_control + effect
    data$y <- ifelse(data$group == 1, y_treat, y_control)

  } else if (scenario_index == 4) {
    # U-shape: effect largest at extremes (in either biomarker)
    effect <- effect_size * u_avg * multiplier
    y_control <- baseline + amplitude * u_avg + noise
    y_treat   <- y_control + effect
    data$y <- ifelse(data$group == 1, y_treat, y_control)

  } else if (scenario_index == 5) {
    # Inverted-U: effect largest near middle (both biomarkers)
    effect <- effect_size * inv_avg * multiplier
    y_control <- baseline + amplitude * inv_avg + noise
    y_treat   <- y_control + effect
    data$y <- ifelse(data$group == 1, y_treat, y_control)
  }

  data$group <- factor(data$group, levels = c(0, 1), labels = c("Control", "Treatment"))
  data$trt   <- as.numeric(data$group == "Treatment")

  if (design_type == "GSE") {
    result <- gse(data, max_screen, max_sample_size, info_rates,
                  alpha, beta, type_alpha, type_beta, epsilon, d, g)

  } else if (design_type == "AE") {
    result <- ae(data, max_screen, max_sample_size, info_rates,
                 alpha, beta, type_alpha, type_beta, epsilon, d, g)

  } else if (design_type == "GSD") {
    result <- gsd(data, max_screen, max_sample_size, info_rates,
                  alpha, beta, type_alpha, type_beta)

  } else if (design_type == "GSE-F") {
    result <- gse_plus(data, max_screen, max_sample_size, info_rates,
                       alpha, beta, type_alpha, type_beta, epsilon, d, g, delta)
  }
  result
}
