# =============================================================================
# Design functions: 1 binary biomarker (x1) + 1 continuous biomarker (x2)
# Designs compared: GSE, AE, GSD, GSE-F
#
# Key features:
# - Adaptive enrichment via PBI from a Bayesian GAM
# - Group sequential testing (GSE / GSD), plus feasibility-based futility (GSE-F)
#
# v4 changes (Jan 2026):
#   1) OUTPUT: Joint (x1, x2-bin) counts only.
#      - Marginal x1 counts and marginal x2-bin counts are removed (derivable from joint counts).
#   2) FIX: Ensure consistency between totals and interim-stage counts.
#      - If the trial stops before the next planned interim milestone (e.g., screening limit reached),
#        the remaining current-stage buffers are committed as the next interim stage.
#
# Joint count output columns (b=1..20, k=1..K):
#   - enrollments_x1neg_bin{b}_interim_{k}
#   - enrollments_x1pos_bin{b}_interim_{k}
#   - screenings_x1neg_bin{b}_interim_{k}
#   - screenings_x1pos_bin{b}_interim_{k}
#
# Stage totals output columns:
#   - interim_enrollments_k
#   - interim_screenings_k
#
# Notes
# - x1 is binary {0,1} (1 = biomarker-positive)
# - x2 is continuous on [0,1], binned into 20 bins for tracking
# - Outcome model (Bayesian GAM):
#     y ~ x1 + group + group:x1 + s(x2) + s(x2, by = trt) + s(x2, by = trt_x1)
#   where trt = 1(group=="Treatment"), trt_x1 = trt*x1
#   This includes a treatment-by-x1-by-x2 interaction via s(x2, by=trt_x1),
#   without including a baseline x1:x2 interaction.
#
# Required packages (load in your simulation script): rstanarm, rpact
# =============================================================================

# ---- Enrollment based on PBI -------------------------------------------------
enroll_subjects <- function(model, subjects_to_screen, threshold,
                            enrolled_indices, enrolled_indices_current_stage,
                            max_enrollments, epsilon) {
  pbis <- calculate_pbi(model, subjects_to_screen, epsilon)
  eligible_indices <- which(pbis > threshold)

  n_take <- min(length(eligible_indices), max_enrollments)
  if (n_take > 0) {
    take <- if (n_take < length(eligible_indices)) sample(eligible_indices, n_take) else eligible_indices
    enrolled_subjects <- subjects_to_screen[take, , drop = FALSE]
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

# ---- Model update (Bayesian GAM) --------------------------------------------
update_model <- function(cumulative_data) {
  dat <- cumulative_data
  dat$group <- factor(dat$group, levels = c("Control", "Treatment"))
  dat$trt    <- as.numeric(dat$group == "Treatment")
  dat$trt_x1 <- dat$trt * dat$x1

  rstanarm::stan_gamm4(
    formula = y ~ x1 + group + group:x1 +
      s(x2) + s(x2, by = trt) + s(x2, by = trt_x1),
    data    = dat,
    prior   = rstanarm::normal(0, 10^4),
    family  = gaussian(),
    chains  = 4,
    iter    = 2000,
    cores   = 1
  )
}

# ---- PBI ---------------------------------------------------------------------
calculate_pbi <- function(model, new_subject_data, epsilon) {
  if (is.null(model) || is.null(new_subject_data) || nrow(new_subject_data) == 0) {
    return(numeric(0))
  }

  make_newdata <- function(df, grp_label) {
    out <- df
    out$group <- factor(grp_label, levels = c("Control", "Treatment"))
    out$trt   <- as.numeric(out$group == "Treatment")
    out$trt_x1 <- out$trt * out$x1
    out
  }

  posterior_treatment <- rstanarm::posterior_predict(model,
    newdata = make_newdata(new_subject_data, "Treatment")
  )
  posterior_control <- rstanarm::posterior_predict(model,
    newdata = make_newdata(new_subject_data, "Control")
  )

  colMeans((posterior_treatment - posterior_control) > epsilon)
}

# ---- Adaptive threshold ------------------------------------------------------
adaptive_threshold <- function(n_k, N, d, g) {
  d * (n_k / N)^g
}

# ---- Incremental t-statistic -------------------------------------------------
calculate_t_statistic <- function(interim_data) {
  treatment_y <- interim_data$y[interim_data$group == "Treatment"]
  control_y   <- interim_data$y[interim_data$group == "Control"]
  n_t <- length(treatment_y)
  n_c <- length(control_y)

  if (n_t < 2 || n_c < 2) return(0)

  mean_t <- mean(treatment_y); mean_c <- mean(control_y)
  var_t  <- stats::var(treatment_y); var_c <- stats::var(control_y)
  den <- sqrt(var_t / n_t + var_c / n_c)
  if (!is.finite(den) || den == 0) return(0)

  (mean_t - mean_c) / den
}

# ---- Stratified statistic ----------------------------------------------------
calculate_stratified_statistic <- function(test_statistics, incremental_sample_sizes, k) {
  if (k <= 0) return(0)
  idx <- seq_len(k)
  n_l <- incremental_sample_sizes[idx]
  n_l_sum <- sum(n_l)
  if (n_l_sum <= 0) return(0)
  weights <- sqrt(n_l / n_l_sum)
  sum(weights * test_statistics[idx])
}

# ---- Group sequential bounds -------------------------------------------------
stopping_bounds <- function(alpha, beta, type_alpha, type_beta, informationRates) {
  gsd <- rpact::getDesignGroupSequential(
    typeOfDesign = type_alpha,
    alpha = alpha,
    typeBetaSpending = type_beta,
    bindingFutility = FALSE,
    informationRates = informationRates,
    sided = 1,
    beta = beta
  )

  list(
    efficacy_bound = gsd$criticalValues,
    futility_bound = c(gsd$futilityBounds, gsd$criticalValues[length(gsd$criticalValues)])
  )
}

# ---- Utilities for x2 bin indices -------------------------------------------
get_bin_indices <- function(x_values, bin_breaks) {
  as.numeric(cut(x_values, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE))
}

# ---- Joint (x1, x2-bin) counts ----------------------------------------------
# Returns a 2 x n_bins integer matrix:
#   row 1 = x1neg (x1==0)
#   row 2 = x1pos (x1==1)
safe_joint_counts_x1_x2bin <- function(trial_data, id_vec, bin_breaks) {
  n_bins <- length(bin_breaks) - 1L
  if (length(id_vec) == 0L) {
    return(matrix(0L, nrow = 2L, ncol = n_bins,
                  dimnames = list(c("x1neg", "x1pos"), paste0("bin", seq_len(n_bins)))))
  }

  x1 <- trial_data[id_vec, "x1", drop = TRUE]
  x2 <- trial_data[id_vec, "x2", drop = TRUE]
  bin_idx <- get_bin_indices(x2, bin_breaks)

  x1neg_counts <- tabulate(bin_idx[x1 == 0], nbins = n_bins)
  x1pos_counts <- tabulate(bin_idx[x1 == 1], nbins = n_bins)

  out <- rbind(x1neg_counts, x1pos_counts)
  rownames(out) <- c("x1neg", "x1pos")
  colnames(out) <- paste0("bin", seq_len(n_bins))
  out
}

# ---- Flatten joint lists to named scalars -----------------------------------
flatten_joint <- function(joint_list, max_interims, n_bins, prefix) {
  out <- list()
  for (k in seq_len(max_interims)) {
    for (b in seq_len(n_bins)) {
      nm_neg <- paste0(prefix, "_x1neg_bin", b, "_interim_", k)
      nm_pos <- paste0(prefix, "_x1pos_bin", b, "_interim_", k)
      if (k <= length(joint_list)) {
        out[[nm_neg]] <- joint_list[[k]]["x1neg", b]
        out[[nm_pos]] <- joint_list[[k]]["x1pos", b]
      } else {
        out[[nm_neg]] <- NA
        out[[nm_pos]] <- NA
      }
    }
  }
  out
}

# ---- Feasibility futility helpers (for GSE-F) --------------------------------
estimate_acceptance_rate <- function(model, epsilon, threshold, new_subjects) {
  if (is.null(model) || is.null(new_subjects) || nrow(new_subjects) == 0) return(NA_real_)
  pbi <- calculate_pbi(model, new_subjects, epsilon)
  if (length(pbi) == 0) return(NA_real_)
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
# DESIGN 1: GSE (adaptive enrichment + group sequential early stop)
# =============================================================================
gse <- function(trial_data, max_screen, max_sample_size, info_rates,
                alpha, beta, type_alpha, type_beta, epsilon, d, g) {

  bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)

  total_screened <- 0L
  total_enrolled <- 0L

  enrolled_indices <- integer(0)
  enrolled_indices_current_stage <- integer(0)
  screened_indices_current_stage <- integer(0)

  model <- NULL
  test_statistics <- numeric(0)
  incremental_sample_sizes <- integer(0)
  interim_looks <- integer(0)

  early_stop <- FALSE
  conclusion <- NA_character_
  stop_loop <- FALSE

  n_bins <- 20L
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)

  interim_enrollments <- integer(0)
  interim_screenings  <- integer(0)

  enroll_joint_list <- list()
  screen_joint_list <- list()

  start_idx <- 1L
  threshold <- 0
  screen_cohort <- ceiling(max_sample_size * info_rates[1])

  # Commit current-stage buffers as an interim stage
  commit_stage <- function() {
    interim_enrollments <<- c(interim_enrollments, length(enrolled_indices_current_stage))
    interim_screenings  <<- c(interim_screenings,  length(screened_indices_current_stage))
    enroll_joint_list[[length(interim_enrollments)]] <<-
      safe_joint_counts_x1_x2bin(trial_data, enrolled_indices_current_stage, bin_breaks)
    screen_joint_list[[length(interim_screenings)]]  <<-
      safe_joint_counts_x1_x2bin(trial_data, screened_indices_current_stage, bin_breaks)
    enrolled_indices_current_stage <<- integer(0)
    screened_indices_current_stage <<- integer(0)
  }

  while (total_enrolled < max_sample_size && total_screened < max_screen && !stop_loop) {

    next_milestone_index <- length(interim_looks) + 1L
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    remaining_to_milestone <- next_milestone - total_enrolled

    screen_count <- min(max_screen - total_screened, screen_cohort)
    end_idx <- min(start_idx + screen_count - 1L, nrow(trial_data))
    new_subjects <- trial_data[start_idx:end_idx, , drop = FALSE]

    total_screened <- total_screened + nrow(new_subjects)
    screened_indices_current_stage <- c(screened_indices_current_stage, new_subjects$subject_id)

    # enroll
    if (total_enrolled < info_rates[1] * max_sample_size) {
      num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
      subjects_to_enroll <- new_subjects[seq_len(num_to_enroll), , drop = FALSE]
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

    # interim look?
    if (total_enrolled >= next_milestone) {

      interim_looks <- c(interim_looks, total_enrolled)

      cumulative_data <- trial_data[enrolled_indices, , drop = FALSE]
      model <- update_model(cumulative_data)
      threshold <- adaptive_threshold(total_enrolled, max_sample_size, d, g)

      # incremental test statistic uses only current-stage enrolled subjects
      interim_data <- trial_data[enrolled_indices_current_stage, , drop = FALSE]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))

      k <- length(test_statistics)
      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)

      # record stage totals + joint distributions
      commit_stage()

      # stopping
      if (total_enrolled < max_sample_size) {
        if (stratified_stat > bounds$efficacy_bound[k]) {
          early_stop <- TRUE; conclusion <- "Efficacy"
        } else if (stratified_stat < bounds$futility_bound[k]) {
          early_stop <- TRUE; conclusion <- "Futility"
        }
      } else {
        conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
      }

      if (early_stop) stop_loop <- TRUE
    }

    start_idx <- end_idx + 1L
    if (start_idx > nrow(trial_data)) break
  }

  # FIX: commit leftover buffered stage (if any) so totals match sums of stage counts
  if (length(enrolled_indices_current_stage) > 0L || length(screened_indices_current_stage) > 0L) {
    commit_stage()
  }

  # Ensure conclusion if still missing (use last computed stratified stat if possible)
  if (is.na(conclusion) && length(test_statistics) > 0) {
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
  }

  status <- if (total_enrolled == max_sample_size) {
    "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    "Screening limit reached"
  } else if (early_stop) {
    "Early stopping"
  }

  max_interims <- length(info_rates)
  pad_int <- function(x) c(x, rep(NA, max_interims - length(x)))

  interim_enrollments_padded <- pad_int(interim_enrollments)
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)

  interim_screenings_padded <- pad_int(interim_screenings)
  names(interim_screenings_padded) <- paste0("interim_screenings_", 1:max_interims)

  c(
    list(
      total_screened = total_screened,
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(interim_screenings_padded),
    flatten_joint(enroll_joint_list, max_interims, n_bins, "enrollments"),
    flatten_joint(screen_joint_list, max_interims, n_bins, "screenings")
  )
}

# =============================================================================
# DESIGN 2: AE (adaptive enrichment; no early stopping)
# =============================================================================
ae <- function(trial_data, max_screen, max_sample_size, info_rates,
               alpha, beta, type_alpha, type_beta, epsilon, d, g) {

  bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)

  total_screened <- 0L
  total_enrolled <- 0L

  enrolled_indices <- integer(0)
  enrolled_indices_current_stage <- integer(0)
  screened_indices_current_stage <- integer(0)

  model <- NULL
  test_statistics <- numeric(0)
  incremental_sample_sizes <- integer(0)
  interim_looks <- integer(0)

  conclusion <- NA_character_

  n_bins <- 20L
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)

  interim_enrollments <- integer(0)
  interim_screenings  <- integer(0)

  enroll_joint_list <- list()
  screen_joint_list <- list()

  start_idx <- 1L
  threshold <- 0
  screen_cohort <- ceiling(max_sample_size * info_rates[1])

  commit_stage <- function() {
    interim_enrollments <<- c(interim_enrollments, length(enrolled_indices_current_stage))
    interim_screenings  <<- c(interim_screenings,  length(screened_indices_current_stage))
    enroll_joint_list[[length(interim_enrollments)]] <<-
      safe_joint_counts_x1_x2bin(trial_data, enrolled_indices_current_stage, bin_breaks)
    screen_joint_list[[length(interim_screenings)]]  <<-
      safe_joint_counts_x1_x2bin(trial_data, screened_indices_current_stage, bin_breaks)
    enrolled_indices_current_stage <<- integer(0)
    screened_indices_current_stage <<- integer(0)
  }

  while (total_enrolled < max_sample_size && total_screened < max_screen) {

    next_milestone_index <- length(interim_looks) + 1L
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    remaining_to_milestone <- next_milestone - total_enrolled

    screen_count <- min(max_screen - total_screened, screen_cohort)
    end_idx <- min(start_idx + screen_count - 1L, nrow(trial_data))
    new_subjects <- trial_data[start_idx:end_idx, , drop = FALSE]

    total_screened <- total_screened + nrow(new_subjects)
    screened_indices_current_stage <- c(screened_indices_current_stage, new_subjects$subject_id)

    # enroll
    if (total_enrolled < info_rates[1] * max_sample_size) {
      num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
      subjects_to_enroll <- new_subjects[seq_len(num_to_enroll), , drop = FALSE]
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

    # interim look?
    if (total_enrolled >= next_milestone) {

      interim_looks <- c(interim_looks, total_enrolled)

      cumulative_data <- trial_data[enrolled_indices, , drop = FALSE]
      model <- update_model(cumulative_data)
      threshold <- adaptive_threshold(total_enrolled, max_sample_size, d, g)

      # record incremental test statistic at planned looks
      interim_data <- trial_data[enrolled_indices_current_stage, , drop = FALSE]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))

      # record stage totals + joint distributions
      commit_stage()
    }

    start_idx <- end_idx + 1L
    if (start_idx > nrow(trial_data)) break
  }

  # FIX: commit leftover buffered stage (if any)
  if (length(enrolled_indices_current_stage) > 0L || length(screened_indices_current_stage) > 0L) {
    commit_stage()
  }

  # final decision using last available stratified statistic (based on recorded increments)
  k <- length(test_statistics)
  if (k > 0) {
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
  } else {
    conclusion <- NA_character_
  }

  status <- if (total_enrolled == max_sample_size) {
    "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    "Screening limit reached"
  }

  max_interims <- length(info_rates)
  pad_int <- function(x) c(x, rep(NA, max_interims - length(x)))

  interim_enrollments_padded <- pad_int(interim_enrollments)
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)

  interim_screenings_padded <- pad_int(interim_screenings)
  names(interim_screenings_padded) <- paste0("interim_screenings_", 1:max_interims)

  c(
    list(
      total_screened = total_screened,
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(interim_screenings_padded),
    flatten_joint(enroll_joint_list, max_interims, n_bins, "enrollments"),
    flatten_joint(screen_joint_list, max_interims, n_bins, "screenings")
  )
}

# =============================================================================
# DESIGN 3: GSD (group sequential, no enrichment)
# =============================================================================
gsd <- function(trial_data, max_screen, max_sample_size, info_rates,
                alpha, beta, type_alpha, type_beta) {

  bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)

  total_screened <- 0L
  total_enrolled <- 0L

  enrolled_indices <- integer(0)
  enrolled_indices_current_stage <- integer(0)
  screened_indices_current_stage <- integer(0)

  test_statistics <- numeric(0)
  incremental_sample_sizes <- integer(0)
  interim_looks <- integer(0)

  early_stop <- FALSE
  conclusion <- NA_character_
  stop_loop <- FALSE

  n_bins <- 20L
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)

  interim_enrollments <- integer(0)
  interim_screenings  <- integer(0)

  enroll_joint_list <- list()
  screen_joint_list <- list()

  start_idx <- 1L
  screen_cohort <- ceiling(max_sample_size * info_rates[1])

  commit_stage <- function() {
    interim_enrollments <<- c(interim_enrollments, length(enrolled_indices_current_stage))
    interim_screenings  <<- c(interim_screenings,  length(screened_indices_current_stage))
    enroll_joint_list[[length(interim_enrollments)]] <<-
      safe_joint_counts_x1_x2bin(trial_data, enrolled_indices_current_stage, bin_breaks)
    screen_joint_list[[length(interim_screenings)]]  <<-
      safe_joint_counts_x1_x2bin(trial_data, screened_indices_current_stage, bin_breaks)
    enrolled_indices_current_stage <<- integer(0)
    screened_indices_current_stage <<- integer(0)
  }

  while (total_enrolled < max_sample_size && total_screened < max_screen && !stop_loop) {

    next_milestone_index <- length(interim_looks) + 1L
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    remaining_to_milestone <- next_milestone - total_enrolled

    screen_count <- min(max_screen - total_screened, screen_cohort)
    end_idx <- min(start_idx + screen_count - 1L, nrow(trial_data))
    new_subjects <- trial_data[start_idx:end_idx, , drop = FALSE]

    total_screened <- total_screened + nrow(new_subjects)
    screened_indices_current_stage <- c(screened_indices_current_stage, new_subjects$subject_id)

    # enroll everyone (no enrichment)
    num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
    subjects_to_enroll <- new_subjects[seq_len(num_to_enroll), , drop = FALSE]
    enrolled_indices <- c(enrolled_indices, subjects_to_enroll$subject_id)
    enrolled_indices_current_stage <- c(enrolled_indices_current_stage, subjects_to_enroll$subject_id)
    total_enrolled <- length(enrolled_indices)

    if (total_enrolled >= next_milestone) {

      interim_looks <- c(interim_looks, total_enrolled)

      interim_data <- trial_data[enrolled_indices_current_stage, , drop = FALSE]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))

      k <- length(test_statistics)
      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)

      commit_stage()

      if (total_enrolled < max_sample_size) {
        if (stratified_stat > bounds$efficacy_bound[k]) {
          early_stop <- TRUE; conclusion <- "Efficacy"
        } else if (stratified_stat < bounds$futility_bound[k]) {
          early_stop <- TRUE; conclusion <- "Futility"
        }
      } else {
        conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
      }

      if (early_stop) stop_loop <- TRUE
    }

    start_idx <- end_idx + 1L
    if (start_idx > nrow(trial_data)) break
  }

  # FIX: commit leftover buffered stage (if any)
  if (length(enrolled_indices_current_stage) > 0L || length(screened_indices_current_stage) > 0L) {
    commit_stage()
  }

  if (is.na(conclusion) && length(test_statistics) > 0) {
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
  }

  status <- if (total_enrolled == max_sample_size) {
    "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    "Screening limit reached"
  } else if (early_stop) {
    "Early stopping"
  }

  max_interims <- length(info_rates)
  pad_int <- function(x) c(x, rep(NA, max_interims - length(x)))

  interim_enrollments_padded <- pad_int(interim_enrollments)
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)

  interim_screenings_padded <- pad_int(interim_screenings)
  names(interim_screenings_padded) <- paste0("interim_screenings_", 1:max_interims)

  c(
    list(
      total_screened = total_screened,
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(interim_screenings_padded),
    flatten_joint(enroll_joint_list, max_interims, n_bins, "enrollments"),
    flatten_joint(screen_joint_list, max_interims, n_bins, "screenings")
  )
}

# =============================================================================
# DESIGN 4: GSE-F (GSE + feasibility-based futility)
# =============================================================================
gse_f <- function(trial_data, max_screen, max_sample_size, info_rates,
                  alpha, beta, type_alpha, type_beta, epsilon, d, g, delta = 0) {

  bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)

  total_screened <- 0L
  total_enrolled <- 0L

  enrolled_indices <- integer(0)
  enrolled_indices_current_stage <- integer(0)
  screened_indices_current_stage <- integer(0)

  model <- NULL
  test_statistics <- numeric(0)
  incremental_sample_sizes <- integer(0)
  interim_looks <- integer(0)

  early_stop <- FALSE
  conclusion <- NA_character_

  n_bins <- 20L
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)

  interim_enrollments <- integer(0)
  interim_screenings  <- integer(0)

  enroll_joint_list <- list()
  screen_joint_list <- list()

  start_idx <- 1L
  threshold <- 0
  screen_cohort <- ceiling(max_sample_size * info_rates[1])

  commit_stage <- function() {
    interim_enrollments <<- c(interim_enrollments, length(enrolled_indices_current_stage))
    interim_screenings  <<- c(interim_screenings,  length(screened_indices_current_stage))
    enroll_joint_list[[length(interim_enrollments)]] <<-
      safe_joint_counts_x1_x2bin(trial_data, enrolled_indices_current_stage, bin_breaks)
    screen_joint_list[[length(interim_screenings)]]  <<-
      safe_joint_counts_x1_x2bin(trial_data, screened_indices_current_stage, bin_breaks)
    enrolled_indices_current_stage <<- integer(0)
    screened_indices_current_stage <<- integer(0)
  }

  while (total_enrolled < max_sample_size && total_screened < max_screen) {

    next_milestone_index <- length(interim_looks) + 1L
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    remaining_to_milestone <- next_milestone - total_enrolled

    screen_count <- min(max_screen - total_screened, screen_cohort)
    end_idx <- min(start_idx + screen_count - 1L, nrow(trial_data))
    new_subjects <- trial_data[start_idx:end_idx, , drop = FALSE]

    # feasibility futility check (only after first look, when model exists)
    if (!is.null(model) && total_enrolled >= info_rates[1] * max_sample_size) {
      thr_now  <- adaptive_threshold(total_enrolled, max_sample_size, d, g)
      acc_rate <- estimate_acceptance_rate(model, epsilon, thr_now, new_subjects)
      proj_tot <- projected_total_screens_to_finish(total_enrolled, total_screened, max_sample_size, acc_rate)

      if (proj_tot > (max_screen + delta)) {
        early_stop <- TRUE
        conclusion <- "Futility"
        # FIX: commit partial stage so stage totals sum to total_*
        if (length(enrolled_indices_current_stage) > 0L || length(screened_indices_current_stage) > 0L) {
          commit_stage()
        }
        break
      }
    }

    total_screened <- total_screened + nrow(new_subjects)
    screened_indices_current_stage <- c(screened_indices_current_stage, new_subjects$subject_id)

    # enroll
    if (total_enrolled < info_rates[1] * max_sample_size) {
      num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
      subjects_to_enroll <- new_subjects[seq_len(num_to_enroll), , drop = FALSE]
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

      cumulative_data <- trial_data[enrolled_indices, , drop = FALSE]
      model <- update_model(cumulative_data)
      threshold <- adaptive_threshold(total_enrolled, max_sample_size, d, g)

      interim_data <- trial_data[enrolled_indices_current_stage, , drop = FALSE]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))

      k <- length(test_statistics)
      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)

      commit_stage()

      if (total_enrolled < max_sample_size) {
        if (stratified_stat > bounds$efficacy_bound[k]) {
          early_stop <- TRUE; conclusion <- "Efficacy"; break
        } else if (stratified_stat < bounds$futility_bound[k]) {
          early_stop <- TRUE; conclusion <- "Futility"; break
        }
      } else {
        conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
        break
      }
    }

    start_idx <- end_idx + 1L
    if (start_idx > nrow(trial_data)) break
  }

  # FIX: commit leftover buffered stage (if any), unless already committed at stop
  if (length(enrolled_indices_current_stage) > 0L || length(screened_indices_current_stage) > 0L) {
    commit_stage()
  }

  if (is.na(conclusion) && length(test_statistics) > 0) {
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
  }

  status <- if (total_enrolled == max_sample_size) {
    "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    "Screening limit reached"
  } else if (early_stop) {
    "Early stopping"
  }

  max_interims <- length(info_rates)
  pad_int <- function(x) c(x, rep(NA, max_interims - length(x)))

  interim_enrollments_padded <- pad_int(interim_enrollments)
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)

  interim_screenings_padded <- pad_int(interim_screenings)
  names(interim_screenings_padded) <- paste0("interim_screenings_", 1:max_interims)

  c(
    list(
      total_screened = total_screened,
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(interim_screenings_padded),
    flatten_joint(enroll_joint_list, max_interims, n_bins, "enrollments"),
    flatten_joint(screen_joint_list, max_interims, n_bins, "screenings")
  )
}

# =============================================================================
# Trial simulation wrapper for the mixed biomarker setting
# =============================================================================
trial_simulation <- function(sim_index,
                             design_type,
                             prev_rate_x1,
                             max_screen,
                             max_sample_size,
                             info_rates,
                             effect_size,
                             x1_neg_fraction = 0.5,
                             alpha, beta,
                             type_alpha, type_beta,
                             epsilon, d, g,
                             delta = 0,
                             scenario_index) {

  set.seed(sim_index)

  dat <- data.frame(
    subject_id = 1:max_screen,
    x1 = rbinom(max_screen, 1, prev_rate_x1),
    x2 = runif(max_screen, 0, 1),
    group = rbinom(max_screen, 1, 0.5)
  )

  beta0   <- 5
  beta_x1 <- 0.25
  beta_x2 <- 0.5

  sigma <- 0.1

  baseline  <- 3
  amplitude <- 2
  multiplier <- 1.2

  noise <- rnorm(max_screen, 0, sigma)

  # baseline mean under control
  if (scenario_index %in% c(1, 2, 3)) {
    mu_c <- beta0 + beta_x1 * dat$x1 + beta_x2 * dat$x2
  } else if (scenario_index == 4) {
    mu_c <- baseline + amplitude * (4 * (dat$x2 - 0.5)^2) + beta_x1 * dat$x1
  } else if (scenario_index == 5) {
    mu_c <- baseline + amplitude * (1 - 4 * (dat$x2 - 0.5)^2) + beta_x1 * dat$x1
  } else {
    stop("scenario_index must be 1..5")
  }

  # treatment effect delta(x1,x2): x1=0 gets a fraction of x1=1 effect
  x1_weight <- x1_neg_fraction + (1 - x1_neg_fraction) * dat$x1

  if (scenario_index == 1) {
    delta_x <- 0
  } else if (scenario_index == 2) {
    delta_x <- effect_size
  } else if (scenario_index == 3) {
    delta_x <- effect_size * dat$x2 * multiplier * x1_weight
  } else if (scenario_index == 4) {
    delta_x <- effect_size * 4 * (dat$x2 - 0.5)^2 * multiplier * x1_weight
  } else if (scenario_index == 5) {
    delta_x <- effect_size * (4 * dat$x2 * (1 - dat$x2)) * multiplier * x1_weight
  }

  mu_t <- mu_c + delta_x
  y <- ifelse(dat$group == 1, mu_t + noise, mu_c + noise)

  dat$y <- y
  dat$group <- factor(dat$group, levels = c(0, 1), labels = c("Control", "Treatment"))

  if (design_type == "GSE") {
    gse(dat, max_screen, max_sample_size, info_rates,
        alpha, beta, type_alpha, type_beta, epsilon, d, g)
  } else if (design_type == "AE") {
    ae(dat, max_screen, max_sample_size, info_rates,
       alpha, beta, type_alpha, type_beta, epsilon, d, g)
  } else if (design_type == "GSD") {
    gsd(dat, max_screen, max_sample_size, info_rates,
        alpha, beta, type_alpha, type_beta)
  } else if (design_type == "GSE-F") {
    gse_f(dat, max_screen, max_sample_size, info_rates,
          alpha, beta, type_alpha, type_beta, epsilon, d, g, delta = delta)
  } else {
    stop("design_type must be one of: 'GSE', 'AE', 'GSD', 'GSE-F'")
  }
}
