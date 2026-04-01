# Function to enroll subjects based on PBI
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

# Function to update regression model
update_model <- function(cumulative_data) {
  
  # Fit spline model
  updated_model <- stan_gamm4(
    formula = y ~ group + s(x) + s(x, by = group),
    data    = cumulative_data,
    prior = normal(0, 10),
    family  = gaussian(),      
    chains  = 4,               
    iter    = 2000,            
    cores   = 1                
  )
  
  return(updated_model)
}

# Function to calculate personalized benefit index (PBI) for a subject
# epsilon: minimal clinically significant difference
calculate_pbi <- function(model, new_subject_data, epsilon) {
  # Predict outcome for new subject based on updated model
  # posterior_treatment and posterior_control are matrices of dimension:
  # nrow = 4000 (4 chains and 2000 iterations, 1000 iterations for warmup, 1000 iterations for sampling)
  # ncol = number of subjects in new_subjects_data
  posterior_treatment <- posterior_predict(model, 
                                           newdata = transform(new_subject_data, group = factor("Treatment", levels = c("Control","Treatment"))))
  posterior_control <- posterior_predict(model,
                                         newdata = transform(new_subject_data, group = factor("Control", levels = c("Control","Treatment"))))
  
  # PBI
  pbi <- colMeans((posterior_treatment - posterior_control) > epsilon)
  
  return(pbi)
}

# Function to define adaptive threshold
adaptive_threshold <- function(n_k, N, d, g) {
  # N -- Total sample size of the trial
  # n_k -- Total number of subjects enrolled up to the current cohort
  # d,g -- Tuning parameters calibrated by simulation
  # if g is 0, the threshold is constant, Pocock type of threshold
  # if g is 1, the threshold is linearly increasing, O'Brien-Fleming type of threshold
  # d is a scaling factor that adjusts the base level of the threshold
  return(d*(n_k/N)^g)
}

# Function to calculate t-statistic at each interim stage
calculate_t_statistic <- function(interim_data) {
  treatment_y <- interim_data$y[interim_data$group=="Treatment"]
  control_y   <- interim_data$y[interim_data$group=="Control"]
  n_t <- length(treatment_y)
  n_c <- length(control_y)
  # if too few subjects, return a neutral t‐value (e.g. 0)
  if (n_t < 2 || n_c < 2) return(0)
  
  # otherwise do the usual
  mean_t  <- mean(treatment_y)
  mean_c  <- mean(control_y)
  var_t   <- var(treatment_y)
  var_c   <- var(control_y)
  den <- sqrt(var_t/n_t + var_c/n_c)
  if (!is.finite(den) || den == 0) return(0)
  (mean_t - mean_c) / den
}

# Function to calculate stratified test statistic (Simon & Simon, 2013)
calculate_stratified_statistic <- function(test_statistics, incremental_sample_sizes, k) {
  if (k <= 0) return(0)  # safe default when nothing has been analyzed yet
  idx <- seq_len(k)
  n_l <- incremental_sample_sizes[idx]
  n_l_sum <- sum(n_l)
  if (n_l_sum <= 0) return(0)  # guard against pathological zero sum
  weights <- sqrt(n_l / n_l_sum)
  sum(weights * test_statistics[idx])
}

# Function to calculate stopping boundaries for efficacy and futility
stopping_bounds <- function(alpha, beta, type_alpha, type_beta, informationRates) {
  gsd <- getDesignGroupSequential(typeOfDesign = type_alpha, alpha=alpha,
                                  typeBetaSpending = type_beta,
                                  bindingFutility = FALSE, informationRates = informationRates, 
                                  sided=1, beta=beta)
  
  return(list(efficacy_bound = gsd$criticalValues,
              futility_bound = c(gsd$futilityBounds, gsd$criticalValues[length(gsd$criticalValues)])))
}

# Function to map a vector of x values to bin indices
get_bin_indices <- function(x_values, bin_breaks) {
  # returns an integer vector, e.g. 1,2,...,n_bins
  as.numeric(cut(x_values, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE))
}

# Function to count subjects in each bin
count_subjects_in_bins <- function(x_values, bin_breaks) {
  # bin_indices will be integers in {1, 2, ..., n_bins}
  bin_indices <- get_bin_indices(x_values, bin_breaks)
  # Use tabulate to count how many values fall in each bin
  bin_counts <- tabulate(bin_indices, nbins = length(bin_breaks) - 1)
  return(bin_counts)
}

# Function to estimate empirical acceptance rate from current screening batch (new_subjects)
estimate_acceptance_rate <- function(model, epsilon, threshold, new_subjects) {
  if (is.null(model) || is.null(new_subjects) || nrow(new_subjects) == 0) return(NA_real_)
  pbi <- calculate_pbi(model, new_subjects, epsilon)
  mean(pbi > threshold)
}

# Function to project total number of screened subjects needed to finish enrollment
projected_total_screens_to_finish <- function(total_enrolled, total_screened,
                                              max_sample_size, acceptance_rate) {
  remaining_enroll <- max(0, max_sample_size - total_enrolled)
  if (is.na(acceptance_rate) || acceptance_rate <= 0) return(Inf)
  expected_additional_screens <- ceiling(remaining_enroll / acceptance_rate)
  total_screened + expected_additional_screens
}

# Function for adaptive enrichment with early stop
gse <- function(trial_data, max_screen, max_sample_size, info_rates, 
                          alpha, beta, type_alpha, type_beta, epsilon, d, g) {
  total_screened <- 0
  total_enrolled <- 0
  enrolled_indices <- c()
  enrolled_indices_current_stage <- c()
  screened_indices_current_stage <- c()
  model <- NULL
  test_statistics <- c()
  incremental_sample_sizes <- c()
  interim_looks <- c()
  early_stop <- FALSE
  conclusion <- NA
  stop_loop <- FALSE
  
  # Define bins for x
  n_bins <- 20
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)
  
  # Returns a full-length vector of bin counts; zeros if no indices
  safe_bin_counts <- function(id_vec) {
    if (length(id_vec) == 0L) {
      integer(length(bin_breaks) - 1L)
    } else {
      x_vals <- trial_data[id_vec, "x", drop = TRUE]
      count_subjects_in_bins(x_vals, bin_breaks)
    }
  }
  
  # Initialize interim result variables
  interim_enrollments <- c()
  interim_screenings <- c()
  enrollment_bin_counts_list <- list()
  screening_bin_counts_list  <- list()
  
  # Start index for screening
  start_idx <- 1
  threshold <- 0  # Initial threshold
  
  # Set screen cohort size
  screen_cohort <- max_sample_size*info_rates[1]
  
  while (total_enrolled < max_sample_size && total_screened < max_screen && !stop_loop) {
    # Calculate the next interim milestone
    next_milestone_index <- length(interim_looks) + 1
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    
    remaining_to_milestone <- next_milestone - total_enrolled
    
    # Determine the number of subjects to screen
    screen_count <- min(max_screen - total_screened, screen_cohort)
    
    # Get the new subjects to screen
    end_idx <- start_idx + screen_count - 1
    if (end_idx > nrow(trial_data)) {
      end_idx <- nrow(trial_data)
    }
    
    new_subjects <- trial_data[start_idx:end_idx, ]
    
    # Update total_screened and screened indices
    total_screened <- total_screened + nrow(new_subjects)
    screened_indices_current_stage <- c(screened_indices_current_stage, new_subjects$subject_id)
    
    # Enroll subjects
    if (total_enrolled < info_rates[1] * max_sample_size) {
      # Enroll all subjects before first interim milestone
      num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
      subjects_to_enroll <- new_subjects[1:num_to_enroll, ]
      
      enrolled_subjects <- subjects_to_enroll
      enrolled_indices <- c(enrolled_indices, subjects_to_enroll$subject_id)
      total_enrolled <- length(enrolled_indices)
      
      # Update enrolled_indices_current_stage (incremental)
      enrolled_indices_current_stage <- c(enrolled_indices_current_stage, subjects_to_enroll$subject_id)
    } else {
      # Enroll based on PBI and threshold, limited by remaining_to_milestone
      enrollment_result <- enroll_subjects(model, new_subjects, threshold, enrolled_indices, enrolled_indices_current_stage, remaining_to_milestone, epsilon)
      enrolled_subjects <- enrollment_result$enrolled_subjects
      enrolled_indices <- enrollment_result$enrolled_indices
      enrolled_indices_current_stage <- enrollment_result$enrolled_indices_current_stage
      total_enrolled <- length(enrolled_indices)
    }
    
    # Check if next interim milestone is reached
    if (total_enrolled >= next_milestone) {
      
      # Update interim_looks
      interim_looks <- c(interim_looks, total_enrolled)
      
      # Update the model using cumulative data
      cumulative_data <- trial_data[enrolled_indices, ]
      model <- update_model(cumulative_data)
      
      # Calculate adaptive threshold
      threshold <- adaptive_threshold(total_enrolled, max_sample_size, d, g)
      
      # Use incremental data for t-statistic calculation
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      
      # Calculate t-statistic based on incremental data
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      
      # Store sample size of incremental data
      n_k <- nrow(interim_data)
      incremental_sample_sizes <- c(incremental_sample_sizes, n_k)
      
      # Calculate stratified statistic using incremental sample sizes
      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, length(test_statistics))
      bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
      efficacy_bound <- bounds$efficacy_bound
      futility_bound <- bounds$futility_bound
      
      # Check for early stopping
      if (total_enrolled < max_sample_size) {
        if (stratified_stat > efficacy_bound[length(interim_looks)]) {
          early_stop <- TRUE
          conclusion <- "Efficacy"
        } else if (stratified_stat < futility_bound[length(interim_looks)]) {
          early_stop <- TRUE
          conclusion <- "Futility"
        }
      } else {
        # At final analysis, set conclusion without setting early_stop to TRUE
        if (stratified_stat > efficacy_bound[length(interim_looks)]) {
          conclusion <- "Efficacy"
        } else {
          conclusion <- "Futility"
        }
      }
      
      # Update interim result variables
      interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
      interim_screenings  <- c(interim_screenings,  length(screened_indices_current_stage))
      
      enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
      screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)
      
      # Reset current stage indices
      screened_indices_current_stage <- c()
      enrolled_indices_current_stage <- c()
      
      if (early_stop) {
        stop_loop <- TRUE
      }
    }
    
    start_idx <- end_idx + 1
    if (start_idx > nrow(trial_data)) {
      break
    }
    
  }
  
  # Final‐look analysis if we never early‐stopped and haven’t done the last info‐rate
  max_interims <- length(info_rates)
  if (!early_stop && total_enrolled > 0 && length(interim_looks) < max_interims) {
    # record this “last” look
    interim_looks <- c(interim_looks, total_enrolled)
    if (length(enrolled_indices_current_stage) > 0) {
      # do the usual incremental‐data analysis
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
    }
    
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
    # capture counts for this final partial look
    interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
    interim_screenings  <- c(interim_screenings,  length(screened_indices_current_stage))
    enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
    screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)
    # clear the “current_stage” buffers
    enrolled_indices_current_stage  <- integer(0)
    screened_indices_current_stage  <- integer(0)
  }
  
  # Define trial status
  if (total_enrolled == max_sample_size) {
    status <- "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    status <- "Screening limit reached"
  } else if (early_stop) {
    status <- "Early stopping"
  }
  
  # --- GUARANTEE a final analysis & conclusion if we exited the loop without one ---
  # (This covers the "Full enrollment" path where the last look wasn't finalized inside the loop.)
  if (is.na(conclusion) && total_enrolled > 0) {
    max_interims <- length(info_rates)
    have_final_look_at_N <- length(interim_looks) > 0 && tail(interim_looks, 1) >= max_sample_size
    
    # If we don't already have a recorded final look at N, add one using the current-stage buffer.
    if (!have_final_look_at_N) {
      # Create a final look at the current total enrolled
      interim_looks <- c(interim_looks, total_enrolled)
      
      if (length(enrolled_indices_current_stage) > 0) {
        # Use the incremental batch since the prior look
        interim_data <- trial_data[enrolled_indices_current_stage, ]
        t_stat <- calculate_t_statistic(interim_data)
        test_statistics <- c(test_statistics, t_stat)
        incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
        
        # Record per-bin counts for this final (partial) stage
        enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
        screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)
        
        # Clear stage buffers now that they’ve been committed
        enrolled_indices_current_stage  <- integer(0)
        screened_indices_current_stage  <- integer(0)
        
      } else if (length(test_statistics) == 0) {
        # Edge case: no looks recorded at all — analyze the full enrolled set as one look
        interim_data <- trial_data[enrolled_indices, ]
        t_stat <- calculate_t_statistic(interim_data)
        test_statistics <- c(t_stat)
        incremental_sample_sizes <- c(nrow(interim_data))
        
        # Bin counts based on all enrolled to keep lengths consistent
        enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices)
        # We don't have a well-defined "incremental screening" set here; record zeros
        screening_bin_counts_list [[length(interim_looks)]] <- integer(length(bin_breaks) - 1)
      }
    }
    
    # Now make the final decision using the last available look
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
  }
  
  # Pad interim results to length of info_rates
  max_interims <- length(info_rates)
  
  # Pad interim_enrollments and interim_screenings
  interim_enrollments_padded <- c(interim_enrollments, rep(NA, max_interims - length(interim_enrollments)))
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)
  interim_screenings_padded <- c(interim_screenings, rep(NA, max_interims - length(interim_screenings)))
  names(interim_screenings_padded) <- paste0("interim_screenings_", 1:max_interims)
  
  # Define bin names
  enroll_bin_count_names <- paste0("enrollments_bin", 1:n_bins, "_interim_")
  screen_bin_count_names <- paste0("screenings_bin", 1:n_bins, "_interim_")
  
  # Flattening enrollment counts
  flattened_enroll_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(enrollment_bin_counts_list)) {
      # have actual data
      for (b in seq_len(n_bins)) {
        flattened_enroll_counts[[ paste0(enroll_bin_count_names[b], k) ]] <- 
          enrollment_bin_counts_list[[k]][b]
      }
    } else {
      # pad with NA
      for (b in seq_len(n_bins)) {
        flattened_enroll_counts[[ paste0(enroll_bin_count_names[b], k) ]] <- NA
      }
    }
  }
  
  # Flattening screening counts
  flattened_screen_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(screening_bin_counts_list)) {
      # have actual data
      for (b in seq_len(n_bins)) {
        flattened_screen_counts[[ paste0(screen_bin_count_names[b], k) ]] <- 
          screening_bin_counts_list[[k]][b]
      }
    } else {
      # pad with NA
      for (b in seq_len(n_bins)) {
        flattened_screen_counts[[ paste0(screen_bin_count_names[b], k) ]] <- NA
      }
    }
  }
  
  # Return results
  return(c(
    list(
      total_screened = total_screened,
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(interim_screenings_padded),
    as.list(flattened_enroll_counts),
    as.list(flattened_screen_counts)
  ))
}

# Function for adaptive enrichment with early stop (additional futility check)
gse_plus <- function(trial_data, max_screen, max_sample_size, info_rates, 
                               alpha, beta, type_alpha, type_beta, epsilon, d, g, delta) {
  total_screened <- 0
  total_enrolled <- 0
  enrolled_indices <- c()
  enrolled_indices_current_stage <- c()
  screened_indices_current_stage <- c()
  model <- NULL
  test_statistics <- c()
  incremental_sample_sizes <- c()
  interim_looks <- c()
  early_stop <- FALSE
  conclusion <- NA
  
  # bins
  n_bins <- 20
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)
  
  # Returns a full-length vector of bin counts; zeros if no indices
  safe_bin_counts <- function(id_vec) {
    if (length(id_vec) == 0L) {
      integer(length(bin_breaks) - 1L)
    } else {
      x_vals <- trial_data[id_vec, "x", drop = TRUE]
      count_subjects_in_bins(x_vals, bin_breaks)
    }
  }
  
  # interim trackers
  interim_enrollments <- c()
  interim_screenings  <- c()
  enrollment_bin_counts_list <- list()
  screening_bin_counts_list  <- list()
  
  # small local helper: commit current (partial) batch as one interim snapshot
  commit_current_stage <- function() {
    # tag a new "look" at the current total enrolled
    interim_looks <<- c(interim_looks, total_enrolled)
    
    # sizes
    k_enr <- length(enrolled_indices_current_stage)
    k_scr <- length(screened_indices_current_stage)
    interim_enrollments <<- c(interim_enrollments, k_enr)
    interim_screenings  <<- c(interim_screenings,  k_scr)
    
    # per-bin counts (match exactly what's in the current-stage buffers)
    if (k_enr > 0) {
      enr_data <- trial_data[enrolled_indices_current_stage, ]
      enrollment_bin_counts_list[[length(interim_looks)]] <<- count_subjects_in_bins(enr_data$x, bin_breaks)
    } else {
      enrollment_bin_counts_list[[length(interim_looks)]] <<- integer(length(bin_breaks) - 1)
    }
    if (k_scr > 0) {
      scr_data <- trial_data[screened_indices_current_stage, ]
      screening_bin_counts_list[[length(interim_looks)]] <<- count_subjects_in_bins(scr_data$x, bin_breaks)
    } else {
      screening_bin_counts_list[[length(interim_looks)]] <<- integer(length(bin_breaks) - 1)
    }
  }
  
  # screening loop
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
        # record this partial stage so counts remain consistent
        early_stop <- TRUE
        conclusion <- "Futility"
        # add what we are ABOUT to screen/enroll in this batch? No—because we haven't
        # added them to the running totals yet. We commit exactly what is in the
        # current-stage buffers so far.
        if (length(enrolled_indices_current_stage) > 0 || length(screened_indices_current_stage) > 0) {
          commit_current_stage()
          # do NOT clear buffers; we are exiting anyway
        }
        break
      }
    }
    
    # update screened totals/buffer
    total_screened <- total_screened + nrow(new_subjects)
    screened_indices_current_stage <- c(screened_indices_current_stage, new_subjects$subject_id)
    
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
    
    # hit an info-rate? do a standard look
    if (total_enrolled >= next_milestone) {
      # update model and threshold
      cumulative_data <- trial_data[enrolled_indices, ]
      model <- update_model(cumulative_data)
      threshold <- adaptive_threshold(total_enrolled, max_sample_size, d, g)
      
      # compute t for incremental data (safe for very small n)
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
      
      # group-seq decision
      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, length(test_statistics))
      bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
      
      # record this interim (counts + bins) BEFORE clearing buffers
      interim_looks <- c(interim_looks, total_enrolled)
      interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
      interim_screenings  <- c(interim_screenings,  length(screened_indices_current_stage))
      enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
      screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)

      # clear stage buffers
      screened_indices_current_stage <- c()
      enrolled_indices_current_stage <- c()
      
      # check stopping
      if (total_enrolled < max_sample_size) {
        if (stratified_stat > bounds$efficacy_bound[length(interim_looks)]) {
          early_stop <- TRUE; conclusion <- "Efficacy"; break
        } else if (stratified_stat < bounds$futility_bound[length(interim_looks)]) {
          early_stop <- TRUE; conclusion <- "Futility"; break
        }
      } else {
        conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[length(interim_looks)],
                             "Efficacy","Futility")
        break
      }
    }
    
    start_idx <- end_idx + 1
    if (start_idx > nrow(trial_data)) break
  }
  
  # Final‐look analysis if we never early‐stopped and haven’t done the last info‐rate
  max_interims <- length(info_rates)
  if (!early_stop && total_enrolled > 0 && length(interim_looks) < max_interims) {
    # record this “last” look
    interim_looks <- c(interim_looks, total_enrolled)
    if (length(enrolled_indices_current_stage) > 0) {
      # do the usual incremental‐data analysis
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
    }
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
    # capture counts for this final partial look
    interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
    interim_screenings  <- c(interim_screenings,  length(screened_indices_current_stage))
    enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
    screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)
    # clear the “current_stage” buffers
    enrolled_indices_current_stage  <- integer(0)
    screened_indices_current_stage  <- integer(0)
  }
  
  # Define trial status
  if (total_enrolled == max_sample_size) {
    status <- "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    status <- "Screening limit reached"
  } else if (early_stop) {
    status <- "Early stopping"
  }
  
  # --- GUARANTEE a final analysis & conclusion if we exited the loop without one ---
  # (This covers the "Full enrollment" path where the last look wasn't finalized inside the loop.)
  if (is.na(conclusion) && total_enrolled > 0) {
    max_interims <- length(info_rates)
    have_final_look_at_N <- length(interim_looks) > 0 && tail(interim_looks, 1) >= max_sample_size
    
    # If we don't already have a recorded final look at N, add one using the current-stage buffer.
    if (!have_final_look_at_N) {
      # Create a final look at the current total enrolled
      interim_looks <- c(interim_looks, total_enrolled)
      
      if (length(enrolled_indices_current_stage) > 0) {
        # Use the incremental batch since the prior look
        interim_data <- trial_data[enrolled_indices_current_stage, ]
        t_stat <- calculate_t_statistic(interim_data)
        test_statistics <- c(test_statistics, t_stat)
        incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
        
        # Record per-bin counts for this final (partial) stage
        enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
        screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)
        
        # Clear stage buffers now that they’ve been committed
        enrolled_indices_current_stage  <- integer(0)
        screened_indices_current_stage  <- integer(0)
        
      } else if (length(test_statistics) == 0) {
        # Edge case: no looks recorded at all — analyze the full enrolled set as one look
        interim_data <- trial_data[enrolled_indices, ]
        t_stat <- calculate_t_statistic(interim_data)
        test_statistics <- c(t_stat)
        incremental_sample_sizes <- c(nrow(interim_data))
        
        # Bin counts based on all enrolled to keep lengths consistent
        enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices)
        # We don't have a well-defined "incremental screening" set here; record zeros
        screening_bin_counts_list [[length(interim_looks)]] <- integer(length(bin_breaks) - 1)
      }
    }
    
    # Now make the final decision using the last available look
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
  }
  
  # pad to |info_rates|
  max_interims <- length(info_rates)
  interim_enrollments_padded <- c(interim_enrollments, rep(NA, max_interims - length(interim_enrollments)))
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)
  interim_screenings_padded  <- c(interim_screenings,  rep(NA, max_interims - length(interim_screenings)))
  names(interim_screenings_padded)  <- paste0("interim_screenings_", 1:max_interims)
  
  # flatten bin counts
  enroll_bin_count_names <- paste0("enrollments_bin", 1:n_bins, "_interim_")
  screen_bin_count_names <- paste0("screenings_bin",  1:n_bins, "_interim_")
  
  flattened_enroll_counts <- list()
  flattened_screen_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(enrollment_bin_counts_list)) {
      for (b in seq_len(n_bins)) {
        flattened_enroll_counts[[paste0(enroll_bin_count_names[b], k)]] <- enrollment_bin_counts_list[[k]][b]
      }
    } else {
      for (b in seq_len(n_bins)) {
        flattened_enroll_counts[[paste0(enroll_bin_count_names[b], k)]] <- NA
      }
    }
    if (k <= length(screening_bin_counts_list)) {
      for (b in seq_len(n_bins)) {
        flattened_screen_counts[[paste0(screen_bin_count_names[b], k)]] <- screening_bin_counts_list[[k]][b]
      }
    } else {
      for (b in seq_len(n_bins)) {
        flattened_screen_counts[[paste0(screen_bin_count_names[b], k)]] <- NA
      }
    }
  }
  
  c(
    list(
      total_screened = total_screened,
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(interim_screenings_padded),
    as.list(flattened_enroll_counts),
    as.list(flattened_screen_counts)
  )
}

# Function for adaptive enrichment without early stop
ae <- function(trial_data, max_screen, max_sample_size, info_rates, 
                             alpha, beta, type_alpha, type_beta, epsilon, d, g) {
  total_screened <- 0
  total_enrolled <- 0
  enrolled_indices <- c()
  enrolled_indices_current_stage <- c()
  screened_indices_current_stage <- c()
  model <- NULL
  test_statistics <- c()
  incremental_sample_sizes <- c()
  interim_looks <- c()
  early_stop <- FALSE
  conclusion <- NA
  #stop_loop <- FALSE
  
  # Define bins for x
  n_bins <- 20
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)
  
  # Returns a full-length vector of bin counts; zeros if no indices
  safe_bin_counts <- function(id_vec) {
    if (length(id_vec) == 0L) {
      integer(length(bin_breaks) - 1L)
    } else {
      x_vals <- trial_data[id_vec, "x", drop = TRUE]
      count_subjects_in_bins(x_vals, bin_breaks)
    }
  }
  
  # Initialize interim result variables
  interim_enrollments <- c()
  interim_screenings <- c()
  enrollment_bin_counts_list <- list()
  screening_bin_counts_list  <- list()
  
  # Start index for screening
  start_idx <- 1
  threshold <- 0  # Initial threshold
  
  # Set screen cohort size
  screen_cohort <- max_sample_size*info_rates[1]
  
  while (total_enrolled < max_sample_size && total_screened < max_screen) {
    # Calculate the next interim milestone
    next_milestone_index <- length(interim_looks) + 1
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    
    remaining_to_milestone <- next_milestone - total_enrolled
    
    # Determine the number of subjects to screen
    screen_count <- min(max_screen - total_screened, screen_cohort)
    
    # Get the new subjects to screen
    end_idx <- start_idx + screen_count - 1
    if (end_idx > nrow(trial_data)) {
      end_idx <- nrow(trial_data)
    }
    
    new_subjects <- trial_data[start_idx:end_idx, ]
    
    # Update total_screened and screened indices
    total_screened <- total_screened + nrow(new_subjects)
    screened_indices_current_stage <- c(screened_indices_current_stage, new_subjects$subject_id)
    
    # Enroll subjects
    if (total_enrolled < info_rates[1] * max_sample_size) {
      # Enroll all subjects before first interim milestone
      num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
      subjects_to_enroll <- new_subjects[1:num_to_enroll, ]
      
      enrolled_subjects <- subjects_to_enroll
      enrolled_indices <- c(enrolled_indices, subjects_to_enroll$subject_id)
      total_enrolled <- length(enrolled_indices)
      
      # Update enrolled_indices_current_stage (incremental)
      enrolled_indices_current_stage <- c(enrolled_indices_current_stage, subjects_to_enroll$subject_id)
    } else {
      # Enroll based on PBI and threshold, limited by remaining_to_milestone
      enrollment_result <- enroll_subjects(model, new_subjects, threshold, enrolled_indices, enrolled_indices_current_stage, remaining_to_milestone, epsilon)
      enrolled_subjects <- enrollment_result$enrolled_subjects
      enrolled_indices <- enrollment_result$enrolled_indices
      enrolled_indices_current_stage <- enrollment_result$enrolled_indices_current_stage
      total_enrolled <- length(enrolled_indices)
    }
    
    # Check if next interim milestone is reached
    if (total_enrolled >= next_milestone) {
      
      # Update interim_looks
      interim_looks <- c(interim_looks, total_enrolled)
      
      # Update the model using cumulative data
      cumulative_data <- trial_data[enrolled_indices, ]
      model <- update_model(cumulative_data)
      
      # Calculate adaptive threshold
      threshold <- adaptive_threshold(total_enrolled, max_sample_size, d, g)
      
      # Use incremental data for t-statistic calculation
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      
      # Calculate t-statistic based on incremental data
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      
      # Store sample size of incremental data
      n_k <- nrow(interim_data)
      incremental_sample_sizes <- c(incremental_sample_sizes, n_k)
      
      # Update interim result variables
      interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
      interim_screenings  <- c(interim_screenings,  length(screened_indices_current_stage))
      enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
      screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)
      
      # Reset current stage indices
      screened_indices_current_stage <- c()
      enrolled_indices_current_stage <- c()
      
    }
    
    start_idx <- end_idx + 1
    if (start_idx > nrow(trial_data)) {
      break
    }
  }
  
  # Final‐look analysis if we never did the last info‐rate
  max_interims <- length(info_rates)
  if (total_enrolled > 0 && length(interim_looks) < max_interims) {
    interim_looks <- c(interim_looks, total_enrolled)
    if (length(enrolled_indices_current_stage) > 0) {
      # do the usual incremental‐data analysis
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
    }
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
    # capture counts
    interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
    interim_screenings  <- c(interim_screenings,  length(screened_indices_current_stage))
    enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
    screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)
    
    enrolled_indices_current_stage <- integer(0)
    screened_indices_current_stage <- integer(0)
  }
  
  # Define trial status
  if (total_enrolled == max_sample_size) {
    status <- "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    status <- "Screening limit reached"
  } else if (early_stop) {
    status <- "Early stopping"
  }
  
  # --- GUARANTEE a final analysis & conclusion if we exited the loop without one ---
  # (This covers the "Full enrollment" path where the last look wasn't finalized inside the loop.)
  if (is.na(conclusion) && total_enrolled > 0) {
    max_interims <- length(info_rates)
    have_final_look_at_N <- length(interim_looks) > 0 && tail(interim_looks, 1) >= max_sample_size
    
    # If we don't already have a recorded final look at N, add one using the current-stage buffer.
    if (!have_final_look_at_N) {
      # Create a final look at the current total enrolled
      interim_looks <- c(interim_looks, total_enrolled)
      
      if (length(enrolled_indices_current_stage) > 0) {
        # Use the incremental batch since the prior look
        interim_data <- trial_data[enrolled_indices_current_stage, ]
        t_stat <- calculate_t_statistic(interim_data)
        test_statistics <- c(test_statistics, t_stat)
        incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
        
        # Record per-bin counts for this final (partial) stage
        enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
        screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)
        
        # Clear stage buffers now that they’ve been committed
        enrolled_indices_current_stage  <- integer(0)
        screened_indices_current_stage  <- integer(0)
        
      } else if (length(test_statistics) == 0) {
        # Edge case: no looks recorded at all — analyze the full enrolled set as one look
        interim_data <- trial_data[enrolled_indices, ]
        t_stat <- calculate_t_statistic(interim_data)
        test_statistics <- c(t_stat)
        incremental_sample_sizes <- c(nrow(interim_data))
        
        # Bin counts based on all enrolled to keep lengths consistent
        enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices)
        # We don't have a well-defined "incremental screening" set here; record zeros
        screening_bin_counts_list [[length(interim_looks)]] <- integer(length(bin_breaks) - 1)
      }
    }
    
    # Now make the final decision using the last available look
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
  }
  
  # Pad interim results to length of info_rates
  max_interims <- length(info_rates)
  
  # Pad interim_enrollments and interim_screenings
  interim_enrollments_padded <- c(interim_enrollments, rep(NA, max_interims - length(interim_enrollments)))
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)
  interim_screenings_padded <- c(interim_screenings, rep(NA, max_interims - length(interim_screenings)))
  names(interim_screenings_padded) <- paste0("interim_screenings_", 1:max_interims)
  
  # Define bin names
  enroll_bin_count_names <- paste0("enrollments_bin", 1:n_bins, "_interim_")
  screen_bin_count_names <- paste0("screenings_bin", 1:n_bins, "_interim_")
  
  # Flattening enrollment counts
  flattened_enroll_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(enrollment_bin_counts_list)) {
      # have actual data
      for (b in seq_len(n_bins)) {
        flattened_enroll_counts[[ paste0(enroll_bin_count_names[b], k) ]] <- 
          enrollment_bin_counts_list[[k]][b]
      }
    } else {
      # pad with NA
      for (b in seq_len(n_bins)) {
        flattened_enroll_counts[[ paste0(enroll_bin_count_names[b], k) ]] <- NA
      }
    }
  }
  
  # Flattening screening counts
  flattened_screen_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(screening_bin_counts_list)) {
      # have actual data
      for (b in seq_len(n_bins)) {
        flattened_screen_counts[[ paste0(screen_bin_count_names[b], k) ]] <- 
          screening_bin_counts_list[[k]][b]
      }
    } else {
      # pad with NA
      for (b in seq_len(n_bins)) {
        flattened_screen_counts[[ paste0(screen_bin_count_names[b], k) ]] <- NA
      }
    }
  }
  
  # Return results
  return(c(
    list(
      total_screened = total_screened,
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(interim_screenings_padded),
    as.list(flattened_enroll_counts),
    as.list(flattened_screen_counts)
  ))
}

# Function for group sequential design without enrichment
gsd <- function(trial_data, max_screen, max_sample_size, info_rates, 
                          alpha, beta, type_alpha, type_beta) {
  total_screened <- 0
  total_enrolled <- 0
  enrolled_indices <- c()
  enrolled_indices_current_stage <- c()
  screened_indices_current_stage <- c()
  model <- NULL
  test_statistics <- c()
  incremental_sample_sizes <- c()
  interim_looks <- c()
  early_stop <- FALSE
  conclusion <- NA
  stop_loop <- FALSE
  
  # Define bins for x
  n_bins <- 20
  bin_breaks <- seq(0, 1, length.out = n_bins + 1)
  
  # Returns a full-length vector of bin counts; zeros if no indices
  safe_bin_counts <- function(id_vec) {
    if (length(id_vec) == 0L) {
      integer(length(bin_breaks) - 1L)
    } else {
      x_vals <- trial_data[id_vec, "x", drop = TRUE]
      count_subjects_in_bins(x_vals, bin_breaks)
    }
  }
  
  # Initialize interim result variables
  interim_enrollments <- c()
  interim_screenings <- c()
  enrollment_bin_counts_list <- list()
  screening_bin_counts_list  <- list()
  
  # Start index for screening
  start_idx <- 1
  threshold <- 0  # Initial threshold
  
  # Set screen cohort size
  screen_cohort <- max_sample_size*info_rates[1]
  
  while (total_enrolled < max_sample_size && total_screened < max_screen && !stop_loop) {
    # Calculate the next interim milestone
    next_milestone_index <- length(interim_looks) + 1
    next_milestone <- min(info_rates[next_milestone_index] * max_sample_size, max_sample_size)
    
    remaining_to_milestone <- next_milestone - total_enrolled
    
    # Determine the number of subjects to screen
    screen_count <- min(max_screen - total_screened, screen_cohort)
    
    # Get the new subjects to screen
    end_idx <- start_idx + screen_count - 1
    if (end_idx > nrow(trial_data)) {
      end_idx <- nrow(trial_data)
    }
    
    new_subjects <- trial_data[start_idx:end_idx, ]
    
    # Update total_screened and screened indices
    total_screened <- total_screened + nrow(new_subjects)
    screened_indices_current_stage <- c(screened_indices_current_stage, new_subjects$subject_id)
    
    # Enroll subjects without enrichment
    num_to_enroll <- min(nrow(new_subjects), remaining_to_milestone)
    subjects_to_enroll <- new_subjects[1:num_to_enroll, ]
    
    enrolled_subjects <- subjects_to_enroll
    enrolled_indices <- c(enrolled_indices, subjects_to_enroll$subject_id)
    total_enrolled <- length(enrolled_indices)
    
    # Update enrolled_indices_current_stage (incremental)
    enrolled_indices_current_stage <- c(enrolled_indices_current_stage, subjects_to_enroll$subject_id)
    
    # Check if next interim milestone is reached
    if (total_enrolled >= next_milestone) {
      
      # Update interim_looks
      interim_looks <- c(interim_looks, total_enrolled)
      
      # Update the model using cumulative data
      cumulative_data <- trial_data[enrolled_indices, ]
      
      # Use incremental data for t-statistic calculation
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      
      # Calculate t-statistic based on incremental data
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      
      # Store sample size of incremental data
      n_k <- nrow(interim_data)
      incremental_sample_sizes <- c(incremental_sample_sizes, n_k)
      
      # Calculate stratified statistic using incremental sample sizes
      stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, length(test_statistics))
      bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
      efficacy_bound <- bounds$efficacy_bound
      futility_bound <- bounds$futility_bound
      
      # Check for early stopping
      if (total_enrolled < max_sample_size) {
        if (stratified_stat > efficacy_bound[length(interim_looks)]) {
          early_stop <- TRUE
          conclusion <- "Efficacy"
        } else if (stratified_stat < futility_bound[length(interim_looks)]) {
          early_stop <- TRUE
          conclusion <- "Futility"
        }
      } else {
        # At final analysis, set conclusion without setting early_stop to TRUE
        if (stratified_stat > efficacy_bound[length(interim_looks)]) {
          conclusion <- "Efficacy"
        } else {
          conclusion <- "Futility"
        }
      }
      
      # Update interim result variables
      interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
      interim_screenings  <- c(interim_screenings,  length(screened_indices_current_stage))
      enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
      screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)

      # Reset current stage indices
      screened_indices_current_stage <- c()
      enrolled_indices_current_stage <- c()
      
      if (early_stop) {
        stop_loop <- TRUE
      }
    }
    
    start_idx <- end_idx + 1
    if (start_idx > nrow(trial_data)) {
      break
    }
  }
  
  # Final‐look analysis if we never did the last info‐rate
  max_interims <- length(info_rates)
  if (total_enrolled > 0 && length(interim_looks) < max_interims) {
    interim_looks <- c(interim_looks, total_enrolled)
    if (length(enrolled_indices_current_stage) > 0) {
      # do the usual incremental‐data analysis
      interim_data <- trial_data[enrolled_indices_current_stage, ]
      t_stat <- calculate_t_statistic(interim_data)
      test_statistics <- c(test_statistics, t_stat)
      incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
    }
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
    # capture counts
    interim_enrollments <- c(interim_enrollments, length(enrolled_indices_current_stage))
    interim_screenings  <- c(interim_screenings,  length(screened_indices_current_stage))
    enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
    screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)
    
    enrolled_indices_current_stage <- integer(0)
    screened_indices_current_stage <- integer(0)
  }
  
  # Define trial status
  if (total_enrolled == max_sample_size) {
    status <- "Full enrollment"
  } else if (total_enrolled < max_sample_size && total_screened == max_screen) {
    status <- "Screening limit reached"
  } else if (early_stop) {
    status <- "Early stopping"
  }
  
  # --- GUARANTEE a final analysis & conclusion if we exited the loop without one ---
  # (This covers the "Full enrollment" path where the last look wasn't finalized inside the loop.)
  if (is.na(conclusion) && total_enrolled > 0) {
    max_interims <- length(info_rates)
    have_final_look_at_N <- length(interim_looks) > 0 && tail(interim_looks, 1) >= max_sample_size
    
    # If we don't already have a recorded final look at N, add one using the current-stage buffer.
    if (!have_final_look_at_N) {
      # Create a final look at the current total enrolled
      interim_looks <- c(interim_looks, total_enrolled)
      
      if (length(enrolled_indices_current_stage) > 0) {
        # Use the incremental batch since the prior look
        interim_data <- trial_data[enrolled_indices_current_stage, ]
        t_stat <- calculate_t_statistic(interim_data)
        test_statistics <- c(test_statistics, t_stat)
        incremental_sample_sizes <- c(incremental_sample_sizes, nrow(interim_data))
        
        # Record per-bin counts for this final (partial) stage
        enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices_current_stage)
        screening_bin_counts_list [[length(interim_looks)]] <- safe_bin_counts(screened_indices_current_stage)
        
        # Clear stage buffers now that they’ve been committed
        enrolled_indices_current_stage  <- integer(0)
        screened_indices_current_stage  <- integer(0)
        
      } else if (length(test_statistics) == 0) {
        # Edge case: no looks recorded at all — analyze the full enrolled set as one look
        interim_data <- trial_data[enrolled_indices, ]
        t_stat <- calculate_t_statistic(interim_data)
        test_statistics <- c(t_stat)
        incremental_sample_sizes <- c(nrow(interim_data))
        
        # Bin counts based on all enrolled to keep lengths consistent
        enrollment_bin_counts_list[[length(interim_looks)]] <- safe_bin_counts(enrolled_indices)
        # We don't have a well-defined "incremental screening" set here; record zeros
        screening_bin_counts_list [[length(interim_looks)]] <- integer(length(bin_breaks) - 1)
      }
    }
    
    # Now make the final decision using the last available look
    bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
    k <- length(test_statistics)
    stratified_stat <- calculate_stratified_statistic(test_statistics, incremental_sample_sizes, k)
    conclusion <- ifelse(stratified_stat > bounds$efficacy_bound[k], "Efficacy", "Futility")
  }
  
  # Pad interim results to length of info_rates
  max_interims <- length(info_rates)
  
  # Pad interim_enrollments and interim_screenings
  interim_enrollments_padded <- c(interim_enrollments, rep(NA, max_interims - length(interim_enrollments)))
  names(interim_enrollments_padded) <- paste0("interim_enrollments_", 1:max_interims)
  interim_screenings_padded <- c(interim_screenings, rep(NA, max_interims - length(interim_screenings)))
  names(interim_screenings_padded) <- paste0("interim_screenings_", 1:max_interims)
  
  # Define bin names
  enroll_bin_count_names <- paste0("enrollments_bin", 1:n_bins, "_interim_")
  screen_bin_count_names <- paste0("screenings_bin", 1:n_bins, "_interim_")
  
  # Flattening enrollment counts
  flattened_enroll_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(enrollment_bin_counts_list)) {
      # have actual data
      for (b in seq_len(n_bins)) {
        flattened_enroll_counts[[ paste0(enroll_bin_count_names[b], k) ]] <- 
          enrollment_bin_counts_list[[k]][b]
      }
    } else {
      # pad with NA
      for (b in seq_len(n_bins)) {
        flattened_enroll_counts[[ paste0(enroll_bin_count_names[b], k) ]] <- NA
      }
    }
  }
  
  # Flattening screening counts
  flattened_screen_counts <- list()
  for (k in seq_len(max_interims)) {
    if (k <= length(screening_bin_counts_list)) {
      # have actual data
      for (b in seq_len(n_bins)) {
        flattened_screen_counts[[ paste0(screen_bin_count_names[b], k) ]] <- 
          screening_bin_counts_list[[k]][b]
      }
    } else {
      # pad with NA
      for (b in seq_len(n_bins)) {
        flattened_screen_counts[[ paste0(screen_bin_count_names[b], k) ]] <- NA
      }
    }
  }
  
  # Return results
  return(c(
    list(
      total_screened = total_screened,
      total_enrolled = total_enrolled,
      status = status,
      conclusion = conclusion
    ),
    as.list(interim_enrollments_padded),
    as.list(interim_screenings_padded),
    as.list(flattened_enroll_counts),
    as.list(flattened_screen_counts)
  ))
}


# Function to simulate trial design

# Parameters:
# design_type: type of trial design ("AE with early stop", "AE without early stop", "GSD without enrichment")
# max_screen: maximum number of subjects for screening
# max_sample_size: number of subjects planned to enroll
# info_rates: vector of information rates for interim analyses
# effect_size: effect size for treatment effect
# alpha: type-I error rate
# beta: type-II error rate
# type_alpha: type of alpha spending function (Pocock or O'Brien-Fleming)
# type_beta: type of beta spending function (Pocock or O'Brien-Fleming)
# epsilon: minimal clinically significant difference
# d: parameter for adaptive threshold
# g: parameter for adaptive threshold
# scenario_index: index of the scenario to simulate

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
  
  # Generate continuous biomarker x and treatment assignment group
  data <- data.frame(
    subject_id = 1:max_screen,
    x = runif(max_screen, min = 0, max = 1),         # continuous biomarker in [0,1]
    group = rbinom(max_screen, size = 1, prob = 0.5) # equal randomization
  )
  
  # Define parameters
  # Linear case
  beta0 <- 5         # baseline intercept
  beta1 <- 0.5       # baseline slope
  sigma <- 0.1       # residual SD
  # Non-linear case
  baseline  <- 3     # lowest value at x=0 or x=1
  amplitude <- 2     # difference between peak at x=0.5 and lowest value at x=0 or x=1 for control arm
  multiplier <- 1.2  # multiply effect size by 1.2 for scenario 3,4,5
  # Random noise
  noise <- rnorm(max_screen, mean = 0, sd = sigma)
  
  if (scenario_index == 1) {
    # Null case: effect size=0
    # Linear relationship between x and y
    # No treatment effect for all patients
    effect <- 0
    
    y_control <- beta0 + beta1 * data$x + noise
    y_treat   <- y_control + effect
    
    data$y <- ifelse(data$group == 1, y_treat, y_control)
    
  } else if (scenario_index == 2) {
    # Linear relationship between x and y
    # Constant treatment effect
    effect <- effect_size
    
    y_control <- beta0 + beta1 * data$x + noise
    y_treat   <- y_control + effect
    
    data$y <- ifelse(data$group == 1, y_treat, y_control)
    
  } else if (scenario_index == 3) {
    # Linear relationship between x and y
    # Treatment effect increases with x
    effect <- effect_size * data$x * multiplier
    
    y_control <- beta0 + beta1 * data$x + noise
    y_treat   <- y_control + effect
    
    data$y <- ifelse(data$group == 1, y_treat, y_control)
    
  } else if (scenario_index == 4) {
    # Non-linear relationship between x and y (U shape)
    # Treatment effect is largest at two ends, minimal at x=0.5 (U shape)
    effect <- effect_size * 4 * (data$x - 0.5)^2 * multiplier
    
    # Control outcome: U-shaped baseline curve
    y_control <- baseline + amplitude * (4 * (data$x - 0.5)^2) + noise
    y_treat <- y_control + effect
    
    data$y <- ifelse(data$group == 1, y_treat, y_control)
    
  } else if (scenario_index == 5) {
    # Non-linear relationship between x and y (Inverted U shape)
    # Treatment effect is 0 at x=0 and x=1, max at x=0.5 (Inverted U shape)
    effect <- effect_size * (4 * data$x * (1 - data$x)) * multiplier
    
    y_control <- baseline + amplitude * (1 - 4 * (data$x - 0.5)^2) + noise
    y_treat <- y_control + effect
    
    data$y <- ifelse(data$group == 1, y_treat, y_control)
    
  }
  
  # Convert group to factor
  data$group <- factor(data$group, levels = c(0,1), labels = c("Control", "Treatment"))
  
  # Apply design-specific simulation
  if (design_type == "GSE") {
    result <- gse(
      trial_data      = data, 
      max_screen      = max_screen, 
      max_sample_size = max_sample_size, 
      info_rates      = info_rates, 
      alpha           = alpha, 
      beta            = beta, 
      type_alpha      = type_alpha, 
      type_beta       = type_beta, 
      epsilon         = epsilon, 
      d               = d, 
      g               = g
    )
  } else if (design_type == "AE") {
    result <- ae(
      trial_data      = data, 
      max_screen      = max_screen, 
      max_sample_size = max_sample_size, 
      info_rates      = info_rates, 
      alpha           = alpha, 
      beta            = beta, 
      type_alpha      = type_alpha, 
      type_beta       = type_beta, 
      epsilon         = epsilon, 
      d               = d, 
      g               = g
    )
  } else if (design_type == "GSD") {
    result <- gsd(
      trial_data      = data, 
      max_screen      = max_screen, 
      max_sample_size = max_sample_size, 
      info_rates      = info_rates, 
      alpha           = alpha, 
      beta            = beta, 
      type_alpha      = type_alpha, 
      type_beta       = type_beta
    )
  } else if (design_type == "GSE plus") {
    result <- gse_plus(
      trial_data      = data, 
      max_screen      = max_screen, 
      max_sample_size = max_sample_size, 
      info_rates      = info_rates, 
      alpha           = alpha, 
      beta            = beta, 
      type_alpha      = type_alpha, 
      type_beta       = type_beta, 
      epsilon         = epsilon, 
      d               = d, 
      g               = g,
      delta           = delta
    )
  }
  
  # Return the results
  return(result)
}