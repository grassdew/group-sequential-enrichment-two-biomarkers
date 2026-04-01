library(foreach)
library(doParallel)
library(doRNG)
library(data.table)   # for fast bind/write

# ---- cluster setup ----------------------------------------------------------
num_cores <- max(1L, parallel::detectCores() - 1L)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# load pkgs + source ONCE per worker; set threading to 1
clusterEvalQ(cl, {
  library(rstanarm)
  library(rpact)
  library(tibble)
  # If you add a "fast mode", also: library(mgcv); library(mvtnorm)
  Sys.setenv(OMP_NUM_THREADS = "1")  # avoid nested threading
  options(mc.cores = 1)
  source("/Users/emily/Documents/Adaptive Enrichment Trial Design/Continuous biomarker/One biomarker/Code/Design functions.R")
})

# ---- params -----------------------------------------------------------------
num_sim       <- 1000
chunk_size    <- 10        # run 10 sims per worker task (tune 10–50)
design_type   <- "GSE"
max_screen    <- 5000
max_sample_size <- 600
info_rates    <- c(1/3, 2/3, 1)
effect_size   <- 0.2
alpha         <- 0.025
beta          <- 0.2
type_alpha    <- "asOF"
type_beta     <- "bsOF"
epsilon       <- 0.1
d             <- 0.5
g             <- 0.1
delta         <- 0
scenario_indices <- 4

# reproducibility across parallel workers
registerDoRNG(303)

# ---- helper: ensure output dir exists --------------------------------------
out_dir <- file.path("/Users/emily/Documents/Adaptive Enrichment Trial Design/Continuous biomarker/One biomarker/Simulation results",
                     sprintf("Prior N(0, 10^4), d=%s, g=%s", d, g),
                     sprintf("%d interim analysis", length(info_rates) - 1),
                     design_type)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- run sims ---------------------------------------------------------------
for (scenario_index in scenario_indices) {
  
  n_chunks <- ceiling(num_sim / chunk_size)
  
  # OPTIONAL: precompute group-seq bounds once and pass them in (requires trivial edits to your design fns)
  # bounds <- stopping_bounds(alpha, beta, type_alpha, type_beta, info_rates)
  
  res_list <- foreach(
    chunk_id = 1:n_chunks,
    .packages = c("rstanarm", "tibble", "rpact", "data.table"),
    .multicombine = TRUE,
    .maxcombine = 50
  ) %dorng% {
    
    # local accumulator
    rows <- vector("list", length = min(chunk_size, num_sim - (chunk_id - 1L) * chunk_size))
    i0   <- (chunk_id - 1L) * chunk_size
    
    for (j in seq_along(rows)) {
      sim_index <- i0 + j
      
      # IMPORTANT: don't set.seed() here; doRNG already did.
      
      trial_result <- trial_simulation(
        sim_index    = sim_index,
        design_type  = design_type,
        max_screen   = max_screen,
        max_sample_size = max_sample_size,
        info_rates   = info_rates,
        effect_size  = effect_size,
        alpha        = alpha,
        beta         = beta,
        type_alpha   = type_alpha,
        type_beta    = type_beta,
        epsilon      = epsilon,
        d            = d,
        g            = g,
        delta        = delta,
        scenario_index = scenario_index
        # If you precomputed bounds:
        # bounds = bounds
      )
      
      # coerce once; data.table is fast for rbind later
      rows[[j]] <- as.data.table(as.list(trial_result))
      # optionally add sim metadata:
      rows[[j]][, `:=`(sim_index = sim_index, scenario = scenario_index)]
    }
    
    rbindlist(rows, use.names = TRUE, fill = TRUE)
  }
  
  sim_results_dt <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  
  # write CSV quickly
  file_name <- file.path(
    out_dir,
    sprintf("efsize_%s_scenario_%s.csv", effect_size, scenario_index)
  )
  fwrite(sim_results_dt, file_name)
}

stopCluster(cl)
