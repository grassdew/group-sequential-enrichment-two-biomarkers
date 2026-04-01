# =============================================================================
# Simulation driver: 1 binary biomarker (x1) + 1 continuous biomarker (x2)
#
# Grid:
#   prev_rate_x1 in {0.3, 0.5, 0.7}
#   effect_size  in {0.1, 0.2, 0.3}
#   scenario_index in {1,2,3,4,5}
#   design_type in {"GSE","AE","GSD","GSE-F"}
#
# IMPORTANT:
# - Update `source()` path below to where you save the design-functions script.
# - Update `root_out_dir` to your preferred output folder.
# =============================================================================

library(foreach)
library(doParallel)
library(doRNG)
library(data.table)

# pkgs used inside workers
pkgs <- c("rstanarm", "rpact", "tibble", "data.table")

# ---- user paths --------------------------------------------------------------
root_out_dir <- "/Users/emily/Documents/Adaptive Enrichment Trial Design/Binary and continuous biomarker/Simulation results"

# ---- parallel setup ----------------------------------------------------------
num_cores <- max(1L, parallel::detectCores() - 1L)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  library(rstanarm)
  library(rpact)
  library(tibble)
  library(data.table)
  Sys.setenv(OMP_NUM_THREADS = "1")
  options(mc.cores = 1)
  source("/Users/emily/Documents/Adaptive Enrichment Trial Design/Binary and continuous biomarker/Code/Design functions.R")
})

registerDoRNG(303)

# ---- global parameters -------------------------------------------------------
num_sim         <- 1000
chunk_size      <- 10      # tune 10–50
max_screen      <- 5000
max_sample_size <- 600
info_rates      <- c(1/3, 2/3, 1)

alpha      <- 0.025
beta       <- 0.2
type_alpha <- "asOF"
type_beta  <- "bsOF"

epsilon <- 0.1
d <- 0.5
g <- 0.1

# GSE-F only: allow up to delta extra projected screens beyond max_screen
delta <- 0

design_types <- "GSE"
prev_grid    <- 0.5
ef_grid      <- 0.2
scenarios    <- 4

# ---- output dir helper -------------------------------------------------------
make_out_dir <- function(design_type) {
  out_dir <- file.path(
    root_out_dir,
    sprintf("Prior N(0, 10^4), d=%s, g=%s", d, g),
    sprintf("%d interim analysis", length(info_rates) - 1),
    design_type
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

# ---- run grid ----------------------------------------------------------------
for (design_type in design_types) {

  out_dir <- make_out_dir(design_type)

  for (scenario_index in scenarios) {
    for (prev_rate_x1 in prev_grid) {
      for (effect_size in ef_grid) {

        n_chunks <- ceiling(num_sim / chunk_size)

        res_list <- foreach(
          chunk_id = 1:n_chunks,
          .packages = pkgs,
          .multicombine = TRUE,
          .maxcombine = 50
        ) %dorng% {

          rows <- vector("list", length = min(chunk_size, num_sim - (chunk_id - 1L) * chunk_size))
          i0 <- (chunk_id - 1L) * chunk_size

          for (j in seq_along(rows)) {
            sim_index <- i0 + j

            trial_result <- trial_simulation(
              sim_index = sim_index,
              design_type = design_type,
              prev_rate_x1 = prev_rate_x1,
              max_screen = max_screen,
              max_sample_size = max_sample_size,
              info_rates = info_rates,
              effect_size = effect_size,
              alpha = alpha,
              beta = beta,
              type_alpha = type_alpha,
              type_beta = type_beta,
              epsilon = epsilon,
              d = d,
              g = g,
              delta = delta,
              scenario_index = scenario_index
            )

            dt <- as.data.table(as.list(trial_result))
            rows[[j]] <- dt
          }

          rbindlist(rows, use.names = TRUE, fill = TRUE)
        }

        sim_results_dt <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

        file_name <- file.path(
          out_dir,
          sprintf("prev_%s_efsize_%s_scenario_%s.csv", prev_rate_x1, effect_size, scenario_index)
        )
        fwrite(sim_results_dt, file_name)

      }
    }
  }
}

stopCluster(cl)
