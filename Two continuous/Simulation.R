# =============================================================================
# Simulation driver: 2 continuous biomarkers (x1, x2 ~ Uniform(0,1))
#
# Grid:
#   effect_size     in ef_grid
#   scenario_index  in scenario_indices
#   design_type     in design_types ("GSE", "AE", "GSD", "GSE-F")
#
# Notes:
# - Update `design_fn_path` to point to your Design functions_2biomarkers.R file.
# - Update `root_out_dir` to your preferred output folder.
# - Output: one CSV per (design_type, effect_size, scenario_index).
# =============================================================================

library(foreach)
library(doParallel)
library(doRNG)
library(data.table)

# ---- user paths --------------------------------------------------------------
# Where the 2-biomarker design functions live
design_fn_path <- "/Users/emily/Documents/Adaptive Enrichment Trial Design/Continuous biomarker/Two biomarkers/Code/Design functions.R"

# Base output folder
root_out_dir <- "/Users/emily/Documents/Adaptive Enrichment Trial Design/Continuous biomarker/Two biomarkers/Simulation results"

# ---- parallel setup ----------------------------------------------------------
num_cores <- max(1L, parallel::detectCores() - 1L)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Export the design-functions path to workers
clusterExport(cl, varlist = c("design_fn_path"), envir = environment())

# Load pkgs + source ONCE per worker; set threading to 1
clusterEvalQ(cl, {
  library(rstanarm)
  library(rpact)
  library(tibble)
  library(data.table)

  Sys.setenv(OMP_NUM_THREADS = "1")  # avoid nested threading
  options(mc.cores = 1)

  source(design_fn_path)
})

# ---- packages used inside foreach workers -----------------------------------
pkgs <- c("rstanarm", "rpact", "tibble", "data.table")

# ---- global parameters -------------------------------------------------------
num_sim         <- 1000
chunk_size      <- 10      # tune 10–50
max_screen      <- 5000
max_sample_size <- 600
# Vary number of interim analyses via equally-spaced information rates
#   1 interim: info_rates = c(0.5, 1)
#   2 interims: info_rates = c(1/3, 2/3, 1)
#   3 interims: info_rates = c(0.25, 0.5, 0.75, 1)
info_rates_list <- list(
  c(0.5, 1),
  c(1/3, 2/3, 1),
  c(0.25, 0.5, 0.75, 1)
)


alpha      <- 0.025
beta       <- 0.2
type_alpha <- "asOF"
type_beta  <- "bsOF"

epsilon <- 0.1
d <- 0.5
g <- 0.1

# GSE-F only: allow up to delta extra projected screens beyond max_screen
delta <- 0

# ---- grids -------------------------------------------------------------------
design_types    <- c("GSE", "AE", "GSD", "GSE-F")
ef_grid         <- c(0.1, 0.2, 0.3)
scenario_indices <- 1:5

# reproducibility across parallel workers
base_seed <- 303

# ---- output dir helper -------------------------------------------------------
make_out_dir <- function(design_type, info_rates) {
  out_dir <- file.path(
    root_out_dir,
    sprintf("Prior N(0, 10), d=%s, g=%s", d, g),
    sprintf("%d interim analysis", length(info_rates) - 1),
    design_type
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

# ---- run grid ----------------------------------------------------------------
for (info_rates in info_rates_list) {

  n_interim <- length(info_rates) - 1L

  for (design_type in design_types) {

    out_dir <- make_out_dir(design_type, info_rates)

    for (scenario_index in scenario_indices) {

      # Scenario 1 is the null: force effect_size = 0 (no looping over ef_grid)
      ef_this <- if (scenario_index == 1) 0 else ef_grid

      for (effect_size in ef_this) {

        n_chunks <- ceiling(num_sim / chunk_size)

        # Unique seed per (n_interim, design_type, scenario, effect_size)
        grid_seed <- base_seed +
          1000000L * as.integer(n_interim) +
          10000L   * match(design_type, design_types) +
          100L     * as.integer(scenario_index) +
          1L       * as.integer(round(effect_size * 100))

        res_list <- foreach(
          chunk_id = 1:n_chunks,
          .packages = pkgs,
          .multicombine = TRUE,
          .maxcombine = 50,
          .options.RNG = grid_seed
        ) %dorng% {

          # Safety: ensure trial_simulation exists on worker
          if (!exists("trial_simulation")) {
            source(design_fn_path)
          }

          rows <- vector("list", length = min(chunk_size, num_sim - (chunk_id - 1L) * chunk_size))
          i0   <- (chunk_id - 1L) * chunk_size

          for (j in seq_along(rows)) {
            sim_index <- i0 + j

            trial_result <- trial_simulation(
              sim_index       = sim_index,
              design_type     = design_type,
              max_screen      = max_screen,
              max_sample_size = max_sample_size,
              info_rates      = info_rates,
              effect_size     = effect_size,
              alpha           = alpha,
              beta            = beta,
              type_alpha      = type_alpha,
              type_beta       = type_beta,
              epsilon         = epsilon,
              d               = d,
              g               = g,
              delta           = delta,
              scenario_index  = scenario_index
            )

            dt <- as.data.table(as.list(trial_result))
            dt[, `:=`(
              sim_index   = sim_index,
              scenario    = scenario_index,
              effect_size = effect_size,
              design_type = design_type,
              n_interim   = as.integer(n_interim),
              info_rates  = paste(info_rates, collapse = ",")
            )]
            rows[[j]] <- dt
          }

          rbindlist(rows, use.names = TRUE, fill = TRUE)
        }

        sim_results_dt <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

        ef_label <- format(effect_size, trim = TRUE, scientific = FALSE)

        file_name <- file.path(
          out_dir,
          sprintf("efsize_%s_scenario_%s.csv", ef_label, scenario_index)
        )
        fwrite(sim_results_dt, file_name)
      }
    }
  }
}

stopCluster(cl)

