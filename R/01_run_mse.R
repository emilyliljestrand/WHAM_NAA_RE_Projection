# Run a reproducible SPASAM.MSE experiment.
# Install once with: remotes::install_github("lichengxue/SPASAM.MSE", dependencies = TRUE)

suppressPackageStartupMessages({
  library(wham)
  library(SPASAM.MSE)
})

# ---- Experiment settings -----------------------------------------------------
source(file.path("R", "functions", "project_paths.R"))
source(file.path("config", "mse_settings.R"))

project_dir <- project_path()
results_dir <- create_output_dir("mse_results")

seed <- mse_settings$seed
n_replicates <- mse_settings$n_replicates
percent_fxspr <- mse_settings$percent_fxspr
strategy_names <- paste0("FXSPR_", percent_fxspr)

year_start <- mse_settings$year_start
year_end <- mse_settings$year_end
n_feedback_years <- mse_settings$n_feedback_years
assessment_interval <- mse_settings$assessment_interval

if (!requireNamespace("SPASAM.MSE", quietly = TRUE)) {
  stop(
    "SPASAM.MSE is not installed. Run remotes::install_github(\"lichengxue/SPASAM.MSE\", dependencies = TRUE) first.",
    call. = FALSE
  )
}

if (length(percent_fxspr) != length(strategy_names) || anyDuplicated(strategy_names)) {
  stop("strategy_names must be unique and match percent_fxspr in length.", call. = FALSE)
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

# This is a compact two-stock, two-region example OM. Replace these settings
# with the Black Sea Bass configuration once it has been converted to SPASAM.MSE.
info <- generate_basic_info(
  n_stocks = 2,
  n_regions = 2,
  n_indices = 2,
  n_fleets = 2,
  n_seasons = 1,
  base.years = year_start:year_end,
  n_feedback_years = n_feedback_years,
  life_history = "short",
  n_ages = 10
)

n_ages <- 10L
naa_sigma <- array(0.2, dim = c(2, 2, n_ages))
naa_sigma[, , 1] <- 0.5
naa_re <- list(
  N1_model = rep("equilibrium", 2),
  sigma = rep("rec+1", 2),
  cor = rep("iid", 2),
  recruit_model = 2,
  sigma_vals = naa_sigma
)

om_input <- prepare_wham_input(
  basic_info = info$basic_info,
  NAA_re = naa_re,
  catch_info = info$catch_info,
  index_info = info$index_info,
  F = info$F
)
random <- om_input$random
om_input$random <- NULL
om <- fit_wham(om_input, do.fit = FALSE, do.brps = TRUE, MakeADFun.silent = TRUE)

base_years <- year_start:year_end
assessment_years <- seq(
  from = year_end,
  to = tail(om$years, 1) - assessment_interval,
  by = assessment_interval
)

metadata <- list(
  seed = seed,
  n_replicates = n_replicates,
  strategies = stats::setNames(percent_fxspr, strategy_names),
  base_years = base_years,
  assessment_years = assessment_years,
  assessment_interval = assessment_interval,
  package_versions = list(
    wham = as.character(utils::packageVersion("wham")),
    SPASAM.MSE = as.character(utils::packageVersion("SPASAM.MSE"))
  ),
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(metadata, file.path(results_dir, "run_metadata.rds"))

for (strategy_index in seq_along(percent_fxspr)) {
  strategy_dir <- file.path(results_dir, strategy_names[[strategy_index]])
  dir.create(strategy_dir, recursive = TRUE, showWarnings = FALSE)

  for (replicate_id in seq_len(n_replicates)) {
    replicate_seed <- seed + strategy_index * 10000L + replicate_id
    output_file <- file.path(strategy_dir, sprintf("replicate_%03d.rds", replicate_id))

    message(sprintf("Running %s, replicate %d of %d", strategy_names[[strategy_index]], replicate_id, n_replicates))
    tryCatch({
      om_with_data <- update_om_fn(om, seed = replicate_seed, random = random)
      mse_result <- loop_through_fn(
        om = om_with_data,
        em_info = info,
        random = random,
        NAA_re_em = naa_re,
        em.opt = list(separate.em = FALSE, separate.em.type = 1, do.move = FALSE, est.move = FALSE),
        assess_years = assessment_years,
        assess_interval = assessment_interval,
        base_years = base_years,
        year.use = length(base_years),
        hcr = list(hcr.type = 1, hcr.opts = list(use_FXSPR = TRUE, percentFXSPR = percent_fxspr[[strategy_index]])),
        seed = replicate_seed,
        save.last.em = FALSE
      )
      saveRDS(mse_result, output_file)
    }, error = function(error) {
      saveRDS(
        list(error = conditionMessage(error), seed = replicate_seed),
        sub("\\.rds$", "_error.rds", output_file)
      )
      message(sprintf("Failed %s replicate %d: %s", strategy_names[[strategy_index]], replicate_id, conditionMessage(error)))
    })
  }
}

message("MSE run complete. Results saved in: ", results_dir)