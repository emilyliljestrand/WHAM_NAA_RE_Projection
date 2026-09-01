#' @title Generate MSE Summary Report and Plots
#' @description Loads results from completed SPASAM.MSE simulations (`R/01_run_mse.R`) 
#'   and compiles an HTML report comparing management strategies.
#' @details Uses `SPASAM.MSE::plot_mse_output()` to aggregate replicate runs across target 
#'   F%XSPR levels, building visual summaries of spawning stock biomass (SSB), recruitment, catch, 
#'   and relative performance metrics.
#' @author Emily Liljestrand
#' @name plot_mse
NULL

#' # ==============================================================================
#' # 1) Environment Setup & Verification
#' # ==============================================================================

# Attach required package quietly
suppressPackageStartupMessages(library(SPASAM.MSE))

# Resolve absolute paths for input results and target report output directory
project_dir <- normalizePath(file.path(getwd()), mustWork = TRUE)
results_dir <- file.path(project_dir, "output", "mse_results")
report_dir <- file.path(project_dir, "output", "mse_report")

# Check dependency installation
if (!requireNamespace("SPASAM.MSE", quietly = TRUE)) {
  stop(
    "SPASAM.MSE is not installed. Run remotes::install_github(\"lichengxue/SPASAM.MSE\", dependencies = TRUE) first.",
    call. = FALSE
  )
}
if (!file.exists(file.path(results_dir, "run_metadata.rds"))) {
  stop("No MSE results found. Run R/01_run_mse.R before plotting.", call. = FALSE)
}

# Load run metadata and extract strategy names
metadata <- readRDS(file.path(results_dir, "run_metadata.rds"))
strategy_names <- names(metadata$strategies)

#' # ==============================================================================
#' # 2) Replicate Data Loader Function
#' # ==============================================================================

#' Load Individual MSE Replicate Result
#'
#' @description Helper function to read a specific RDS replicate file for a given management strategy.
#' @param strategy_name Character name of the harvest strategy folder (e.g., "FXSPR_60").
#' @param replicate_id Integer ID of the simulation replicate to load.
#' @return A list object containing the full simulation output for that replicate.
#' @export
load_replicate <- function(strategy_name, replicate_id) {
  # Build target file path
  result_file <- file.path(results_dir, strategy_name, sprintf("replicate_%03d.rds", replicate_id))
  if (!file.exists(result_file)) {
    stop(sprintf("Missing result: %s", result_file), call. = FALSE)
  }
  readRDS(result_file)
}

#' # ==============================================================================
#' # 3) Aggregate Replicates across Strategies
#' # ==============================================================================

# Build nested list of model results: outer list over replicates, inner list over strategies
mods <- lapply(seq_len(metadata$n_replicates), function(replicate_id) {
  strategy_results <- lapply(strategy_names, load_replicate, replicate_id = replicate_id)
  names(strategy_results) <- strategy_names
  strategy_results
})

#' # ==============================================================================
#' # 4) Render Performance Plots and Report
#' # ==============================================================================

# Create report folder
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

# Generate HTML performance report using SPASAM.MSE
plot_mse_output(
  mods,
  main_dir = project_dir,
  output_dir = "output/mse_report",
  output_format = "html",
  width = 10,
  height = 7,
  dpi = 300,
  col.opt = "D",
  new_model_names = strategy_names,
  base.model = strategy_names[[1]],
  start.years = min(metadata$assessment_years) + 1,
  use.n.years.first = 5,
  use.n.years.last = 5
)

message("MSE report created in: ", report_dir)