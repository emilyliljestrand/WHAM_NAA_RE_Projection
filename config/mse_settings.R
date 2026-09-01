#' @title Management Strategy Evaluation (MSE) Settings Configuration
#' @description Central configuration list containing seed values, simulation timeline control, 
#'   replicate counts, and fishing target options for SPASAM.MSE experiments.
#' @details Sourced by `R/01_run_mse.R` to establish reproducible simulation parameters 
#'   for spatial stock assessment projections.
#' @name mse_settings
NULL

#' Configuration settings list for SPASAM.MSE runs
#'
#' @field seed Integer seed for global pseudo-random number generation ensuring reproducibility.
#' @field n_replicates Integer number of simulation replicates per harvest control rule strategy.
#' @field percent_fxspr Numeric vector of target F_X%SPR percentages (e.g., 60% F_XSPR, 80% F_XSPR).
#' @field year_start Starting year (integer) of the historical base period.
#' @field year_end Ending year (integer) of the historical base period.
#' @field n_feedback_years Number of projection/feedback years in the MSE horizon.
#' @field assessment_interval Frequency (in years) between successive stock assessment updates.
mse_settings <- list(
  # Global random seed for reproducible simulation replicates
  seed = 20260901L,
  
  # Number of Monte Carlo simulation replicates per strategy
  n_replicates = 10L,
  
  # Target F%XSPR levels to evaluate in harvest control rules
  percent_fxspr = c(60, 80),
  
  # Base historical period start and end years
  year_start = 1L,
  year_end = 30L,
  
  # Duration of MSE projection feedback period in years
  n_feedback_years = 10L,
  
  # Assessment frequency interval (every N years)
  assessment_interval = 2L
)