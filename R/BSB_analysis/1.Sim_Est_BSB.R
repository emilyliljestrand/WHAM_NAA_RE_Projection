#' @title Black Sea Bass Simulation-Estimation Projection Test
#' @description Simulates 1989-2024 catch and index observations from BSB.EM.Y,
#'   fits full (1989-2024) and truncated (1989-2021) estimation models, projects
#'   the truncated model through 2024, and compares projections with the full fit.
#' @name sim_est_bsb

# WHAM and SPASAM.MSE supply the model-fitting and simulation functions used here.
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  wham_lib <- "C:/Users/emily.liljestrand/AppData/Local/R/win-library/4.4/wham_2.1.0.9011"
  if (!file.exists(file.path(wham_lib, "wham", "DESCRIPTION"))) {
    stop(
      paste("WHAM 2.1.0.9011 was not found at", wham_lib),
      call. = FALSE
    )
  }
  # rm(list = ls()) does not unload a WHAM namespace already loaded from another library.
  if ("wham" %in% loadedNamespaces()) {
    if ("package:wham" %in% search()) {
      detach("package:wham", unload = TRUE, character.only = TRUE)
    } else {
      unloadNamespace("wham")
    }
  }
  library(wham, lib.loc = wham_lib)
  library(SPASAM.MSE)
  library(tidyverse)
})

#Read in existing seeds file, or make one if it doesn't exist
if (!file.exists("data/raw/seeds/Sim_Est_BSB_seeds.csv")) {
  set.seed(799291)
  nseeds <- 5000
  r.seed.set <- trunc(1e7 * runif(n = nseeds), 7) +
    trunc(1e3 * runif(n = nseeds), 3)
  write.table(
    r.seed.set,
    file = "data/raw/seeds/Sim_Est_BSB_seeds.csv",
    quote = F,
    row.names = F,
    col.names = F,
    sep = ","
  )
} else {
  r.seed.set <- read.table("data/raw/seeds/Sim_Est_BSB_seeds.csv")
}
# ==============================================================================
# 1) User controls
# ==============================================================================
# Full OM/EM window is 1989-2024. The training EM stops in 2021 and is then
# projected through 2022-2024 for comparison with the full fit.

nreps <- 20
converged_reps <- integer(0)
ssb_bias_results <- list()

extract_projection_comparison <- function(
  quantity,
  full_model,
  projected_model,
  projection_years,
  labels = NULL
) {
  full_values <- full_model$rep[[quantity]]
  projected_values <- projected_model$rep[[quantity]]

  if (is.null(dim(full_values))) {
    full_values <- matrix(full_values, ncol = 1L)
    projected_values <- matrix(projected_values, ncol = 1L)
  }

  full_rows <- match(projection_years, full_model$years)
  projected_rows <- match(projection_years, projected_model$years)
  if (anyNA(full_rows) || anyNA(projected_rows)) {
    stop(
      "Projection years were not found in both fitted model objects.",
      call. = FALSE
    )
  }

  component_labels <- labels %||%
    paste0("Component_", seq_len(ncol(full_values)))
  tidyr::expand_grid(year = projection_years, component = component_labels) %>%
    mutate(
      quantity = quantity,
      full_fit = as.vector(t(full_values[full_rows, , drop = FALSE])),
      projection = as.vector(t(projected_values[
        projected_rows,
        ,
        drop = FALSE
      ])),
      difference = projection - full_fit,
      relative_difference = difference / full_fit
    ) %>%
    select(quantity, component, year, everything())
}

for (i in 1:nreps) {
  simulation_seed <- r.seed.set[[1]][i]
  full_years <- 1989:2024
  training_years <- 1989:2021
  projection_years <- 2022:2024

  # Continue recruitment and NAA random effects (opt = 1)
  projection_options1 <- list(
    proj_R_opt = 1,
    proj_NAA_opt = 1,
    proj_F_opt = rep(3, length(projection_years))
  )
  # Turn off random effects on both Recruitment and NAA (opt = 4 for R, opt = 3 for NAA)
  projection_options2 <- list(
    proj_R_opt = 4,
    proj_NAA_opt = 3,
    proj_F_opt = rep(3, length(projection_years))
  )
  # Average recruitment and NAA random effects over last 5 years (opt = 3 for R, opt = 2 for NAA)
  projection_options3 <- list(
    proj_R_opt = 3,
    proj_NAA_opt = 2,
    proj_F_opt = rep(3, length(projection_years))
  )

  output_dir <- here("output", "BSB_analysis", "1.Sim_Est_BSB")

  stopifnot(
    identical(training_years, full_years[full_years <= max(training_years)]),
    length(projection_years) == 3L,
    min(projection_years) == max(training_years) + 1L,
    file.exists(here("models", "BSB.EM.Y.RDS"))
  )

  # ==============================================================================
  # 2) Simulate one 1989-2024 catch and index dataset from BSB.EM.Y
  # ==============================================================================
  # Treat the fitted BSB assessment as the operating-model truth and draw one
  # catch/index realization from it.

  operating_model <- readRDS(here("models", "BSB.EM.Y.RDS"))

  if (!identical(as.integer(operating_model$years), full_years)) {
    stop("Operating model does not span expected years.", call. = FALSE)
  }

  random_effects <- operating_model$input$random
  # random_effects <- NULL

  # Simulate one OM draw: process error (recruitment/NAA) and observation error
  # (catch and indices). This is the dataset the estimation models will see.
  simulated_model <- SPASAM.MSE::update_om_fn(
    operating_model,
    seed = simulation_seed,
    random = random_effects
  )

  # Also save just the simulated catch/index observations for inspection without
  # loading the full WHAM object.
  simulated_data <- list(
    years = full_years,
    catch = simulated_model$input$data[c(
      "agg_catch",
      "catch_paa",
      "use_agg_catch",
      "use_catch_paa",
      "agg_catch_sigma",
      "catch_Neff"
    )],
    index = simulated_model$input$data[c(
      "agg_index",
      "index_paa",
      "use_indices",
      "use_index_paa",
      "agg_index_sigma",
      "index_Neff"
    )]
  )
  saveRDS(
    simulated_data,
    file.path(
      output_dir,
      paste("simulated_catch_index_1989_2024.sim", i, ".RDS", sep = "")
    )
  )

  # ==============================================================================
  # 3) Create and fit the full and truncated estimation models
  # ==============================================================================

  # make_em_input() needs the list returned by generate_basic_info(). The fitted
  # OM does not store that object (input$basic_info / catch_info / index_info / F
  # are NULL), so rebuild it from OM data. Copy observation-error settings here;
  # simulated catch and index values are filled in later by make_em_input().
  om_data <- simulated_model$input$data

  # BSB configuration: 2 stocks, 2 regions, 4 fleets, 4 indices, 11 seasons, 8 ages.
  n_stocks <- om_data$n_stocks
  n_regions <- om_data$n_regions
  n_fleets <- om_data$n_fleets
  n_indices <- om_data$n_indices
  n_seasons <- om_data$n_seasons
  n_ages <- om_data$n_ages

  # Stock x year x age maturity from the OM (avoids the default life-history ogive).
  user_maturity <- om_data$mature

  # WAA is source x year x age (4 fleets + 4 indices + 2 stocks). Pointers map
  # those slices to fleets, indices, SSB, and M. totcatch is not stored on the OM,
  # so it uses the SSB pointer. generate_basic_info() requires user_waa$waa to be a
  # vector, matrix, or 3-D array.
  user_waa <- list(
    waa = om_data$waa,
    waa_pointer_fleets = om_data$waa_pointer_fleets,
    waa_pointer_indices = om_data$waa_pointer_indices,
    waa_pointer_totcatch = om_data$waa_pointer_ssb,
    waa_pointer_ssb = om_data$waa_pointer_ssb,
    waa_pointer_M = om_data$waa_pointer_M
  )

  # generate_basic_info() wants year-fractions from January 1; WHAM stores
  # within-season fractions plus season pointers, so convert spawn and index timing.
  season_starts <- cumsum(c(0, om_data$fracyr_seasons))
  fracyr_spawn <- om_data$fracyr_SSB[1, 1] +
    season_starts[om_data$spawn_seasons[1]]
  fracyr_indices <- colMeans(om_data$fracyr_indices)

  # Convert WHAM lognormal SDs to CVs for generate_basic_info().
  sigma_to_cv <- function(sigma) {
    sqrt(exp(sigma^2) - 1)
  }

  # Observation-error settings only. make_em_input() later overwrites catch/index
  # observations and year-specific CVs/Neff from the simulated OM draw.
  catch_info <- list(
    catch_cv = sigma_to_cv(colMeans(om_data$agg_catch_sigma)),
    catch_Neff = colMeans(om_data$catch_Neff),
    use_agg_catch = 1,
    use_catch_paa = 1
  )

  # Survey CVs, Neff, catchability, units, and timing from the OM. Simulated
  # index observations are not copied here.
  index_info <- list(
    index_cv = sigma_to_cv(colMeans(om_data$agg_index_sigma)),
    index_Neff = colMeans(om_data$index_Neff),
    fracyr_indices = fracyr_indices,
    q = operating_model$rep$q[1, ],
    use_indices = 1,
    use_index_paa = 1,
    units_indices = om_data$units_indices,
    units_index_paa = om_data$units_index_paa
  )

  # OM F_config = 1, so fully selected fleet F is exp(cumsum(F_pars)). That
  # matches max FAA-at-age and is the year x fleet matrix expected as user_F.
  # Fbar is age-averaged (with extra region/total columns) and should not be used.
  om_F <- t(operating_model$rep$FAA[1:4, , 8])

  # F.year1 / Fhist / Fmax / Fmin are ignored when user_F is supplied, but still
  # stored on par_inputs. Fill them from the OM F series for consistency.
  # F_feedback is unused because n_feedback_years = 0.
  em_info <- SPASAM.MSE::generate_basic_info(
    n_stocks = n_stocks,
    n_regions = n_regions,
    n_indices = n_indices,
    n_fleets = n_fleets,
    n_seasons = n_seasons,
    base.years = full_years,
    n_feedback_years = 0,
    n_ages = n_ages,
    F_info = list(
      F.year1 = om_F[1, 1],
      Fhist = "constant",
      Fmax = max(om_F),
      Fmin = min(om_F),
      change_time = 0.5,
      user_F = om_F,
      F_feedback = NULL
    ),
    catch_info = catch_info,
    index_info = index_info,
    user_waa = user_waa,
    user_maturity = user_maturity,
    fracyr_spawn = fracyr_spawn,
    fracyr_seasons = om_data$fracyr_seasons,
    fleet_regions = om_data$fleet_regions,
    index_regions = om_data$index_regions,
    XSPR_R_opt = om_data$XSPR_R_opt,
    mig_type = 0,
    apply_re_trend = 0,
    apply_mu_trend = 0
  )

  # generate_basic_info() does not set occupancy or names. Copy OM NAA_where
  # (north age-1 cannot recruit in the south; south stock never occupies the north)
  # and stock/region labels.
  em_info$basic_info$NAA_where <- om_data$NAA_where
  em_info$basic_info$region_names <- operating_model$input$region_names
  em_info$basic_info$stock_names <- operating_model$input$stock_names

  stopifnot(
    length(em_info$par_inputs$n_indices) == 1L,
    em_info$par_inputs$n_indices > 0L
  )

  # Joint 2-stock EM with movement estimated, matching the OM spatial structure.
  em_options <- list(
    separate.em = FALSE,
    separate.em.type = NULL,
    do.move = TRUE,
    est.move = TRUE
  )

  # Year-specific catch/index CVs and Neff from the simulated OM. make_em_input()
  # uses these after prepare_wham_input() because generate_basic_info() only
  # carries constant CVs/Neff.
  observation_updates <- list(
    catch = list(
      agg_catch_sigma = simulated_model$input$data$agg_catch_sigma,
      catch_Neff = simulated_model$input$data$catch_Neff
    ),
    index = list(
      agg_index_sigma = simulated_model$input$data$agg_index_sigma,
      index_Neff = simulated_model$input$data$index_Neff
    )
  )

  estimation_model_config <- readRDS(here("models", "BSB.EM.Y.Config.RDS"))

  # Build an EM input for any year window. Biological structure comes from em_info
  # / the OM; catch and index data come from the simulated model.
  make_bsb_em_input <- function(years) {
    SPASAM.MSE::make_em_input(
      om = simulated_model,
      em_info = em_info,
      M_em = NULL,
      sel_em = estimation_model_config$sel,
      NAA_re_em = estimation_model_config$NAA_re,
      move_em = estimation_model_config$move,
      catchability_em = NULL,
      ecov_em = estimation_model_config$ecov,
      em.opt = em_options,
      em_years = years,
      year.use = length(years),
      age_comp_em = estimation_model_config$age_comp,

      aggregate_catch_info = NULL,
      aggregate_index_info = NULL,
      aggregate_weights_info = NULL,
      filter_indices = NULL,
      reduce_region_info = NULL,

      update_catch_info = observation_updates$catch,
      update_index_info = observation_updates$index
    )
  }

  # Full-period EM: benchmark that sees 2022-2024 observations.
  message("Building and fitting the 1989-2024 estimation model")
  full_input <- make_bsb_em_input(full_years)

  full_model <- fit_wham(
    full_input,
    do.sdrep = TRUE,
    do.retro = FALSE,
    do.osa = FALSE,
    do.brps = TRUE,
    MakeADFun.silent = FALSE
  )
  # Training EM: withhold 2022-2024, then project those years in section 4.
  message("Building and fitting the 1989-2021 estimation model")
  training_input <- make_bsb_em_input(training_years)
  training_model <- fit_wham(
    training_input,
    do.sdrep = TRUE,
    do.retro = FALSE,
    do.osa = FALSE,
    do.brps = TRUE,
    MakeADFun.silent = FALSE
  )
  # Match WHAM self_test()'s convergence criterion: successful optimizer
  # convergence, successful sdreport, and an invertible Hessian.
  full_converged <- isTRUE(full_model$opt$convergence == 0) &&
    isTRUE(full_model$is_sdrep) &&
    isTRUE(!full_model$na_sdrep)
  training_converged <- isTRUE(training_model$opt$convergence == 0) &&
    isTRUE(training_model$is_sdrep) &&
    isTRUE(!training_model$na_sdrep)
  both_models_converged <- full_converged && training_converged

  convergence <- tibble(
    model = c("Full fit: 1989-2024", "Training fit: 1989-2021"),
    convergence_code = c(
      full_model$opt$convergence,
      training_model$opt$convergence
    ),
    converged = c(full_converged, training_converged),
    pd_hessian = c(
      !isTRUE(full_model$na_sdrep),
      !isTRUE(training_model$na_sdrep)
    )
  )
  write_csv(
    convergence,
    file.path(output_dir, paste("model_convergence.sim", i, ".csv", sep = ""))
  )

  if (both_models_converged) {
    converged_reps <- c(converged_reps, i)
  }

  if (both_models_converged) {
    # ==============================================================================
    # 4) Project the 1989-2021 fit through 2024
    # ==============================================================================
    # Continue NAA/R random effects and apply F40% in 2022-2024 (see projection_options).
    # Compare this forecast with the full EM, which includes those years of data.

    message("Projecting the 1989-2021 estimation model through 2024")
    projected_model1 <- project_wham(
      training_model,
      proj.opts = projection_options1,
      check.version = FALSE
    )
    projected_model2 <- project_wham(
      training_model,
      proj.opts = projection_options2,
      check.version = FALSE
    )
    projected_model3 <- project_wham(
      training_model,
      proj.opts = projection_options3,
      check.version = FALSE
    )

    ssb_bias_results[[length(ssb_bias_results) + 1L]] <- bind_rows(
      extract_projection_comparison(
        "SSB",
        full_model,
        projected_model1,
        projection_years,
        c("North", "South")
      ) %>%
        mutate(projection_model = "Option 1"),
      extract_projection_comparison(
        "SSB",
        full_model,
        projected_model2,
        projection_years,
        c("North", "South")
      ) %>%
        mutate(projection_model = "Option 2"),
      extract_projection_comparison(
        "SSB",
        full_model,
        projected_model3,
        projection_years,
        c("North", "South")
      ) %>%
        mutate(projection_model = "Option 3")
    ) %>%
      group_by(projection_model, component) %>%
      summarise(
        mean_relative_bias = mean(relative_difference, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(simulation = i, .before = 1)
  }
}

convergence_rate <- length(converged_reps) / nreps
message(
  "Convergence rate: ",
  length(converged_reps),
  " of ",
  nreps,
  " simulations (",
  formatC(100 * convergence_rate, format = "f", digits = 1),
  "%; fraction = ",
  formatC(convergence_rate, format = "f", digits = 3),
  ")"
)

# ==============================================================================
# 5) Save terminal three-year SSB relative bias results
# ==============================================================================

ssb_relative_bias <- bind_rows(ssb_bias_results)
write_csv(
  ssb_relative_bias,
  file.path(output_dir, "ssb_relative_bias_2022_2024.csv")
)

message(
  "Simulation-estimation comparison complete. Results saved in: ",
  output_dir
)
