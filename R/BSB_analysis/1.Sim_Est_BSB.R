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
# Full OM/EM window is 1989-2024. The training EM is then projected through the
# user-specified projection years for comparison with the full fit.

nreps <- 100
converged_reps <- integer(0)
ssb_bias_results <- list()
save_simulated_data <- FALSE
save_model_outputs <- FALSE

full_years <- 1989:2024
training_years <- 1989:2021
projection_years <- 2022:2024
output_dir <- here(
  "output",
  "BSB_analysis",
  "1.Sim_Est_BSB",
  paste0("projection_years_", length(projection_years))
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Resume interrupted runs. The convergence CSV acts as the completion marker for
# a replicate, whether or not the full or training model converged.
convergence_files <- list.files(
  output_dir,
  pattern = "^model_convergence[.]sim[0-9]+[.]csv$",
  full.names = TRUE
)
completed_ids <- as.integer(sub(
  "^.*[.]sim([0-9]+)[.]csv$",
  "\\1",
  convergence_files
))
completed_ids <- completed_ids[completed_ids >= 1L & completed_ids <= nreps]

completed_reps <- unique(completed_ids)

if (length(completed_reps) > 0L) {
  converged_reps <- completed_reps[vapply(
    completed_reps,
    function(rep_id) {
      convergence <- read_csv(
        file.path(output_dir, paste0("model_convergence.sim", rep_id, ".csv")),
        show_col_types = FALSE
      )
      all(convergence$converged)
    },
    logical(1)
  )]
}

replicates_to_run <- setdiff(seq_len(nreps), completed_reps)
message(
  "Resuming simulation-estimation run: ",
  length(completed_reps),
  " completed; ",
  length(replicates_to_run),
  " remaining."
)

ssb_bias_file <- file.path(output_dir, "ssb_relative_bias_2022_2024.csv")
replicate_bias_files <- list.files(
  output_dir,
  pattern = "^ssb_relative_bias[.]sim[0-9]+[.]csv$",
  full.names = TRUE
)
if (length(replicate_bias_files) > 0L) {
  ssb_bias_results <- lapply(
    replicate_bias_files,
    read_csv,
    show_col_types = FALSE
  )
} else if (file.exists(ssb_bias_file)) {
  ssb_bias_results[[1]] <- read_csv(ssb_bias_file, show_col_types = FALSE)
}

if (length(r.seed.set[[1]]) < nreps) {
  stop("The seed file does not contain enough seeds for nreps.", call. = FALSE)
}

projection_options1 <- list(
  proj_R_opt = 1,
  proj_NAA_opt = 1,
  proj_F_opt = rep(3, length(projection_years))
)
projection_options2 <- list(
  proj_R_opt = 4,
  proj_NAA_opt = 3,
  proj_F_opt = rep(3, length(projection_years))
)
projection_options3 <- list(
  proj_R_opt = 3,
  proj_NAA_opt = 2,
  proj_F_opt = rep(3, length(projection_years))
)

stopifnot(
  identical(training_years, full_years[full_years <= max(training_years)]),
  length(projection_years) >= 1L,
  min(projection_years) == max(training_years) + 1L,
  file.exists(here("models", "BSB.EM.Y.RDS")),
  file.exists(here("models", "BSB.EM.Y.Config.RDS"))
)

operating_model <- readRDS(here("models", "BSB.EM.Y.RDS"))
if (!identical(as.integer(operating_model$years), full_years)) {
  stop("Operating model does not span expected years.", call. = FALSE)
}
random_effects <- operating_model$input$random

estimation_model_config <- readRDS(here("models", "BSB.EM.Y.Config.RDS"))
em_options <- list(
  separate.em = FALSE,
  separate.em.type = NULL,
  do.move = TRUE,
  est.move = TRUE
)
sigma_to_cv <- function(sigma) {
  sqrt(exp(sigma^2) - 1)
}

# Prepare invariant EM structure once; replicate-specific observations and
# observation-error values are supplied inside the loop.
om_data <- operating_model$input$data
n_stocks <- om_data$n_stocks
n_regions <- om_data$n_regions
n_fleets <- om_data$n_fleets
n_indices <- om_data$n_indices
n_seasons <- om_data$n_seasons
n_ages <- om_data$n_ages
user_maturity <- om_data$mature
user_waa <- list(
  waa = om_data$waa,
  waa_pointer_fleets = om_data$waa_pointer_fleets,
  waa_pointer_indices = om_data$waa_pointer_indices,
  waa_pointer_totcatch = om_data$waa_pointer_ssb,
  waa_pointer_ssb = om_data$waa_pointer_ssb,
  waa_pointer_M = om_data$waa_pointer_M
)
season_starts <- cumsum(c(0, om_data$fracyr_seasons))
fracyr_spawn <- om_data$fracyr_SSB[1, 1] +
  season_starts[om_data$spawn_seasons[1]]
fracyr_indices <- colMeans(om_data$fracyr_indices)
catch_info <- list(
  catch_cv = sigma_to_cv(colMeans(om_data$agg_catch_sigma)),
  catch_Neff = colMeans(om_data$catch_Neff),
  use_agg_catch = 1,
  use_catch_paa = 1
)
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
om_F <- t(operating_model$rep$FAA[1:4, , 8])
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
em_info$basic_info$NAA_where <- om_data$NAA_where
em_info$basic_info$region_names <- operating_model$input$region_names
em_info$basic_info$stock_names <- operating_model$input$stock_names
stopifnot(
  length(em_info$par_inputs$n_indices) == 1L,
  em_info$par_inputs$n_indices > 0L
)

extract_projection_ssb_comparison <- function(
  full_model,
  projected_model,
  projection_years
) {
  full_ssb <- full_model$rep$SSB
  projected_ssb <- projected_model$rep$SSB
  comparison_years <- tail(full_model$years, length(projection_years))

  full_rows <- match(comparison_years, full_model$years)
  projected_rows <- match(comparison_years, projected_model$years)
  if (anyNA(full_rows) || anyNA(projected_rows)) {
    stop(
      "The comparison years were not found in both fitted model objects.",
      call. = FALSE
    )
  }

  tibble(
    region = c("North", "South"),
    mean_relative_error = colMeans(
      (projected_ssb[projected_rows, , drop = FALSE] -
        full_ssb[full_rows, , drop = FALSE]) /
        full_ssb[full_rows, , drop = FALSE],
      na.rm = TRUE
    )
  )
}

for (i in replicates_to_run) {
  run_number <- match(i, replicates_to_run)
  message(
    "Starting replicate ",
    i,
    " (",
    run_number,
    " of ",
    length(replicates_to_run),
    ")"
  )
  simulation_seed <- r.seed.set[[1]][i]

  # ==============================================================================
  # 2) Simulate one 1989-2024 catch and index dataset from BSB.EM.Y
  # ==============================================================================
  # Treat the fitted BSB assessment as the operating-model truth and draw one
  # catch/index realization from it.

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
  if (save_simulated_data) {
    saveRDS(
      simulated_data,
      file.path(
        output_dir,
        paste("simulated_catch_index_1989_2024.sim", i, ".RDS", sep = "")
      )
    )
  }

  # ==============================================================================
  # 3) Create and fit the full and truncated estimation models
  # ==============================================================================

  # Year-specific observation-error values vary by replicate; the biological
  # and structural EM inputs below are prepared once outside this loop.
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
  message("Replicate ", i, ": fitting full 1989-2024 model")
  full_input <- make_bsb_em_input(full_years)

  full_model <- fit_wham(
    full_input,
    do.sdrep = TRUE,
    do.retro = FALSE,
    do.osa = FALSE,
    do.brps = FALSE,
    MakeADFun.silent = TRUE
  )

  # Match WHAM self_test()'s convergence criterion before fitting the training
  # model, since projection comparisons require a converged full fit.
  full_converged <- isTRUE(full_model$opt$convergence == 0) &&
    isTRUE(full_model$is_sdrep) &&
    isTRUE(!full_model$na_sdrep)
  if (!full_converged) {
    convergence <- tibble(
      model = "Full fit: 1989-2024",
      convergence_code = full_model$opt$convergence,
      converged = FALSE,
      pd_hessian = !isTRUE(full_model$na_sdrep)
    )
    write_csv(
      convergence,
      file.path(output_dir, paste("model_convergence.sim", i, ".csv", sep = ""))
    )
    message(
      "Replicate ",
      i,
      ": full model did not converge; skipping remaining fits"
    )
    next
  }

  # Training EM: withhold 2022-2024, then project those years in section 4.
  message("Replicate ", i, ": fitting training 1989-2021 model")
  training_input <- make_bsb_em_input(training_years)
  training_model <- fit_wham(
    training_input,
    do.sdrep = TRUE,
    do.retro = FALSE,
    do.osa = FALSE,
    do.brps = TRUE,
    MakeADFun.silent = TRUE
  )
  # Match WHAM self_test()'s convergence criterion: successful optimizer
  # convergence, successful sdreport, and an invertible Hessian.
  training_converged <- isTRUE(training_model$opt$convergence == 0) &&
    isTRUE(training_model$is_sdrep) &&
    isTRUE(!training_model$na_sdrep)
  both_models_converged <- full_converged && training_converged
  message(
    "Replicate ",
    i,
    ": convergence full=",
    full_converged,
    ", training=",
    training_converged
  )

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
  if (both_models_converged) {
    converged_reps <- c(converged_reps, i)
  } else {
    message(
      "Replicate ",
      i,
      ": full and/or training model did not converge; skipping remaining work"
    )
  }

  if (both_models_converged) {
    # ==============================================================================
    # 4) Project the 1989-2021 fit through 2024
    # ==============================================================================
    # Continue NAA/R random effects and apply F40% in 2022-2024 (see projection_options).
    # Compare this forecast with the full EM, which includes those years of data.

    message("Replicate ", i, ": running three projection options")
    projected_model1 <- project_wham(
      training_model,
      proj.opts = projection_options1,
      do.sdrep = FALSE,
      check.version = FALSE
    )
    projected_model2 <- project_wham(
      training_model,
      proj.opts = projection_options2,
      do.sdrep = FALSE,
      check.version = FALSE
    )
    projected_model3 <- project_wham(
      training_model,
      proj.opts = projection_options3,
      do.sdrep = FALSE,
      check.version = FALSE
    )

    if (save_model_outputs) {
      saveRDS(
        full_model,
        file.path(output_dir, paste0("full_model.sim", i, ".RDS"))
      )
      saveRDS(
        projected_model1,
        file.path(output_dir, paste0("projected_model1.sim", i, ".RDS"))
      )
      saveRDS(
        projected_model2,
        file.path(output_dir, paste0("projected_model2.sim", i, ".RDS"))
      )
      saveRDS(
        projected_model3,
        file.path(output_dir, paste0("projected_model3.sim", i, ".RDS"))
      )
    }

    replicate_bias <- bind_rows(
      extract_projection_ssb_comparison(
        full_model,
        projected_model1,
        projection_years
      ) %>%
        mutate(projection_model = "Option 1"),
      extract_projection_ssb_comparison(
        full_model,
        projected_model2,
        projection_years
      ) %>%
        mutate(projection_model = "Option 2"),
      extract_projection_ssb_comparison(
        full_model,
        projected_model3,
        projection_years
      ) %>%
        mutate(projection_model = "Option 3")
    ) %>%
      mutate(simulation = i, .before = 1)
    ssb_bias_results[[length(ssb_bias_results) + 1L]] <- replicate_bias
    write_csv(
      replicate_bias,
      file.path(output_dir, paste0("ssb_relative_bias.sim", i, ".csv"))
    )
  }

  # Write the completion marker last. If a run stops before this point, the
  # replicate will be retried the next time the script is started.
  write_csv(
    convergence,
    file.path(output_dir, paste("model_convergence.sim", i, ".csv", sep = ""))
  )
  message(
    "Finished replicate ",
    i,
    " (",
    run_number,
    " of ",
    length(replicates_to_run),
    ")"
  )
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
# 5) Save terminal projection-horizon SSB relative bias results
# ==============================================================================

ssb_relative_bias <- bind_rows(ssb_bias_results)
write_csv(
  ssb_relative_bias,
  ssb_bias_file
)

message(
  "Simulation-estimation comparison complete. Results saved in: ",
  output_dir
)
