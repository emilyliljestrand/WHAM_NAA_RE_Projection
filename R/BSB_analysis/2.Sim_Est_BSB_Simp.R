#' @title Black Sea Bass Simulation-Estimation Projection Test
#' @description Simulates 1989-2024 catch and index observations from BSB.EM.Y,
#'   fits a reduced (1989-2021) estimation model, projects converged fits through
#'   2024, and compares projected SSB with the simulated operating-model truth.
#' @name sim_est_bsb

# WHAM and SPASAM.MSE supply the model-fitting and simulation functions used here.
rm(list = ls())

suppressPackageStartupMessages({
  library(here)
  wham_lib <- "C:/Users/emily.liljestrand/AppData/Local/R/win-library/4.4/wham_2.1.0.9011"
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
# user-specified projection years for comparison with simulated SSB.

nreps <- 50
replicates_to_run <- seq_len(nreps)
converged_reps <- integer(0)
ssb_bias_results <- list()
save_simulation_outputs <- TRUE

om_disable_movement <- FALSE
om_change_naa_sd <- TRUE
om_recruitment_sd <- 0.5
om_older_ages_sd <- 0.1

full_years <- 1989:2024
training_years <- 1989:2021
projection_years <- 2022:2024
stopifnot(
  is.logical(om_disable_movement),
  length(om_disable_movement) == 1L,
  !is.na(om_disable_movement),
  is.logical(om_change_naa_sd),
  length(om_change_naa_sd) == 1L,
  !is.na(om_change_naa_sd)
)
output_name <- paste0("projection_years_", length(projection_years))
if (om_disable_movement) {
  output_name <- paste0(output_name, "_movement_disabled")
}
if (om_change_naa_sd) {
  stopifnot(
    length(om_recruitment_sd) == 1L,
    length(om_older_ages_sd) == 1L,
    is.finite(om_recruitment_sd),
    is.finite(om_older_ages_sd),
    om_recruitment_sd > 0,
    om_older_ages_sd > 0
  )
  output_name <- paste0(
    output_name,
    "_recruitment_sd_",
    om_recruitment_sd,
    "_older_ages_sd_",
    om_older_ages_sd
  )
}
output_dir <- here(
  "output",
  "BSB_analysis",
  "2.Sim_Est_BSB_Simp",
  output_name
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ssb_bias_file <- file.path(output_dir, "ssb_relative_bias_2022_2024.csv")
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
if (om_disable_movement || om_change_naa_sd) {
  om_input <- operating_model$input
  om_input$par <- operating_model$parList
  if (om_disable_movement) {
    om_input$data$can_move[] <- 0
    om_input$data$must_move[] <- 0
    om_input$data$use_mu_prior[] <- 0
    om_input$data$mu_model[] <- 1
    movement_parameters <- c("trans_mu", "mu_prior_re", "mu_re", "mu_repars")
    for (parameter in movement_parameters) {
      om_input$par[[parameter]][] <- 0
      om_input$map[[parameter]][] <- NA
      om_input$map[[parameter]] <- factor(om_input$map[[parameter]])
    }
    om_input$random <- setdiff(om_input$random, movement_parameters)
  }

  if (om_change_naa_sd) {
    om_input$par$log_NAA_sigma[] <- log(om_older_ages_sd)
    om_input$par$log_NAA_sigma[,, 1] <- log(om_recruitment_sd)
    om_input$map$log_NAA_sigma[] <- NA
    om_input$map$log_NAA_sigma <- factor(om_input$map$log_NAA_sigma)
  }
  operating_model <- wham::fit_wham(
    om_input,
    do.fit = FALSE,
    do.brps = FALSE,
    MakeADFun.silent = TRUE
  )
}
random_effects <- operating_model$input$random

estimation_model_config <- readRDS(here("models", "BSB.EM.Y.Config.RDS"))
em_options <- list(
  separate.em = FALSE,
  separate.em.type = NULL,
  do.move = FALSE,
  est.move = FALSE
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
  simulated_model,
  projected_model,
  projection_years
) {
  simulated_ssb <- simulated_model$rep$SSB
  projected_ssb <- projected_model$rep$SSB

  simulated_rows <- match(projection_years, simulated_model$years)
  projected_rows <- match(projection_years, projected_model$years_full)
  if (anyNA(simulated_rows) || anyNA(projected_rows)) {
    stop(
      "The projection years were not found in both simulated and projected model objects.",
      call. = FALSE
    )
  }

  tibble(
    region = c("North", "South"),
    mean_relative_error = colMeans(
      (projected_ssb[projected_rows, , drop = FALSE] -
        simulated_ssb[simulated_rows, , drop = FALSE]) /
        simulated_ssb[simulated_rows, , drop = FALSE],
      na.rm = TRUE
    )
  )
}

replicate_error_file <- file.path(output_dir, "replicate_errors.csv")
record_replicate_error <- function(replicate, seed, error) {
  error_record <- tibble(
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    simulation = replicate,
    seed = seed,
    error_class = paste(class(error), collapse = ";"),
    error_message = conditionMessage(error),
    call = if (is.null(conditionCall(error))) {
      ""
    } else {
      paste(deparse(conditionCall(error)), collapse = " ")
    }
  )
  write.table(
    error_record,
    file = replicate_error_file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(replicate_error_file),
    append = file.exists(replicate_error_file),
    quote = TRUE,
    qmethod = "double"
  )
}

run_replicate <- function(i) {
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

  # Keep all simulated input data, including observations, error values, and the
  # structural arrays needed to interpret the simulated outputs. The full WHAM
  # object is not saved because it contains unnecessary TMB internals.
  simulated_data <- list(
    years = full_years,
    data = simulated_model$input$data
  )

  # ==============================================================================
  # 3) Create and fit the reduced estimation model
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
  message(
    "Replicate ",
    i,
    ": convergence training=",
    training_converged
  )

  convergence <- tibble(
    model = "Training fit: 1989-2021",
    convergence_code = training_model$opt$convergence,
    converged = training_converged,
    pd_hessian = !isTRUE(training_model$na_sdrep)
  )
  if (training_converged) {
    converged_reps <<- c(converged_reps, i)
  } else {
    message(
      "Replicate ",
      i,
      ": training model did not converge; skipping remaining work"
    )
  }

  if (training_converged) {
    # ==============================================================================
    # 4) Project the 1989-2021 fit through 2024
    # ==============================================================================
    # Continue NAA/R random effects and apply F40% in 2022-2024 (see projection_options).
    # Compare this forecast with the simulated SSB for the projection years.

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

    if (save_simulation_outputs) {
      # Save only data and report arrays. This preserves all simulated data and
      # model outputs while excluding large TMB environments, functions, and
      # sdrep objects.
      simulation_outputs <- list(
        simulation = i,
        seed = simulation_seed,
        years = list(
          full = full_years,
          training = training_years,
          projection = projection_years
        ),
        simulated_data = simulated_data,
        simulated_outputs = simulated_model$rep,
        training_model = list(
          years = training_model$years,
          rep = training_model$rep
        ),
        projected_models = list(
          option1 = list(
            years = projected_model1$years,
            rep = projected_model1$rep
          ),
          option2 = list(
            years = projected_model2$years,
            rep = projected_model2$rep
          ),
          option3 = list(
            years = projected_model3$years,
            rep = projected_model3$rep
          )
        )
      )
      saveRDS(
        simulation_outputs,
        file.path(output_dir, paste0("simulation_outputs.sim", i, ".RDS")),
        compress = "xz"
      )
    }

    replicate_bias <- bind_rows(
      extract_projection_ssb_comparison(
        simulated_model,
        projected_model1,
        projection_years
      ) %>%
        mutate(projection_model = "Option 1"),
      extract_projection_ssb_comparison(
        simulated_model,
        projected_model2,
        projection_years
      ) %>%
        mutate(projection_model = "Option 2"),
      extract_projection_ssb_comparison(
        simulated_model,
        projected_model3,
        projection_years
      ) %>%
        mutate(projection_model = "Option 3")
    ) %>%
      mutate(simulation = i, .before = 1)
    ssb_bias_results[[length(ssb_bias_results) + 1L]] <<- replicate_bias
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

for (i in replicates_to_run) {
  simulation_seed <- r.seed.set[[1]][i]
  tryCatch(
    run_replicate(i),
    error = function(error) {
      tryCatch(
        record_replicate_error(i, simulation_seed, error),
        error = function(log_error) {
          message(
            "Could not write the error log for replicate ",
            i,
            ": ",
            conditionMessage(log_error)
          )
        }
      )
      message(
        "Replicate ",
        i,
        " failed: ",
        conditionMessage(error),
        ". Continuing with the next replicate."
      )
    }
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

# Read the saved results before plotting.
ssb_relative_bias <- read_csv(
  ssb_bias_file,
  show_col_types = FALSE
)

# Compare mean terminal-horizon SSB relative bias across projection options and
# regions. The zero line represents no relative bias.
label_offset <- max(
  diff(range(ssb_relative_bias$mean_relative_error, na.rm = TRUE)) * 0.04,
  0.01
)
mean_labels <- ssb_relative_bias |>
  group_by(projection_model, region) |>
  summarise(
    mean_value = mean(mean_relative_error, na.rm = TRUE),
    min_value = min(mean_relative_error, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    label_y = min_value - label_offset,
    label = sprintf("Mean: %.3f", mean_value)
  )

ggplot(
  ssb_relative_bias,
  aes(x = projection_model, y = mean_relative_error, fill = region)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA
  ) +
  geom_point(
    aes(group = region),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    color = "black",
    alpha = 0.45,
    size = 1.2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(
    data = mean_labels,
    aes(
      x = projection_model,
      y = label_y,
      label = label,
      group = region
    ),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.05))) +
  labs(
    x = "Projection model",
    y = "Mean SSB relative bias",
    fill = "Region"
  )
ggsave(
  file.path(output_dir, "ssb_relative_bias_comparison.png"),
  width = 8,
  height = 5,
  units = "in",
  dpi = 300
)

message(
  "Simulation-estimation comparison complete. Results saved in: ",
  output_dir
)
