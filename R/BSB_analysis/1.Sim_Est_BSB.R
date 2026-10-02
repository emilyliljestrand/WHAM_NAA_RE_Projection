#' @title Black Sea Bass Simulation-Estimation Projection Test
#' @description Simulates 1989-2024 catch and index observations from BSB.EM.Y,
#'   fits full (1989-2024) and truncated (1989-2021) estimation models, projects
#'   the truncated model through 2024, and compares projections with the full fit.
#' @name sim_est_bsb

# WHAM is loaded from a pinned project library so the BSB model stays on
# version 2.1.0.9003. SPASAM.MSE supplies update_om_fn() and make_em_input().
rm(list = ls())

suppressPackageStartupMessages({
  library(wham, lib.loc = "C:/Users/emily.liljestrand/AppData/Local/R/win-library/4.4/wham_2.1.0.9003")
  library(SPASAM.MSE)
  library(tidyverse)
  library(here)
})

#Read in existing seeds file, or make one if it doesn't exist
if(!file.exists("data/raw/seeds/Sim_Est_BSB_seeds.csv"))
{
  set.seed(799291)
  nseeds <- 5000
  r.seed.set <- trunc(1e7*runif(n=nseeds),7) + trunc(1e3*runif(n=nseeds),3)
  write.table(r.seed.set, file="data/raw/seeds/Sim_Est_BSB_seeds.csv", quote=F, row.names=F, col.names=F, sep=",")
} else r.seed.set <- read.table("data/raw/seeds/Sim_Est_BSB_seeds.csv")
# ==============================================================================
# 1) User controls
# ==============================================================================
# Full OM/EM window is 1989-2024. The training EM stops in 2021 and is then
# projected through 2022-2024 for comparison with the full fit.

nreps <- 10
convergence_rate <- c()

for(i in 1:nreps)
{
  simulation_seed <- r.seed.set[[1]][i]
  full_years <- 1989:2024
  training_years <- 1989:2021
  projection_years <- 2022:2024
  
  # Continue recruitment and NAA random effects (opt = 1) and set each projected
  # year to F40% (proj_F_opt = 3).
  projection_options1 <- list(
    proj_R_opt = 1,
    proj_NAA_opt = 1,
    proj_F_opt = rep(3, length(projection_years))
  )
  # Turn off random effects on both Recruitment and NAA and set each projected
  # year to F40% (proj_F_opt = 3).
  projection_options2 <- list(
    proj_R_opt = 4,
    proj_NAA_opt = 3,
    proj_F_opt = rep(3, length(projection_years))
  )
  # Average recruitment and NAA random effects over historical reference period
  projection_options3 <- list(
    proj_R_opt = 3,
    proj_NAA_opt = 2,
    proj_F_opt = rep(3, length(projection_years))
  )
  
  output_dir <- here("output", "BSB_analysis", "1.Sim_Est_BSB")
  model_dir <- file.path(output_dir, "models")
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  
  stopifnot(
    identical(training_years, full_years[full_years <= max(training_years)]),
    length(projection_years) == 3L,
    min(projection_years) == max(training_years) + 1L,
    file.exists(here("models", "BSB.EM.Y.RDS"))
  )
  
  # Summarize whether each EM converged and whether the Hessian was positive
  # definite (needed for sdreport / uncertainty).
  check_convergence <- function(model) {
    tibble(
      convergence_code = model$opt$convergence,
      converged = identical(model$opt$convergence, 0L),
      pd_hessian = !isTRUE(model$na_sdrep)
    )
  }
  
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
  
  saveRDS(simulated_model, file.path(model_dir, paste("BSB.simulated.1989.2024.sim",i,".RDS",sep="")))
  
  # Also save just the simulated catch/index observations for inspection without
  # loading the full WHAM object.
  simulated_data <- list(
    years = full_years,
    catch = simulated_model$input$data[c(
      "agg_catch", "catch_paa", "use_agg_catch", "use_catch_paa",
      "agg_catch_sigma", "catch_Neff"
    )],
    index = simulated_model$input$data[c(
      "agg_index", "index_paa", "use_indices", "use_index_paa",
      "agg_index_sigma", "index_Neff"
    )]
  )
  saveRDS(simulated_data, file.path(output_dir, paste("simulated_catch_index_1989_2024.sim",i,".RDS",sep="")))
  
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
  fracyr_spawn <- om_data$fracyr_SSB[1, 1] + season_starts[om_data$spawn_seasons[1]]
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
  om_F <- t(operating_model$rep$FAA[1:4,,8])
  
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
    F_info = list(F.year1 = om_F[1, 1],
                   Fhist = "constant",
                   Fmax = max(om_F),
                   Fmin = min(om_F),
                   change_time = 0.5,
                   user_F = om_F,
                   F_feedback = NULL),
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
  saveRDS(full_model, file.path(model_dir, paste("BSB.simulated.EM.1989.2024.sim",i,".RDS",sep="")))
  
  
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
  saveRDS(training_model, file.path(model_dir, paste("BSB.simulated.EM.1989.2021.sim",i,".RDS",sep="")))
  
  convergence <- bind_rows(
    check_convergence(full_model) %>% mutate(model = "Full fit: 1989-2024", .before = 1),
    check_convergence(training_model) %>% mutate(model = "Training fit: 1989-2021", .before = 1)
  )
  write_csv(convergence, file.path(output_dir, paste("model_convergence.sim",i,".csv",sep="")))
  
  if(all(convergence$converged)) convergence_rate <- c(convergence_rate,i)
  
  if(all(convergence$converged))
  {
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
  saveRDS(projected_model1, file.path(model_dir, paste("BSB.simulated.EM.1989.2021.Proj.2022.2024.Opt1.sim",i,".RDS",sep="")))
  saveRDS(projected_model1, file.path(model_dir, paste("BSB.simulated.EM.1989.2021.Proj.2022.2024.Opt2.sim",i,".RDS",sep="")))
  saveRDS(projected_model1, file.path(model_dir, paste("BSB.simulated.EM.1989.2021.Proj.2022.2024.Opt3.sim",i,".RDS",sep="")))
  }
}

# ==============================================================================
# 5) Compare projected values with the full 1989-2024 fit
# ==============================================================================
# Differences and relative differences for 2022-2024 (years held out of the
# training EM) for SSB, biomass, recruitment, and Fbar.

# Pull one reported quantity from both models and reshape to year x component.
extract_projection_comparison <- function(quantity, labels = NULL) {
  full_values <- full_model$rep[[quantity]]
  projected_values <- projected_model$rep[[quantity]]

  if (is.null(dim(full_values))) {
    full_values <- matrix(full_values, ncol = 1L)
    projected_values <- matrix(projected_values, ncol = 1L)
  }

  full_rows <- match(projection_years, full_model$years)
  projected_rows <- match(projection_years, projected_model$years)
  if (anyNA(full_rows) || anyNA(projected_rows)) {
    stop("Projection years were not found in both fitted model objects.", call. = FALSE)
  }

  component_labels <- labels %||% paste0("Component_", seq_len(ncol(full_values)))
  tidyr::expand_grid(year = projection_years, component = component_labels) %>%
    mutate(
      quantity = quantity,
      full_fit = as.vector(t(full_values[full_rows, , drop = FALSE])),
      projection = as.vector(t(projected_values[projected_rows, , drop = FALSE])),
      difference = projection - full_fit,
      relative_difference = difference / full_fit
    ) %>%
    select(quantity, component, year, everything())
}

comparison <- bind_rows(
  extract_projection_comparison("SSB", c("North", "South")),
  extract_projection_comparison("B", c("North", "South")),
  extract_projection_comparison("R", c("North", "South")),
  extract_projection_comparison("Fbar", operating_model_fit$input$catch_info$fleet_names)
)

write_csv(comparison, file.path(output_dir, "projection_comparison_2022_2024.csv"))

# SSB time series: full EM vs projection, by region.
comparison_plot <- comparison %>%
  filter(quantity == "SSB") %>%
  select(quantity, component, year, full_fit, projection) %>%
  pivot_longer(c(full_fit, projection), names_to = "model", values_to = "value") %>%
  ggplot(aes(year, value, color = model, linetype = model)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(vars(component), scales = "free_y") +
  scale_color_manual(values = c(full_fit = "black", projection = "#0072B2")) +
  labs(x = NULL, y = "Spawning stock biomass", color = NULL, linetype = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_dir, "ssb_projection_comparison_2022_2024.png"),
  comparison_plot,
  width = 8,
  height = 4.5,
  dpi = 300
)

# Standard WHAM comparison plots (SSB, F, recruitment, etc.). AIC is not
# meaningful here because the two models are fit to different year ranges.
compare_wham_models(
  list(Full_1989_2024 = full_model, Projected_from_2021 = projected_model),
  fdir = output_dir,
  calc.aic = FALSE,
  do.table = FALSE,
  plot.opts = list(which = c(1, 6, 7, 8, 9, 10), kobe.yr = 2021)
)

message("Simulation-estimation comparison complete. Results saved in: ", output_dir)