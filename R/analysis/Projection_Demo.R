# ==============================================================================
# BSB whamMSE test script (OM with fixed ecov) + EM convergence test
# + optional full MSE loop with 3 projection options (proj_NAA_opt = 1,2,3)
#
# ==============================================================================

suppressPackageStartupMessages({
  library(wham)
  library(whamMSE)
  library(dplyr)
})

## =============================================================================
## 0) USER CONTROLS (edit here)
## =============================================================================

# ---- Paths
# proj_dir <- "C:/Users/liche/Desktop/Rutgers-MSE"
proj_dir <- getwd()

# ---- Time
year_start <- 1989
year_end   <- 2024

# ---- MSE years (feedback/projection years appended to time-varying inputs)
n_feedback_years <- 3

# ---- Dimensions
n_stocks  <- 2
n_regions <- 2
n_ages    <- 8

# ---- Seeds
seed_main <- 1

# ---- Run controls
do_quick_plots <- TRUE

# (2) EM convergence test vs full MSE loop
EM_fit   <- TRUE    # default TRUE: fit one EM once and check convergence
full_MSE <- FALSE   # default FALSE: do NOT run full MSE unless TRUE

# (3) NAA RE sigma controls
Rec_sigma <- 0.5   # age-1 recruitment sigma
NAA_sigma <- 0.5   # ages >=2 sigma

# (3) Toggle environmental covariate link to recruitment
use_ecov_rec_link <- TRUE   # FALSE => "none" everywhere

# Full MSE: projection options you want to compare
proj_opts_to_run <- c(1, 2, 3)

## =============================================================================
## 1) Setup utilities (keep minimal)
## =============================================================================

setwd(proj_dir)
p <- function(...) file.path(proj_dir, ...)

stopifnot(dir.exists(proj_dir))
stopifnot(file.exists(p("models", "fit.RDS")))
stopifnot(file.exists(p("data", "NORTH.1989.2024.dat")))
stopifnot(file.exists(p("data", "SOUTH.1989.2024.dat")))
stopifnot(file.exists(p("data", "bsb_bt_temp_nmab_1959-2024.csv")))
stopifnot(file.exists(p("data", "bsb_bt_temp_smab_1959-2024.csv")))

# Inverse logit helper for selectivity initialization
gen.invlogit <- function(eta, low, upp, s = 1) {
  low + (upp - low) * plogis(s * eta)
}

# Convergence checker (used both for EM test and loop)
check_conv <- function(em) {
  conv   <- isTRUE(as.logical(1 - em$opt$convergence))
  pdHess <- isTRUE(!isTRUE(em$na_sdrep) && !is.na(em$na_sdrep))
  list(conv = conv, pdHess = pdHess)
}

## =============================================================================
## 2) Load fitted OM + ASAP + ecov time series
## =============================================================================

OMa <- readRDS(p("models", "fit.RDS"))

asap <- wham::read_asap3_dat(p("data", c("NORTH.1989.2024.dat", "SOUTH.1989.2024.dat")))
temp <- wham::prepare_wham_input(asap)

north_bt <- read.csv(p("data", "bsb_bt_temp_nmab_1959-2024.csv"))
south_bt <- read.csv(p("data", "bsb_bt_temp_smab_1959-2024.csv"))

## =============================================================================
## 3) Build ecov object (with ON/OFF rec linkage)
## =============================================================================

ecov <- list(label = c("North_BT", "South_BT"))

# Historical deviations (center by mean), append zeros for future years
ecov$mean <- cbind(north_bt[, "mean"], south_bt[, "mean"])
ecov$mean <- t(t(ecov$mean) - colMeans(ecov$mean))
ecov$mean <- rbind(ecov$mean, matrix(0, n_feedback_years, ncol(ecov$mean)))

# Obs SD (constant here)
ecov$logsigma <- matrix(log(0.2), nrow(ecov$mean), ncol(ecov$mean))

# Years (append)
ecov$year <- c(
  north_bt[, "year"],
  (max(north_bt[, "year"]) + 1):(max(north_bt[, "year"]) + n_feedback_years)
)

ecov$use_obs <- matrix(1, nrow(ecov$mean), ncol(ecov$mean))
ecov$process_model <- "ar1"
ecov$process_mean_vals <- colMeans(ecov$mean)

# Recruitment linkage toggle
if (isTRUE(use_ecov_rec_link)) {
  ecov$recruitment_how <- matrix(
    c("controlling-lag-0-linear", "none",
      "none", "controlling-lag-0-linear"),
    nrow = 2, ncol = 2
  )
} else {
  ecov$recruitment_how <- matrix("none", nrow = 2, ncol = 2)
}

if (do_quick_plots) {
  plot(ecov$mean[, 1], type = "l", main = "Ecov (North) deviations")
  plot(ecov$mean[, 2], type = "l", main = "Ecov (South) deviations")
}

## =============================================================================
## 4) Build maturity + WAA over historical + projection years
## =============================================================================

hist_years  <- length(year_start:year_end)
total_years <- length(year_start:(year_end + n_feedback_years))

# maturity: stock x year x age
user_maturity <- array(NA, dim = c(n_stocks, total_years, n_ages))
user_maturity[, 1:hist_years, ] <- OMa$input$data$mature
for (i in (hist_years + 1):(hist_years + n_feedback_years)) {
  user_maturity[, i, ] <- OMa$input$data$mature[, hist_years, , drop = FALSE]
}

# Ensure M uses same pointer as SSB WAA
OMa$input$data$waa_pointer_M <- OMa$input$data$waa_pointer_ssb

# WAA: source x year x age
user_waa <- list()
user_waa$waa <- array(NA, dim = c(10, hist_years + n_feedback_years, n_ages))
user_waa$waa[, 1:hist_years, ] <- OMa$input$data$waa
for (i in (hist_years + 1):(hist_years + n_feedback_years)) {
  user_waa$waa[, i, ] <- OMa$input$data$waa[, hist_years, , drop = FALSE]
}

user_waa$waa_pointer_fleets   <- OMa$input$data$waa_pointer_fleets
user_waa$waa_pointer_indices  <- OMa$input$data$waa_pointer_indices
user_waa$waa_pointer_totcatch <- OMa$input$data$waa_pointer_ssb
user_waa$waa_pointer_ssb      <- OMa$input$data$waa_pointer_ssb
user_waa$waa_pointer_M        <- OMa$input$data$waa_pointer_M

fracyr_spawn   <- asap[[1]]$dat$fracyr_spawn
fracyr_seasons <- OMa$input$data$fracyr_seasons  # use OM seasons

## =============================================================================
## 5) Catch + index configs (simple defaults)
## =============================================================================

catch_info <- list(
  catch_cv      = c(0.05, 0.15, 0.05, 0.15),
  catch_Neff    = c(50, 50, 50, 50),
  use_agg_catch = 1,
  use_catch_paa = 1
)

# NOTE: if you want to keep your old reconcile logic, you can swap these
# back to your timing_fix workflow. Here I keep it simple as you have.
fracyr_indices <- c(0.5, 0.25, 0.5, 0.25)

index_info <- list(
  index_cv        = rep(0.4, 4),
  index_Neff      = rep(25, 4),
  fracyr_indices  = fracyr_indices,
  q               = OMa$rep$q[1, ],
  use_indices     = rep(1, 4),
  use_index_paa   = rep(1, 4),
  units_indices   = rep(2, 4),
  units_index_paa = rep(2, 4)
)

## =============================================================================
## 6) Generate basic_info, then set NAA_where manually (SIMPLE)
## =============================================================================

info <- whamMSE::generate_basic_info(
  n_stocks         = n_stocks,
  n_regions        = n_regions,
  n_indices        = 4,
  n_fleets         = 4,
  n_seasons        = 11,
  base.years       = year_start:year_end,
  n_feedback_years = n_feedback_years,
  n_ages           = n_ages,
  catch_info       = catch_info,
  index_info       = index_info,
  user_waa         = user_waa,
  user_maturity    = user_maturity,
  fracyr_spawn     = fracyr_spawn,
  fracyr_seasons   = fracyr_seasons
)

basic_info     <- info$basic_info
catch_info_use <- info$catch_info
index_info_use <- info$index_info

# ---- Manual NAA_where rules
basic_info$NAA_where <- array(1, dim = c(n_stocks, n_regions, n_ages))

# Age-1 diagonal only:
basic_info$NAA_where[, , 1] <- matrix(
  c(1, 0,
    0, 1),
  nrow = n_stocks, ncol = n_regions, byrow = TRUE
)

# Stock 2 never in region 1 at any age:
basic_info$NAA_where[2, 1, ] <- 0

info$basic_info <- basic_info

cat("\nNAA_where at age 1 (rows=stock, cols=region):\n")
print(info$basic_info$NAA_where[, , 1])

## =============================================================================
## 7) F setup
## =============================================================================

F_info <- info$F
F_info$F[1:hist_years, ] <- OMa$rep$Fbar[, 1:4]

## =============================================================================
## 8) Selectivity setup
## =============================================================================

sel <- list(
  n_selblocks = 8,
  model       = rep(c("age-specific", "logistic", "age-specific"), c(2, 2, 4)),
  re          = rep("none", 8)
)

sel_vals <- gen.invlogit(OMa$parList$logit_selpars, low = 0, upp = 1)

sel$initial_pars <- list(
  sel_vals[1, 1:8],  # north comm
  sel_vals[2, 1:8],  # north rec
  sel_vals[3, 9:10], # south comm
  sel_vals[4, 9:10], # south rec
  sel_vals[5, 1:8],  # north rec cpa
  sel_vals[6, 1:8],  # north vast
  sel_vals[7, 1:8],  # south rec cpa
  sel_vals[8, 1:8]   # south vast
)

sel$fix_pars <- list(
  4:8,   # north comm
  7:8,   # north rec
  NULL,  # south comm
  NULL,  # south rec
  2:8,   # north rec cpa
  5:8,   # north vast
  3:8,   # south rec cpa
  2:8    # south vast
)

## =============================================================================
## 9) Natural mortality (M)
## =============================================================================

M <- list(
  model = "constant",
  initial_means = array(0.4, dim = c(n_stocks, n_regions, n_ages))
)

## =============================================================================
## 10) NAA random effects (sigma controls at top)
## =============================================================================

sigma_vals <- array(NAA_sigma, dim = c(n_stocks, n_regions, n_ages))
sigma_vals[, , 1] <- Rec_sigma

NAA_re <- list(
  recruit_model = 2,
  recruit_pars  = list(exp(9.032686), exp(9.753513)),
  sigma_vals    = sigma_vals,
  sigma         = list("rec+1", "rec+1"),
  cor           = list("2dar1", "2dar1"),
  N1_model      = rep("equilibrium", 2)
)

## =============================================================================
## 11) Prepare WHAM input for OM (fixed ecov)
## =============================================================================

input_Ecov <- wham::prepare_wham_input(
  basic_info   = info$basic_info,
  selectivity  = sel,
  M            = M,
  NAA_re       = NAA_re,
  ecov         = ecov,
  catch_info   = catch_info_use,
  index_info   = index_info_use,
  F            = F_info
)

# Update WAA pointers (important)
input_Ecov <- update_waa(input_Ecov, waa_info = info$par_inputs$user_waa)

# Initialize ecov process parameters and effects from fitted OM
input_Ecov$par$Ecov_process_pars <- OMa$parList$Ecov_process_pars
input_Ecov$par$Ecov_beta_R       <- OMa$parList$Ecov_beta_R

# Initialize N1 from fitted OM
input_Ecov$par$log_N1 <- OMa$parList$log_N1

# ---- Index sigma/usage: copy hist; repeat last historical year for future
input_Ecov$data$agg_index_sigma[1:hist_years, ] <- OMa$input$data$agg_index_sigma
input_Ecov$data$use_indices[1:hist_years, ]     <- OMa$input$data$use_indices
input_Ecov$data$use_index_paa[1:hist_years, ]   <- OMa$input$data$use_index_paa
for (i in (hist_years + 1):(hist_years + n_feedback_years)) {
  input_Ecov$data$agg_index_sigma[i, ] <- OMa$input$data$agg_index_sigma[hist_years, ]
}

# ---- index_Neff: build from ASAP then append future rows explicitly
idx1 <- which(asap[[1]]$dat$use_index == 1)
idx2 <- which(asap[[2]]$dat$use_index == 1)
Neff1 <- do.call(cbind, lapply(idx1, function(i) asap[[1]]$dat$IAA_mats[[i]][, 12, drop = FALSE]))
Neff2 <- do.call(cbind, lapply(idx2, function(i) asap[[2]]$dat$IAA_mats[[i]][, 12, drop = FALSE]))
index_Neff <- cbind(Neff1, Neff2)
index_Neff <- rbind(index_Neff, index_Neff[rep(hist_years, n_feedback_years), , drop = FALSE])
input_Ecov$data$index_Neff <- index_Neff

input_Ecov <- whamMSE::update_input_index_info(
  input_Ecov,
  agg_index_sigma = input_Ecov$data$agg_index_sigma,
  index_Neff      = input_Ecov$data$index_Neff
)

# ---- Catch sigma/usage: copy hist; repeat last historical year for future
input_Ecov$data$agg_catch_sigma[1:hist_years, ] <- OMa$input$data$agg_catch_sigma
input_Ecov$data$use_agg_catch[1:hist_years, ]   <- OMa$input$data$use_agg_catch
input_Ecov$data$use_catch_paa[1:hist_years, ]   <- OMa$input$data$use_catch_paa
for (i in (hist_years + 1):(hist_years + n_feedback_years)) {
  input_Ecov$data$agg_catch_sigma[i, ] <- OMa$input$data$agg_catch_sigma[hist_years, ]
}

# ---- catch_Neff: from ASAP then append future rows explicitly
catch_Neff <- cbind(asap[[1]]$dat$catch_Neff, asap[[2]]$dat$catch_Neff)
catch_Neff <- rbind(catch_Neff, catch_Neff[rep(hist_years, n_feedback_years), , drop = FALSE])
input_Ecov$data$catch_Neff <- catch_Neff

input_Ecov <- whamMSE::update_input_catch_info(
  input_Ecov,
  agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
  catch_Neff      = input_Ecov$data$catch_Neff
)

# ---- Freeze Ecov_re (NO hard-coded indices)
Ecov_re_hist <- OMa$parList$Ecov_re
hist_n <- nrow(Ecov_re_hist)
fut_n  <- n_feedback_years

stopifnot(ncol(Ecov_re_hist) == length(ecov$label))

# Ensure input_Ecov$par$Ecov_re is long enough
if (nrow(input_Ecov$par$Ecov_re) < (hist_n + fut_n)) {
  n_need <- (hist_n + fut_n) - nrow(input_Ecov$par$Ecov_re)
  input_Ecov$par$Ecov_re <- rbind(
    input_Ecov$par$Ecov_re,
    input_Ecov$par$Ecov_re[rep(nrow(input_Ecov$par$Ecov_re), n_need), , drop = FALSE]
  )
}

# Copy historical
input_Ecov$par$Ecov_re[1:hist_n, ] <- Ecov_re_hist

# Fill future: reuse first fut_n historical rows if available, else repeat last row
if (fut_n > 0) {
  if (hist_n >= fut_n) {
    input_Ecov$par$Ecov_re[(hist_n + 1):(hist_n + fut_n), ] <- Ecov_re_hist[1:fut_n, , drop = FALSE]
  } else {
    input_Ecov$par$Ecov_re[(hist_n + 1):(hist_n + fut_n), ] <- Ecov_re_hist[rep(hist_n, fut_n), , drop = FALSE]
  }
}

input_Ecov$data$do_simulate_Ecov_re <- 0

# Copy NAA AR1 params
input_Ecov$par$trans_NAA_rho <- OMa$parList$trans_NAA_rho

## =============================================================================
## 12) Build OM object (unfitted), remove Ecov_re from random, generate one realization
## =============================================================================

unfitted_om <- fit_wham(input_Ecov, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = FALSE)

if ("Ecov_re" %in% unfitted_om$input$random) {
  input_Ecov$random <- unfitted_om$input$random[unfitted_om$input$random != "Ecov_re"]
}

random <- input_Ecov$random
input_Ecov$random <- NULL

om_ecov <- fit_wham(input_Ecov, do.fit = FALSE, do.brps = TRUE, MakeADFun.silent = FALSE)

om_with_data <- update_om_fn(om_ecov, seed = seed_main, random = random)

if (do_quick_plots) {
  par(mfrow = c(2, 2))
  plot(om_with_data$rep$SSB[, 1], type = "l", main = "OM SSB North")
  plot(om_with_data$rep$SSB[, 2], type = "l", main = "OM SSB South")
  plot(om_with_data$rep$Ecov_x[, 1], type = "l", main = "OM Ecov_x North")
  plot(om_with_data$rep$Ecov_x[, 2], type = "l", main = "OM Ecov_x South")
  par(mfrow = c(1, 1))
}

## =============================================================================
## 13) EM configuration
## =============================================================================

# Movement
seasons <- c(rep(1, 5), 2, rep(1, 5)) / 12
n_seasons <- length(seasons)

move <- list(stock_move = c(TRUE, FALSE), separable = TRUE)
move$must_move <- array(0, dim = c(n_stocks, n_seasons, n_regions))
move$must_move[1, 5, 2] <- 1

move$can_move <- array(0, dim = c(n_stocks, n_seasons, n_regions, n_regions))
move$can_move[1, 1:4, 2, 1]  <- 1
move$can_move[1, 7:11, 1, 2] <- 1
move$can_move[1, 5, 2, ]     <- 1

mus <- array(0, dim = c(n_stocks, n_seasons, n_regions, 1))
mus[1, 1:11, 1, 1] <- 0.02214863
mus[1, 1:11, 2, 1] <- 0.3130358

move$mean_vals  <- mus
move$mean_model <- matrix("stock_constant", 2, 1)

move_em <- move
move_em$use_prior   <- array(0,   dim = c(n_stocks, n_seasons, n_regions, n_regions - 1))
move_em$use_prior[1, 1, , ] <- 1
move_em$prior_sigma <- array(0.2, dim = c(n_stocks, n_seasons, n_regions, n_regions - 1))

# HCR
hcr <- list(
  hcr.type = 1,
  hcr.opts = list(use_FXSPR = TRUE, percentFXSPR = 75)
)

# Assessment schedule
assess.interval <- 3
base.years      <- year_start:year_end
terminal.year   <- tail(base.years, 1)
last.year       <- max(om_with_data$years)
assess.years    <- seq(terminal.year, last.year - assess.interval, by = assess.interval)

# EM NAA RE (equilibrium helps)
NAA_re_em <- NAA_re
NAA_re_em$N1_model[] <- "equilibrium"

# EM ecov
ecov_em <- ecov

# EM options
em.opt <- list(
  separate.em      = FALSE,
  separate.em.type = NULL,
  do.move          = TRUE,
  est.move         = TRUE
)

## =============================================================================
## 14) EM convergence TEST (single fit)
## =============================================================================

em_test <- NULL
em_test_input <- NULL
conv_test <- NULL

if (isTRUE(EM_fit)) {
  
  cat("\n============================================================\n")
  cat("EM convergence TEST (single fit) at year:", terminal.year, "\n")
  cat("============================================================\n")
  
  em_years_test <- base.years[1]:terminal.year
  year.use_test <- length(base.years)
  
  update_catch_info <- list(
    agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
    catch_Neff      = input_Ecov$data$catch_Neff
  )
  update_index_info <- list(
    agg_index_sigma = input_Ecov$data$agg_index_sigma,
    index_Neff      = input_Ecov$data$index_Neff
  )
  
  em_test_input <- make_em_input(
    om = om_with_data,
    em_info = info,
    M_em = M,
    sel_em = sel,
    NAA_re_em = NAA_re_em,
    move_em = move_em,
    catchability_em = NULL,
    ecov_em = ecov_em,
    em.opt = em.opt,
    em_years = em_years_test,
    year.use = year.use_test,
    age_comp_em = "multinomial",
    aggregate_catch_info = NULL,
    aggregate_index_info = NULL,
    aggregate_weights_info = NULL,
    filter_indices = NULL,
    reduce_region_info = NULL,
    update_catch_info = update_catch_info,
    update_index_info = update_index_info
    # ecov_em_opts = NULL
  )
  
  if (isTRUE(em.opt$do.move) && !isTRUE(em.opt$est.move)) {
    em_test_input <- fix_move(em_test_input)
  }
  
  cat("\n[Step] fit_wham() for EM test\n")
  em_test <- fit_wham(
    em_test_input,
    do.retro = FALSE,
    do.osa   = FALSE,
    do.brps  = TRUE,
    MakeADFun.silent = FALSE
  )
  
  conv_test <- check_conv(em_test)
  cat("\nEM TEST results:\n")
  cat("  Converged:", conv_test$conv, "\n")
  cat("  PD Hessian:", conv_test$pdHess, "\n")
  
  if (!conv_test$conv || !conv_test$pdHess) {
    warning("EM TEST did not converge and/or Hessian not PD. Consider adjusting N1_model, priors, sigma, etc.")
  }
}

# If full_MSE is FALSE, stop after EM test
if (!isTRUE(full_MSE)) {
  cat("\n============================================================\n")
  cat("full_MSE = FALSE => stopping after EM convergence test.\n")
  cat("Set full_MSE <- TRUE at the top to run full MSE loops.\n")
  cat("============================================================\n")
  
  out_test <- list(
    settings = list(
      year_start = year_start,
      year_end   = year_end,
      n_feedback_years = n_feedback_years,
      Rec_sigma = Rec_sigma,
      NAA_sigma = NAA_sigma,
      use_ecov_rec_link = use_ecov_rec_link,
      EM_fit   = EM_fit,
      full_MSE = full_MSE
    ),
    om_with_data = om_with_data,
    random = random,
    em_test_input = em_test_input,
    em_test = em_test,
    em_test_convergence = conv_test
  )
  
  # saveRDS(out_test, file = p("models", "em_convergence_test_out.rds"))
  # str(out_test, max.level = 2)
}

## =============================================================================
## 15) FULL MSE LOOP (runs only if full_MSE == TRUE) for proj_NAA_opt = 1,2,3
## =============================================================================
full_MSE = TRUE

out_mse <- vector("list", length(proj_opts_to_run))
names(out_mse) <- paste0("proj_NAA_opt_", proj_opts_to_run)

if (isTRUE(full_MSE)) {
  
  for (k in seq_along(proj_opts_to_run)) {
    
    proj_opt_k <- proj_opts_to_run[k]
    proj.opts <- list(proj_NAA_opt = proj_opt_k)
    
    cat("\n============================================================\n")
    cat("FULL MSE LOOP: proj_NAA_opt =", proj_opt_k, "\n")
    cat("============================================================\n")
    
    # Reset OM at the beginning of each scenario so results are comparable
    # (same starting OM realization; differences arise only from proj_NAA_opt)
    om <- om_with_data
    
    assess_years_use <- assess.years
    year.use <- length(base.years)
    add.years <- TRUE
    
    seed <- seed_main
    
    save.last.em <- TRUE
    save.sdrep   <- FALSE
    do.retro     <- FALSE
    do.osa       <- FALSE
    do.brps      <- FALSE
    by_fleet     <- TRUE
    FXSPR_init   <- NULL
    
    implementation_error <- NULL
    
    update_catch_info <- list(
      agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
      catch_Neff      = input_Ecov$data$catch_Neff
    )
    update_index_info <- list(
      agg_index_sigma = input_Ecov$data$agg_index_sigma,
      index_Neff      = input_Ecov$data$index_Neff
    )
    
    start.time <- Sys.time()
    
    em_list        <- vector("list", length(assess_years_use))
    par.est        <- vector("list", length(assess_years_use))
    par.se         <- vector("list", length(assess_years_use))
    adrep.est      <- vector("list", length(assess_years_use))
    adrep.se       <- vector("list", length(assess_years_use))
    opt_list       <- vector("list", length(assess_years_use))
    converge_list  <- vector("list", length(assess_years_use))
    catch_advice   <- vector("list", length(assess_years_use))
    catch_realized <- vector("list", length(assess_years_use))
    em_full        <- vector("list", length(assess_years_use))
    em_input_list  <- vector("list", length(assess_years_use))
    
    for (y in assess_years_use) {
      
      cat("\n============================================================\n")
      cat("Assessment year:", y, " | proj_NAA_opt:", proj_opt_k, "\n")
      cat("============================================================\n")
      
      i <- which(assess_years_use == y)
      
      em.years <- base.years[1]:y
      if (add.years && i != 1) {
        year.use <- year.use + assess.interval
      }
      
      cat("\n[Step] make_em_input()\n")
      em_input <- make_em_input(
        om = om,
        em_info = info,
        M_em = M,
        sel_em = sel,
        NAA_re_em = NAA_re_em,
        move_em = move_em,
        catchability_em = NULL,
        ecov_em = ecov_em,
        em.opt = em.opt,
        em_years = em.years,
        year.use = year.use,
        age_comp_em = "multinomial",
        aggregate_catch_info = NULL,
        aggregate_index_info = NULL,
        aggregate_weights_info = NULL,
        filter_indices = NULL,
        reduce_region_info = NULL,
        update_catch_info = update_catch_info,
        update_index_info = update_index_info,
        ecov_em_opts = NULL
      )
      
      if (!is.null(FXSPR_init)) {
        cat("\n[Step] set FXSPR_init\n")
        em_input$data$FXSPR_init[] <- FXSPR_init
      }
      
      cat("\n[Step] fit_wham()\n")
      if (isTRUE(em.opt$do.move)) {
        if (isTRUE(em.opt$est.move)) {
          em <- fit_wham(em_input, do.retro = do.retro, do.osa = do.osa,
                         do.brps = TRUE, MakeADFun.silent = FALSE)
        } else {
          em_input_fixed <- fix_move(em_input)
          em <- fit_wham(em_input_fixed, do.retro = do.retro, do.osa = do.osa,
                         do.brps = TRUE, MakeADFun.silent = FALSE)
          em_input <- em_input_fixed
        }
      } else {
        em <- fit_wham(em_input, do.retro = do.retro, do.osa = do.osa,
                       do.brps = TRUE, MakeADFun.silent = FALSE)
      }
      
      cat("\n[Step] convergence\n")
      cc <- check_conv(em)
      cat("Converged:", cc$conv, " | PD Hessian:", cc$pdHess, "\n")
      
      # ---- Advice + OM update
      if (assess.interval != 0) {
        
        cat("\n[Step] advice_fn() (projection)\n")
        advice <- advice_fn(
          em = em,
          pro.yr = assess.interval,
          hcr = hcr,
          proj.opts = proj.opts
        )
        
        # Ensure advice matrix [assess.interval x n_fleets]
        if (is.vector(advice)) {
          if (assess.interval == 1) {
            advice <- as.matrix(t(advice))
          } else {
            advice <- matrix(advice, byrow = TRUE, nrow = assess.interval)
          }
        }
        
        colnames(advice) <- paste0("Fleet_", 1:om$input$data$n_fleets)
        rownames(advice) <- paste0("Year_", y + seq_len(assess.interval))
        
        cat("\nCatch Advice:\n")
        print(advice)
        
        if (!is.null(implementation_error)) {
          cat("\n[Step] add_implementation_error()\n")
          real_catch <- add_implementation_error(
            catch_advice   = advice,
            method         = implementation_error$method,
            mean           = implementation_error$mean,
            cv             = implementation_error$cv,
            sd             = implementation_error$sd,
            min            = implementation_error$min,
            max            = implementation_error$max,
            constant_value = implementation_error$constant_value,
            seed           = seed
          )
        } else {
          real_catch <- advice
        }
        
        cat("\nRealized catch:\n")
        print(real_catch)
        
        interval.info <- list(
          catch = real_catch,
          years = y + seq_len(assess.interval)
        )
        
        cat("\n[Step] update_om_fn()\n")
        om <- update_om_fn(
          om,
          interval.info,
          seed   = seed,
          random = random,
          method = "nlminb",
          by_fleet = by_fleet,
          do.brps  = do.brps
        )
        
      } else {
        advice <- NULL
        real_catch <- NULL
      }
      
      # ---- Save
      cat("\n[Step] Save outputs\n")
      em_list[[i]]        <- em$rep
      par.est[[i]]        <- as.list(em$sdrep, "Estimate")
      par.se[[i]]         <- as.list(em$sdrep, "Std. Error")
      adrep.est[[i]]      <- as.list(em$sdrep, "Estimate", report = TRUE)
      adrep.se[[i]]       <- as.list(em$sdrep, "Std. Error", report = TRUE)
      opt_list[[i]]       <- em$opt
      converge_list[[i]]  <- cc$conv + cc$pdHess
      catch_advice[[i]]   <- advice
      catch_realized[[i]] <- real_catch
      
      if (isTRUE(save.sdrep)) {
        em_full[[i]] <- em
      } else {
        if (y == tail(assess_years_use, 1)) em_full[[1]] <- em
        if (!isTRUE(save.last.em)) em_full[[1]] <- list()
      }
      
      em_input_list[[i]] <- em_input
    }
    
    end.time   <- Sys.time()
    time.taken <- end.time - start.time
    
    cat("\n============================================================\n")
    cat("Runtime for proj_NAA_opt =", proj_opt_k, ":", time.taken, "\n")
    cat("============================================================\n")
    
    out_mse[[k]] <- list(
      proj_NAA_opt = proj_opt_k,
      settings = list(
        year_start = year_start,
        year_end   = year_end,
        n_feedback_years = n_feedback_years,
        Rec_sigma = Rec_sigma,
        NAA_sigma = NAA_sigma,
        use_ecov_rec_link = use_ecov_rec_link,
        EM_fit   = EM_fit,
        full_MSE = full_MSE
      ),
      om_final      = om,
      em_list       = em_list,
      par.est       = par.est,
      par.se        = par.se,
      adrep.est     = adrep.est,
      adrep.se      = adrep.se,
      opt_list      = opt_list,
      converge_list = converge_list,
      catch_advice  = catch_advice,
      catch_realized= catch_realized,
      em_full       = em_full,
      em_input      = em_input_list,
      runtime       = time.taken,
      seed.save     = seed
    )
  }
  
  # Optional: save all 3
  # saveRDS(out_mse, file = p("models", "mse_out_projopts_1_2_3.rds"))
}

## =============================================================================
## 16) Minimal comparison helpers (optional quick checks)
## =============================================================================

if (isTRUE(full_MSE)) {
  cat("\n============================================================\n")
  cat("Quick comparison: total catch advice by scenario (sum over fleets)\n")
  cat("============================================================\n")
  
  for (nm in names(out_mse)) {
    ca <- out_mse[[nm]]$catch_advice
    tot <- sapply(ca, function(x) if (is.null(x)) NA_real_ else sum(x))
    cat("\n", nm, " total advised catch (sum over all assessed intervals):\n", sep = "")
    print(tot)
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)

# --- Build a long dataframe for plotting
df <- bind_rows(
  data.frame(
    Year = out_mse[["proj_NAA_opt_1"]]$om_final$years,
    R1   = out_mse[["proj_NAA_opt_1"]]$om_final$rep$SSB[, 1],
    R2   = out_mse[["proj_NAA_opt_1"]]$om_final$rep$SSB[, 2],
    Case = "proj_NAA_opt_1"
  ),
  data.frame(
    Year = out_mse[["proj_NAA_opt_2"]]$om_final$years,
    R1   = out_mse[["proj_NAA_opt_2"]]$om_final$rep$SSB[, 1],
    R2   = out_mse[["proj_NAA_opt_2"]]$om_final$rep$SSB[, 2],
    Case = "proj_NAA_opt_2"
  ),
  data.frame(
    Year = out_mse[["proj_NAA_opt_3"]]$om_final$years,
    R1   = out_mse[["proj_NAA_opt_3"]]$om_final$rep$SSB[, 1],
    R2   = out_mse[["proj_NAA_opt_3"]]$om_final$rep$SSB[, 2],
    Case = "proj_NAA_opt_3"
  )
) |>
  pivot_longer(cols = c(R1, R2), names_to = "Region", values_to = "SSB")

# --- Plot: 2 facets (regions), lines for 3 cases
ggplot(df, aes(x = Year, y = SSB, col = Case)) +
  geom_line() +
  facet_wrap(~Region, ncol = 1, scales = "free_y") +
  labs(x = "Year", y = "SSB") +
  theme_bw()


# --- Build a long dataframe for plotting Catch
df <- bind_rows(
  data.frame(
    Year = out_mse[["proj_NAA_opt_1"]]$om_final$years,
    Fleet1   = out_mse[["proj_NAA_opt_1"]]$om_final$rep$pred_catch[, 1],
    Fleet2   = out_mse[["proj_NAA_opt_1"]]$om_final$rep$pred_catch[, 2],
    Fleet3   = out_mse[["proj_NAA_opt_1"]]$om_final$rep$pred_catch[, 3],
    Fleet4   = out_mse[["proj_NAA_opt_1"]]$om_final$rep$pred_catch[, 4],
    Case = "proj_NAA_opt_1"
  ),
  data.frame(
    Year = out_mse[["proj_NAA_opt_2"]]$om_final$years,
    Fleet1   = out_mse[["proj_NAA_opt_2"]]$om_final$rep$pred_catch[, 1],
    Fleet2   = out_mse[["proj_NAA_opt_2"]]$om_final$rep$pred_catch[, 2],
    Fleet3   = out_mse[["proj_NAA_opt_2"]]$om_final$rep$pred_catch[, 3],
    Fleet4   = out_mse[["proj_NAA_opt_2"]]$om_final$rep$pred_catch[, 4],
    Case = "proj_NAA_opt_2"
  ),
  data.frame(
    Year = out_mse[["proj_NAA_opt_3"]]$om_final$years,
    Fleet1   = out_mse[["proj_NAA_opt_3"]]$om_final$rep$pred_catch[, 1],
    Fleet2   = out_mse[["proj_NAA_opt_3"]]$om_final$rep$pred_catch[, 2],
    Fleet3   = out_mse[["proj_NAA_opt_3"]]$om_final$rep$pred_catch[, 3],
    Fleet4   = out_mse[["proj_NAA_opt_3"]]$om_final$rep$pred_catch[, 4],
    Case = "proj_NAA_opt_3"
  )
) |>
  pivot_longer(cols = c(Fleet1, Fleet2, Fleet3, Fleet4), names_to = "Fleet", values_to = "Catch")

ggplot(df, aes(x = Year, y = Catch, col = Case)) +
  geom_line() +
  facet_wrap(~Fleet, ncol = 2, scales = "free_y") +
  labs(x = "Year", y = "Catch") +
  theme_bw()



# - If everything is working, then you can use loop_through_fn to run MSE
proj.opts <- list(proj_NAA_opt = 1)  # example

mod <- loop_through_fn(
  om = om_with_data,
  em_info = info,
  random = random,
  
  # EM components
  sel_em = sel,
  M_em = M,
  NAA_re_em = NAA_re_em,
  move_em = move_em,
  ecov_em = ecov_em,
  
  # EM control
  em.opt = list(
    separate.em = FALSE,
    separate.em.type = NULL,
    do.move = TRUE,
    est.move = TRUE
  ),
  
  # Projection control
  proj.opts = proj.opts,
  
  # Update data CV/Neff used by EM
  update_catch_info = list(
    agg_catch_sigma = input_Ecov$data$agg_catch_sigma,
    catch_Neff      = input_Ecov$data$catch_Neff
  ),
  update_index_info = list(
    agg_index_sigma = input_Ecov$data$agg_index_sigma,
    index_Neff      = input_Ecov$data$index_Neff
  ),
  
  # Timeline
  assess_years    = assess.years,
  assess_interval = assess.interval,
  base_years      = base.years,
  year.use        = length(base.years),
  add.years       = TRUE,
  
  # MSE settings
  seed = seed_main,
  hcr  = hcr,
  save.last.em = TRUE
)

