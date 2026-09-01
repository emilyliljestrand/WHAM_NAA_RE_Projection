#' @title Retrospective Forecasting Proof of Concept with Black Sea Bass Data
#' @description Peels back 3 years (2022-2024) from the 2024 Black Sea Bass (BSB) 
#'   stock assessment model (1989-2024) to generate retrospective forecasts under 
#'   three projection options and evaluate bias against full model estimates.
#' @details Evaluates three projection random effect options:
#'   1. Continue random effects on both recruitment (R) and Numbers-at-Age (NAA).
#'   2. Turn off random effects on both R and NAA.
#'   3. Average over specified historical recruitment and NAA years.
#'   Compares relative bias in Spawning Stock Biomass (SSB) estimates across North and South regions.
#' @author Emily Liljestrand
#' @name test_w_bsb_data
NULL

#' # ==============================================================================
#' # 1) Environment Setup & Libraries
#' # ==============================================================================

# Clean workspace environment
rm(list=ls())

# Load required packages for WHAM, TMB optimization, and data visualization
library(wham, lib.loc = "C:/Users/emily.liljestrand/AppData/Local/R/win-library/4.4/wham_EML")
library(tidyverse)

#' # ==============================================================================
#' # 2) Full Estimation / Operating Model (1989-2024)
#' # ==============================================================================
asap <- read_asap3_dat(file.path("data", "raw", "asap", c("NORTH.1989.2024.DAT","SOUTH.1989.2024.DAT")))
north_bt <- read.csv(file.path("data", "raw", "covariates", "bsb_bt_temp_nmab_1959-2024.csv"))
south_bt <- read.csv(file.path("data", "raw", "covariates", "bsb_bt_temp_smab_1959-2024.csv"))
NAA_re = list(sigma = list("rec+1","rec+1"), cor = list("2dar1","2dar1"), N1_model = rep("equilibrium",2))
NAA_re$decouple_recruitment = TRUE
NAA_re$sigma_vals <- array(1,dim = c(2,2,8))
NAA_re$sigma_vals[1,2,2:8] <- 0.05
NAA_re$sigma_map <- array(1,dim = c(2,2,8))
NAA_re$sigma_map[1,1,2:8] <- 2
NAA_re$sigma_map[2,,1] <- 3
NAA_re$sigma_map[2,,2:8] <- 4
NAA_re$sigma_map[1,2,2:8] <- NA
NAA_re$cor_vals <- array(0,dim = c(2,2,3))
NAA_re$cor_map <- array(NA,dim = c(2,2,3))
NAA_re$cor_map[1,1,1] <- 1 #stock 1, region 1, rho_a
NAA_re$cor_map[1,1,2] <- 2 #stock 1, region 1, rho_y 1
NAA_re$cor_map[1,1,3] <- 3 #stock 1, region 1, rho_y 2-8
NAA_re$cor_map[2,,1] <- 4 #stock 2, rho_a
NAA_re$cor_map[2,,2] <- 5 #stock 2, rho_y
NAA_re$cor_map[2,,3] <- 6 #stock 2, rho_y 2-8
ecov <- list(label = c("North_BT","South_BT"))
ecov$mean <- cbind(north_bt[,'mean'], south_bt[,'mean'])
ecov$logsigma <- log(cbind(north_bt[,'se'], south_bt[,'se']))
ecov$year <- north_bt[,'year']
ecov$use_obs <- matrix(1, NROW(ecov$mean),NCOL(ecov$mean))
ecov$process_model <- "ar1"
ecov$process_mean_vals <- apply(ecov$mean, 2, mean)
ecov$recruitment_how <- matrix(c("controlling-lag-0-linear","none","none","none"), 2,2)
temp <- prepare_wham_input(asap, ecov = ecov, NAA_re = NAA_re)
seasons = c(rep(1,5),2,rep(1,5))/12
basic_info <- list(region_names = c("North", "South"), stock_names = paste0("BSB_", c("North", "South"))) #, NAA_where = array(1, dim = c(2,2,6)))
basic_info$fracyr_seasons <- seasons
basic_info$NAA_where <- array(1, dim = c(2,2,8))
basic_info$NAA_where[1,2,1] = 0 #stock 1, age 1 can't be in region 2 
basic_info$NAA_where[2,1,] = 0 #stock 2, any age can't be in region 1 (stock 2 doesn't move) 
basic_info$XSPR_R_avg_yrs <- which(temp$years>1999)
basic_info$XSPR_R_opt <- 2 #use average of recruitments (random effects), not expected/predicted given last time step
move = list(stock_move = c(TRUE,FALSE), separable = TRUE) #north moves, south doesn't
move$must_move = array(0,dim = c(2,length(seasons),2))  
move$must_move[1,5,2] <- 1 
move$can_move = array(0, dim = c(2,length(seasons),2,2))
move$can_move[1,c(1:4),2,1] <- 1 #only north stock can move and in seasons prior to spawning and after spawning
move$can_move[1,c(7:11),1,2] <- 1 #only north stock can move and in seasons prior to spawning and after spawning
move$can_move[1,5,2,] <- 1 #north stock can (and must) move in last season prior to spawning back to north 
mus <- array(0, dim = c(2,length(seasons),2,1))
mus[1,1:11,1,1] <- 0.02214863 #see here("2023.RT.Runs","transform_SS_move_rates.R") for how these numbers are derived. Priors
mus[1,1:11,2,1] <- 0.3130358
move$mean_vals <- mus 
move$mean_model = matrix("stock_constant", 2,1)
move$use_prior <- array(0, dim = c(2,length(seasons),2,1))
move$use_prior[1,1,1,1] <- 1
move$use_prior[1,1,2,1] <- 1
move$prior_sigma <- array(0, dim = c(2,length(seasons),2,1))
move$prior_sigma[1,1,1,1] <- 0.2
move$prior_sigma[1,1,2,1] <- 0.2
catch_info <- list(
  catch_Neff = matrix(1000, length(temp$years), temp$data$n_fleets), 
  selblock_pointer_fleets =  matrix(rep(1:4, each = length(temp$years)), length(temp$years)),
  fleet_names = paste0(rep(c("North_", "South_"),each = 2), temp$fleet_names))
index_info <- list(index_Neff = matrix(1000, length(temp$years), temp$data$n_indices),
                   selblock_pointer_indices =  matrix(rep(5:8, each = length(temp$years)), length(temp$years)),
                   index_names = paste0(rep(c("North_", "South_"),each = 2), temp$index_names),
                   initial_index_sd_scale = c(5,1,5,1),
                   map_index_sd_scale = c(1,NA,2,NA))
age_comp = list(
  fleets = c("dir-mult","logistic-normal-miss0","logistic-normal-ar1-miss0","logistic-normal-ar1-miss0"), 
  indices = c("logistic-normal-miss0","dir-mult","logistic-normal-ar1-miss0","logistic-normal-ar1-miss0"))
sel <- list(n_selblocks = 8, model = rep(c("age-specific","logistic","age-specific","age-specific"),
                                         c(2,2,3,1)))
sel$initial_pars <- list(
  rep(c(0.5,1),c(3,5)), #north comm
  rep(c(0.5,1),c(6,2)), #north rec
  c(5,1), #south comm
  c(5,1), #south rec
  rep(c(0.5,1,1),c(1,1,6)), #north rec cpa
  rep(c(0.5,1),c(4,4)), #north vast
  rep(c(0.5,1,1),c(2,4,2)), #south rec cpa
  rep(c(0.5,1),c(1,7)) #south vast
)
sel$fix_pars <- list(
  4:8, #north comm
  7:8, #north rec
  NULL, #south comm
  NULL, #south rec
  2:8, #north rec cpa
  5:8, #north vast
  3:8, #south rec cpa
  2:8 #south vast
)
sel$re <- rep(c("2dar1","none","ar1_y","2dar1","none"), c(2,2,1,1,2))
temp <- prepare_wham_input(asap, ecov = ecov, NAA_re = NAA_re, basic_info = basic_info, move = move, catch_info = catch_info, index_info = index_info, age_comp = age_comp, selectivity = sel)

# BSB.EM.Y <- fit_wham(temp, do.sdrep = T, do.osa = T, do.retro = T, do.brps = T)
# saveRDS(BSB.EM.Y, "BSB.EM.Y.RDS")
BSB.EM.Y <- readRDS("models/BSB.EM.Y.RDS")
plot_wham_output(BSB.EM.Y)

#' # ==============================================================================
#' # 3) Reduced Estimation / Retrospective Model (1989-2021)
#' # ==============================================================================

# Load truncated ASAP datasets peeling off 3 terminal years (2022-2024)
asap <- read_asap3_dat(file.path("data", "raw", "asap", c("NORTH.1989.2021.DAT","SOUTH.1989.2021.DAT")))
north_bt <- read.csv(file.path("data", "raw", "covariates", "bsb_bt_temp_nmab_1959-2021.csv"))
south_bt <- read.csv(file.path("data", "raw", "covariates", "bsb_bt_temp_smab_1959-2021.csv"))
NAA_re = list(sigma = list("rec+1","rec+1"), cor = list("2dar1","2dar1"), N1_model = rep("equilibrium",2))
NAA_re$decouple_recruitment = TRUE
NAA_re$sigma_vals <- array(1,dim = c(2,2,8))
NAA_re$sigma_vals[1,2,2:8] <- 0.05
NAA_re$sigma_map <- array(1,dim = c(2,2,8))
NAA_re$sigma_map[1,1,2:8] <- 2
NAA_re$sigma_map[2,,1] <- 3
NAA_re$sigma_map[2,,2:8] <- 4
NAA_re$sigma_map[1,2,2:8] <- NA
NAA_re$cor_vals <- array(0,dim = c(2,2,3))
NAA_re$cor_map <- array(NA,dim = c(2,2,3))
NAA_re$cor_map[1,1,1] <- 1 #stock 1, region 1, rho_a
NAA_re$cor_map[1,1,2] <- 2 #stock 1, region 1, rho_y 1
NAA_re$cor_map[1,1,3] <- 3 #stock 1, region 1, rho_y 2-8
NAA_re$cor_map[2,,1] <- 4 #stock 2, rho_a
NAA_re$cor_map[2,,2] <- 5 #stock 2, rho_y
NAA_re$cor_map[2,,3] <- 6 #stock 2, rho_y 2-8
ecov <- list(label = c("North_BT","South_BT"))
ecov$mean <- cbind(north_bt[,'mean'], south_bt[,'mean'])
ecov$logsigma <- log(cbind(north_bt[,'se'], south_bt[,'se']))
ecov$year <- north_bt[,'year']
ecov$use_obs <- matrix(1, NROW(ecov$mean),NCOL(ecov$mean))
ecov$process_model <- "ar1"
ecov$process_mean_vals <- apply(ecov$mean, 2, mean)
ecov$recruitment_how <- matrix(c("controlling-lag-0-linear","none","none","none"), 2,2)
temp <- prepare_wham_input(asap, ecov = ecov, NAA_re = NAA_re)
seasons = c(rep(1,5),2,rep(1,5))/12
basic_info <- list(region_names = c("North", "South"), stock_names = paste0("BSB_", c("North", "South"))) #, NAA_where = array(1, dim = c(2,2,6)))
basic_info$fracyr_seasons <- seasons
basic_info$NAA_where <- array(1, dim = c(2,2,8))
basic_info$NAA_where[1,2,1] = 0 #stock 1, age 1 can't be in region 2 
basic_info$NAA_where[2,1,] = 0 #stock 2, any age can't be in region 1 (stock 2 doesn't move) 
basic_info$XSPR_R_avg_yrs <- which(temp$years>1999)
basic_info$XSPR_R_opt <- 2 #use average of recruitments (random effects), not expected/predicted given last time step
move = list(stock_move = c(TRUE,FALSE), separable = TRUE) #north moves, south doesn't
move$must_move = array(0,dim = c(2,length(seasons),2))  
move$must_move[1,5,2] <- 1 
move$can_move = array(0, dim = c(2,length(seasons),2,2))
move$can_move[1,c(1:4),2,1] <- 1 #only north stock can move and in seasons prior to spawning and after spawning
move$can_move[1,c(7:11),1,2] <- 1 #only north stock can move and in seasons prior to spawning and after spawning
move$can_move[1,5,2,] <- 1 #north stock can (and must) move in last season prior to spawning back to north 
mus <- array(0, dim = c(2,length(seasons),2,1))
mus[1,1:11,1,1] <- 0.02214863 #see here("2023.RT.Runs","transform_SS_move_rates.R") for how these numbers are derived. Priors
mus[1,1:11,2,1] <- 0.3130358
move$mean_vals <- mus 
move$mean_model = matrix("stock_constant", 2,1)
move$use_prior <- array(0, dim = c(2,length(seasons),2,1))
move$use_prior[1,1,1,1] <- 1
move$use_prior[1,1,2,1] <- 1
move$prior_sigma <- array(0, dim = c(2,length(seasons),2,1))
move$prior_sigma[1,1,1,1] <- 0.2
move$prior_sigma[1,1,2,1] <- 0.2
catch_info <- list(
  catch_Neff = matrix(1000, length(temp$years), temp$data$n_fleets), 
  selblock_pointer_fleets =  matrix(rep(1:4, each = length(temp$years)), length(temp$years)),
  fleet_names = paste0(rep(c("North_", "South_"),each = 2), temp$fleet_names))
index_info <- list(index_Neff = matrix(1000, length(temp$years), temp$data$n_indices),
                   selblock_pointer_indices =  matrix(rep(5:8, each = length(temp$years)), length(temp$years)),
                   index_names = paste0(rep(c("North_", "South_"),each = 2), temp$index_names),
                   initial_index_sd_scale = c(5,1,5,1),
                   map_index_sd_scale = c(1,NA,2,NA))
age_comp = list(
  fleets = c("dir-mult","logistic-normal-miss0","logistic-normal-ar1-miss0","logistic-normal-ar1-miss0"), 
  indices = c("logistic-normal-miss0","dir-mult","logistic-normal-ar1-miss0","logistic-normal-ar1-miss0"))
sel <- list(n_selblocks = 8, model = rep(c("age-specific","logistic","age-specific","age-specific"),
                                         c(2,2,3,1)))
sel$initial_pars <- list(
  rep(c(0.5,1),c(3,5)), #north comm
  rep(c(0.5,1),c(6,2)), #north rec
  c(5,1), #south comm
  c(5,1), #south rec
  rep(c(0.5,1,1),c(1,1,6)), #north rec cpa
  rep(c(0.5,1),c(4,4)), #north vast
  rep(c(0.5,1,1),c(2,4,2)), #south rec cpa
  rep(c(0.5,1),c(1,7)) #south vast
)
sel$fix_pars <- list(
  4:8, #north comm
  7:8, #north rec
  NULL, #south comm
  NULL, #south rec
  2:8, #north rec cpa
  5:8, #north vast
  3:8, #south rec cpa
  2:8 #south vast
)
sel$re <- rep(c("2dar1","none","ar1_y","2dar1","none"), c(2,2,1,1,2))
temp <- prepare_wham_input(asap, ecov = ecov, NAA_re = NAA_re, basic_info = basic_info, move = move, catch_info = catch_info, index_info = index_info, age_comp = age_comp, selectivity = sel)

# BSB.EM.Y3 <- fit_wham(temp, do.sdrep = T, do.osa = T, do.retro = T, do.brps = T)
# saveRDS(BSB.EM.Y3, "BSB.EM.Y3.RDS")
BSB.EM.Y3 <- readRDS("models/BSB.EM.Y3.RDS")

#' # ==============================================================================
#' # 4) Three Retrospective Projection Scenarios (2022-2024)
#' # ==============================================================================

# For all options, project at F40% (proj_F_opt = 3) for 3 years (2022-2024)

# Option 1: Continue random effects on both recruitment (R) and Numbers-at-Age (NAA)
# BSB.EM.Y3.Proj.1 <- project_wham(BSB.EM.Y3,proj.opts=list(proj_R_opt=1,proj_NAA_opt=1,proj_F_opt=c(3,3,3)),check.version = F)

# Option 2: Turn off random effects on both recruitment (R) and Numbers-at-Age (NAA)
# BSB.EM.Y3.Proj.2 <- project_wham(BSB.EM.Y3,proj.opts=list(proj_R_opt=4,proj_NAA_opt=3,proj_F_opt=c(3,3,3)),check.version = F)

# Option 3: Average recruitment and NAA random effects over historical reference period
# BSB.EM.Y3.Proj.3 <- project_wham(BSB.EM.Y3,proj.opts=list(proj_R_opt=3,proj_NAA_opt=2,proj_F_opt=c(3,3,3)),check.version = F)

# Save projection objects
# saveRDS(BSB.EM.Y3.Proj.1, "BSB.EM.Y3.Proj.1.RDS")
# saveRDS(BSB.EM.Y3.Proj.2, "BSB.EM.Y3.Proj.2.RDS")
# saveRDS(BSB.EM.Y3.Proj.3, "BSB.EM.Y3.Proj.3.RDS")

BSB.EM.Y3.Proj.1 <- readRDS("models/BSB.EM.Y3.Proj.1.RDS")
BSB.EM.Y3.Proj.2 <- readRDS("models/BSB.EM.Y3.Proj.2.RDS")
BSB.EM.Y3.Proj.3 <- readRDS("models/BSB.EM.Y3.Proj.3.RDS")

#' # ==============================================================================
#' # 5) Forecast Performance Evaluation & Relative Bias Metrics
#' # ==============================================================================

# Compare projection models against the full 2024 estimation model estimates
mods <- list(FullModel = BSB.EM.Y,Option1 = BSB.EM.Y3.Proj.1,Option2 = BSB.EM.Y3.Proj.2,Option3 = BSB.EM.Y3.Proj.3)
compare_wham_models(mods,calc.aic = FALSE, do.table=F,plot.opts=list(which=c(1,6,7,8,9,10),kobe.yr=2021))

# Relative bias in terminal 3-year (2022-2024) SSB estimates in the North Region
mean((BSB.EM.Y3.Proj.1$rep$SSB[34:36,1] - BSB.EM.Y$rep$SSB[34:36,1])/BSB.EM.Y$rep$SSB[34:36,1]) # Option 1 (Full RE)
mean((BSB.EM.Y3.Proj.2$rep$SSB[34:36,1] - BSB.EM.Y$rep$SSB[34:36,1])/BSB.EM.Y$rep$SSB[34:36,1]) # Option 2 (No RE)
mean((BSB.EM.Y3.Proj.3$rep$SSB[34:36,1] - BSB.EM.Y$rep$SSB[34:36,1])/BSB.EM.Y$rep$SSB[34:36,1]) # Option 3 (Avg RE)

# Relative bias in terminal 3-year (2022-2024) SSB estimates in the South Region
mean((BSB.EM.Y3.Proj.1$rep$SSB[34:36,2] - BSB.EM.Y$rep$SSB[34:36,2])/BSB.EM.Y$rep$SSB[34:36,2]) # Option 1 (Full RE)
mean((BSB.EM.Y3.Proj.2$rep$SSB[34:36,2] - BSB.EM.Y$rep$SSB[34:36,2])/BSB.EM.Y$rep$SSB[34:36,2]) # Option 2 (No RE)
mean((BSB.EM.Y3.Proj.3$rep$SSB[34:36,2] - BSB.EM.Y$rep$SSB[34:36,2])/BSB.EM.Y$rep$SSB[34:36,2]) # Option 3 (Avg RE)
