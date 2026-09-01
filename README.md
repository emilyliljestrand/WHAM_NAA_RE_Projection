# WHAM_NAA_RE_Projection :fish:
Simulation-Estimation to test assumed NAA RE in projections on management outcomes

<img src="https://www.fisheries.noaa.gov/s3/2022-08/640x427-Black-Sea-Bass-NOAAFisheries.png" height="450" alt="Black Sea Bass">

# **File Descriptions**

## Files

| Path | Purpose |
| --- | --- |
| `R/01_run_mse.R` | Runs a configurable two-stock, two-region SPASAM.MSE example and saves each strategy/replicate result. |
| `R/02_plot_mse.R` | Reads the saved MSE results and creates an HTML performance report. |
| `output/mse_results/` | Created by the runner; contains `run_metadata.rds` and per-strategy replicate outputs. |
| `output/mse_report/` | Created by the plotting script; contains the rendered report and figures. |

## MSE workflow

The MSE scripts use [SPASAM.MSE](https://github.com/lichengxue/SPASAM.MSE). Install R package dependencies once:

```r
install.packages("remotes")
remotes::install_github("lichengxue/SPASAM.MSE", dependencies = TRUE)
```

From the repository root, run the scripts in order:

```r
source("R/01_run_mse.R")
source("R/02_plot_mse.R")
```

Edit the experiment settings at the top of `R/01_run_mse.R` to choose seeds,
replicate count, and $F_{X\%SPR}$ strategies. The initial configuration is a
small generic two-stock example so the package workflow can be tested before
moving the Black Sea Bass model inputs into the new structure.

## NOAA Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and
Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is
provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of
Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial products, processes, or services by service
mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or
favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a
DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by
DOC or the United States Government.

****************************

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) | [National Oceanographic and Atmospheric Administration](https://www.noaa.gov) | [NOAA Fisheries](https://www.fisheries.noaa.gov/)


