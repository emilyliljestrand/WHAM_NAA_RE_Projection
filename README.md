# WHAM NAA Random-Effects Projection

Simulation-estimation and retrospective-projection analyses for evaluating how
assumptions about recruitment and numbers-at-age (NAA) random effects affect
Black Sea Bass (*Centropristis striata*) projections and spawning-stock biomass
(SSB) estimates.

## Repository contents

| Path | Contents |
| --- | --- |
| [`R/BSB_analysis/0.Test_w_BSB_Data.R`](R/BSB_analysis/0.Test_w_BSB_Data.R) | Fits the full 1989–2024 Black Sea Bass model, fits a truncated 1989–2021 model, projects 2022–2024 under three random-effects options, and calculates relative SSB bias. |
| [`R/BSB_analysis/1.Sim_Est_BSB.R`](R/BSB_analysis/1.Sim_Est_BSB.R) | Runs the 10-replicate simulation-estimation experiment: simulates observations, fits full and truncated models, projects the truncated fits, and saves comparisons. |
| [`R/functions/project_paths.R`](R/functions/project_paths.R) | Utilities for project-relative paths and output directories. |
| [`config/mse_settings.R`](config/mse_settings.R) | General SPASAM.MSE settings list; it is configuration support rather than the current analysis entry point. |
| [`data/raw/asap/`](data/raw/asap/) | North and South ASAP inputs for 1989–2021 and 1989–2024. |
| [`data/raw/covariates/`](data/raw/covariates/) | North and South bottom-temperature covariates for 1959–2021 and 1959–2024. |
| [`data/raw/seeds/`](data/raw/seeds/) | Simulation seeds for the simulation-estimation workflow. |
| [`models/`](models/) | Saved fitted WHAM models and the shared BSB model configuration. |
| [`output/BSB_analysis/`](output/BSB_analysis/) | Generated WHAM reports, relative-bias results, simulation data, convergence summaries, and replicate model files. |
| [`writing/manuscript-first-draft.qmd`](writing/manuscript-first-draft.qmd) | Quarto manuscript scaffold. |
| [`writing/references.bib`](writing/references.bib) | Manuscript bibliography. |

## Analysis overview

The analyses use two spatially structured Black Sea Bass stocks/regions,
North and South, with four fleets, four indices, 11 seasons, and eight ages.
The retrospective analysis withholds 2022–2024 from the estimation model and
compares projections with the full 1989–2024 fit. The three projection options
are:

1. Continue random effects on recruitment and NAA.
2. Turn off random effects on recruitment and NAA.
3. Average recruitment and NAA random effects over a historical reference period.

The simulation-estimation analysis treats `models/BSB.EM.Y.RDS` as the
operating-model truth, simulates catch and index observations, fits full and
1989–2021 estimation models, and projects the truncated model through 2024.

## Requirements

The scripts require R and the following packages:

- `wham`, using the pinned local build referenced in the scripts (`2.1.0.9003`)
- [`SPASAM.MSE`](https://github.com/lichengxue/SPASAM.MSE)
- `tidyverse`
- `here`

The simulation-estimation workflow also requires the saved operating model and
configuration files in `models/`. The fitting calls in the retrospective script
are commented out, so existing `.RDS` files are read from that directory.

## Running the analyses

Run from the repository root so relative paths resolve correctly:

```r
source("R/BSB_analysis/0.Test_w_BSB_Data.R")
source("R/BSB_analysis/1.Sim_Est_BSB.R")
```


## Main outputs

The retrospective analysis writes the relative-bias summary to
[`output/BSB_analysis/0.Test_w_BSB_data/ssb_relative_bias_2022_2024.csv`](output/BSB_analysis/0.Test_w_BSB_data/ssb_relative_bias_2022_2024.csv)
and generates WHAM diagnostic output under the same directory.

The simulation-estimation analysis writes model-convergence summaries,
simulated catch/index data, fitted model objects, and projection-comparison
outputs under [`output/BSB_analysis/1.Sim_Est_BSB/`](output/BSB_analysis/1.Sim_Est_BSB/).

Generated outputs and `.RDS` files are excluded from version control by
[`.gitignore`](.gitignore), although these directories may exist locally after
an analysis run.

## NOAA disclaimer

This repository is a scientific product and is not official communication of
the National Oceanic and Atmospheric Administration or the United States
Department of Commerce. All NOAA GitHub project code is provided on an “as is”
basis, and the user assumes responsibility for its use. References to specific
commercial products, processes, or services do not constitute or imply
endorsement, recommendation, or favoring by the Department of Commerce.
