# Oyster Norovirus GII QMRA

This repository contains the working quantitative microbial risk assessment for
norovirus GII in raw and cooked oysters. The analysis separates uncertainty
from serving-level variability with two-dimensional Monte Carlo simulation.

## Analysis workflow

The master analysis is [`scripts/oyster_qmra_master.qmd`](scripts/oyster_qmra_master.qmd).
It currently connects:

1. Bayesian prevalence and positive-concentration models for initial oyster contamination;
2. genome-copy-to-infectious-unit conversion;
3. posterior uncertainty in thermal survival;
4. 6- and 12-oyster serving simulations; and
5. dose-response and annual-risk calculations under development.

The default simulation dimensions are 10,000 variability iterations and 500
uncertainty iterations. Final dimensions will be selected through Monte Carlo
convergence testing.

## Requirements

- R and RStudio
- Quarto
- CmdStan
- R packages: `here`, `tidyverse`, `mc2d`, `brms`, `cmdstanr`, `posterior`,
  `loo`, `FME`, and `bioinactivation`

Install the CRAN packages with:

```r
install.packages(c(
  "here", "tidyverse", "mc2d", "brms", "posterior", "loo", "FME",
  "bioinactivation"
))
```

Install `cmdstanr` and CmdStan if they are not already available:

```r
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)
cmdstanr::install_cmdstan()
```

## Run the analysis

1. Clone the repository and open `oyster-norovirus-qmra.Rproj` in RStudio.
2. Confirm that CmdStan is available with `cmdstanr::cmdstan_version()`.
3. Render `scripts/oyster_qmra_master.qmd` from RStudio, or run:

```sh
quarto render scripts/oyster_qmra_master.qmd
```

The Module 1 Bayesian models are fitted during rendering. A compact set of
Module 3 posterior draws is included so the master analysis does not require a
separate thermal-model refit. To reproduce those draws, run:

```sh
Rscript scripts/refit_oyster_geeraerd_model.R
```

Generated results and figures are written to `results/` and `outputs/` and are
not tracked, apart from the Module 3 posterior input required by the master
analysis.

## Collaboration

Use a separate branch for each set of revisions. Pull the current default
branch before editing, commit the changes on the revision branch, push it, and
open a pull request for review before merging.
