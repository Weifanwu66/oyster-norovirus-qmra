# Quantitative Microbial Risk Assessment (QMRA) of Human Norovirus GII in Cooked Oysters

This repository contains the complete computational workflow for a two-dimensional Monte Carlo Quantitative Microbial Risk Assessment (2D-QMRA) evaluating human norovirus genogroup II (HuNoV GII) exposure and illness risks from commercially harvested Eastern oysters (*Crassostrea virginica*) subjected to domestic and retail cooking methods (baking and steaming).

The framework couples dynamic thermal inactivation kinetics with empirical oyster core time-temperature heating profiles, a Bayesian hurdle contamination model, an uncertain genome-copy-to-infectious-unit conversion, host genetic susceptibility, and an epidemiological dose-response model.

---

## Model Architecture and Workflow

The quantitative model is structured into four sequential modules:

```
[Module 1.1: Initial Contamination]
  - Bayesian Hurdle Model (Prevalence + Censored Concentration)
  - Retail genome copy (GC) quantification
         │
         ▼
[Module 1.2 & 1.3: Exposure Assessment]
  - Tulane virus surrogate conversion (Q: GC to infectious units)
  - Dynamic non-isothermal thermal inactivation (Geeraerd & Mafart-Weibull)
  - Measured baking & steaming core profiles across small, medium, and large oysters
  - Serving aggregation (Standard 6-oyster servings: homogeneous & mixed-size classes)
         │
         ▼
[Module 2: Hazard Characterization]
  - Teunis et al. (2020) HuNoV GII dose-response model
  - Challenge trials (infection) & acute outbreak data (conditional illness)
  - Host secretor status susceptibility (Se+)
         │
         ▼
[Module 3: Risk Characterization]
  - 2D Monte Carlo integration (10,000 variability × 2,000 uncertainty iterations)
  - Per-serving illness probability across cooking conditions
  - Annual national disease burden:
      * Raw consumption frequency (124.4M events) -> Epidemiological validation benchmark (~58,000 cases)
      * Cooked cohort unheated counterfactual (337.7M events) -> Cooking risk reduction baseline (~157,000 cases)
         │
         ▼
[Module 4: Uncertainty & Sensitivity Analysis]
  - Global Spearman rank correlation Tornado plots (single medium oyster steamed for 4 min)
  - Monte Carlo Filtering (MCF) regional sensitivity using a 10⁻⁴ operational risk benchmark
```

---

## Repository Structure

```
├── README.md                                    # Project documentation
├── oyster-norovirus-qmra.Rproj                  # RStudio project configuration
├── data/
│   ├── experimental_temperature_profiles/       # Measured core temperature curves during baking & steaming
│   ├── literature_norovirus_contamination/      # Retail oyster prevalence and RT-qPCR concentration datasets
│   └── module4_foodnet_oyster_consumption_inputs.csv # CDC FoodNet 7-day recall consumption parameters
├── scripts/
│   ├── oyster_qmra_master.qmd                   # Master Quarto computational document
│   ├── module1_exposure_assessment_thermal_inactivation_functions.R # Contamination & thermal inactivation
│   ├── module2_hazard_characterization_functions.R                  # Dose-response & host susceptibility
│   ├── module3_risk_characterization_functions.R                    # 2D-MC risk integration & annual burden
│   ├── module4_sensitivity_analysis_functions.R                     # Tornado correlation & Monte Carlo filtering
│   └── refit_oyster_geeraerd_model.R            # Bayesian MCMC thermal inactivation kinetic fitting
├── results/                                     # Saved model objects, posterior distributions, and summary CSVs
│   ├── module1_1_initial_contamination/
│   ├── module1_2_thermal_inactivation/
│   ├── module1_3_exposure_dose_calculation/
│   ├── module3_risk_characterization/
│   └── module4_sensitivity_analysis/
├── figures/                                     # High-resolution output figures (TIFF, 300 DPI)
└── tests/
    └── test_rng_reproducibility.R               # RNG reproducibility unit tests
```

---

## Key Methodological Features

1. **Contamination Model (Module 1.1)**:
   - Pooled harvest prevalence: Beta(64, 1641) with a posterior mean of 3.75% [2.91%, 4.72%].
   - Contamination concentration is modeled via a Bayesian left-censored framework with an analytical limit of quantification (LOQ = 100 GC/oyster). A Gamma distribution (shape = 0.152, scale = 3381 GC) was selected as the baseline model to capture right-skewed contamination while avoiding unconstrained heavy-tail inflation in unobserved percentiles.

2. **Infectious Unit Conversion (Q)**:
   - Directly evaluating retail oyster RT-qPCR copies against clinical dose-response curves overestimates raw-oyster risk by ~35-fold because retail qPCR captures degraded, non-viable environmental RNA fragments.
   - Following regulatory QMRA precedents (FDA–Health Canada), genome copies are adjusted to infectious units prior to thermal inactivation using empirical surrogate data: log10 Q ~ Uniform(3.39, 4.22) (Mirmahdi et al., 2025).

3. **Dynamic Thermal Inactivation (Module 1.2)**:
   - Evaluates real-time inactivation kinetics integrated over non-isothermal internal heating trajectories:
     - **Primary Baseline (Geeraerd Model)**: Features a shoulder effect and an asymptotic resistant tail (N_min), bounding maximum inactivation to ~2.56-log reduction (~0.28% residual survival).
     - **Sensitivity Model (Mafart-Weibull Model)**: Features continuous non-linear tailing across holding times.
   - Empirical heating data cover baking (2–10 min) and steaming (2–16 min) across three market size grades: small (<60 g), medium (60–90 g), and large (>90 g).

4. **Serving Size & Dual Population Baseline (Module 3)**:
   - Standard serving size is defined as **6 oysters** (half-dozen), with 12 oysters evaluated to demonstrate linear scaling at low exposure doses.
   - **Epidemiological Validation Benchmark**: Applying the raw-oyster eating frequency (0.007/week, 124.4 million events/year) yields 57,934 [49, 405,436] annual cases, aligning with CDC surveillance estimates and FDA risk assessment benchmarks.
   - **Cooking Risk Reduction Baseline**: Evaluates the cooked-oyster consumer cohort (0.019/week, 337.7 million events/year) under an unheated counterfactual scenario (~157,000 cases). Thermal cooking reduces this burden to 53 annual cases under the Geeraerd plateau, or to 1,079–4,199 cases under Mafart-Weibull.

5. **Sensitivity Analysis Endpoints (Module 4)**:
   - Evaluated at the **single-oyster dose level** (a single medium oyster steamed for 4 min) to isolate parameter uncertainty from multinomial serving variability.
   - **Tornado Plots**: Non-parametric Spearman rank correlation coefficients.
   - **Monte Carlo Filtering (MCF)**: Regional sensitivity partitioning outputs using an operational risk benchmark of 10⁻⁴ per eating occasion, followed by two-sample Kolmogorov-Smirnov tests with Benjamini-Hochberg FDR adjustments.

---

## Reproducibility & RNG Design

- **Random Number Generation**: Explicit seeds are set alongside each stochastic simulation block using `withr::local_seed()` with the `Mersenne-Twister / Inversion / Rejection` algorithm.
- **Common Random Numbers (CRN)**: Shared sampling indices are preserved across baseline and thermal scenarios to ensure that scenario comparisons remain independent of simulation sequence.
- **Monte Carlo Dimensions**: Production runs evaluate **10,000 variability iterations** and **2,000 uncertainty iterations** (20,000,000 iterations per scenario).

To run the RNG reproducibility verification suite:

```sh
Rscript tests/test_rng_reproducibility.R
```

---

## Getting Started

### Prerequisites

- R (version >= 4.3.0)
- Quarto CLI (version >= 1.4.0)
- Required R packages:
  ```r
  install.packages(c(
    "tidyverse", "here", "readr", "dplyr", "ggplot2", "scales",
    "patchwork", "withr", "mc2d", "brms", "rstan", "loo", "extraDistr"
  ))
  ```

### Running the Analysis

1. Clone the repository:
   ```sh
   git clone https://github.com/Weifanwu66/oyster-norovirus-qmra.git
   cd oyster-norovirus-qmra
   ```
2. Open `oyster-norovirus-qmra.Rproj` in RStudio.
3. Open [`scripts/oyster_qmra_master.qmd`](scripts/oyster_qmra_master.qmd) and render the document:
   - In RStudio: Click the **Render** button.
   - Via terminal:
     ```sh
     quarto render scripts/oyster_qmra_master.qmd
     ```

*Note: Large MCMC model refits and 2D-MC simulations have pre-computed outputs saved in `results/`. Set `eval=TRUE` on specific compute chunks in the QMD to execute a full regeneration from scratch.*

---

## Citation

If you use this model or code in your research, please cite:

> Wu, W., et al. (2026). Quantitative microbial risk assessment of human norovirus genogroup II in cooked Eastern oysters (*Crassostrea virginica*): Evaluating dynamic non-isothermal thermal inactivation and model uncertainty. *Risk Analysis* (In review).
