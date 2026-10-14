# Quantitative Microbial Risk Assessment (QMRA) of Human Norovirus in Cooked Oysters

This repository provides the computational workflow for a two-dimensional Monte Carlo Quantitative Microbial Risk Assessment (2D-QMRA) evaluating human norovirus genogroup II (HuNoV GII) exposure and illness risks from cooked Eastern oysters (*Crassostrea virginica*).

The model couples dynamic thermal inactivation kinetics with empirical oyster core time–temperature profiles, a Bayesian contamination model, genome-copy-to-infectious-virus conversion, host susceptibility, and human dose–response modeling across domestic baking and steaming scenarios.

---

## Model Workflow

The QMRA framework is organized into four sequential modules:

1. **Module 1: Exposure Assessment**
   - **1.1 Initial Contamination:** Bayesian hurdle model for oyster-level GII prevalence and censored concentration.
   - **1.2 Thermal Inactivation:** Dynamic non-isothermal inactivation (Geeraerd and Mafart–Weibull models) integrated with measured internal heating curves across small, medium, and large oysters.
   - **1.3 Dose Calculation:** Conversion of genome copies to infectious virus ($Q$) and aggregation into meal-level exposures (6-oyster servings).
2. **Module 2: Hazard Characterization**
   - Host genetic susceptibility (secretor status, $P(\text{Se}+)$).
   - Teunis et al. (2020) dose–response model for infection and illness given infection.
3. **Module 3: Risk Characterization**
   - 2D Monte Carlo integration (10,000 variability $\times$ 2,000 uncertainty iterations).
   - Illness probability per serving and expected annual U.S. disease burden across cooking scenarios.
4. **Module 4: Sensitivity Analysis**
   - Global Spearman rank correlation (tornado sensitivity analysis).
   - Regional sensitivity analysis via Monte Carlo Filtering (MCF) with Benjamini–Hochberg FDR control.

---

## Repository Structure

```
├── README.md                                    # Project overview
├── oyster-norovirus-qmra.Rproj                  # RStudio project file
├── data/
│   ├── experimental_temperature_profiles/       # Core temperature curves for baking and steaming
│   ├── literature_norovirus_contamination/      # Prevalence and RT-qPCR concentration data
│   ├── module1_2_thermal_inactivation/          # Thermal inactivation experimental data
│   └── module4_foodnet_oyster_consumption_inputs.csv # FoodNet consumption parameters
├── scripts/
│   ├── oyster_qmra_master.qmd                   # Master Quarto computational notebook
│   ├── module1_exposure_assessment_thermal_inactivation_functions.R
│   ├── module2_hazard_characterization_functions.R
│   ├── module3_risk_characterization_functions.R
│   ├── module4_sensitivity_analysis_functions.R
│   └── refit_oyster_geeraerd_model.R            # Bayesian MCMC thermal kinetics fitting
├── results/                                     # Simulation draws, posterior objects, and CSV summaries
├── figures/                                     # Output publication figures
└── tables/                                      # Model parameter tables
```
