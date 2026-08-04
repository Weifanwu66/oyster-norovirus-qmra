# Oyster Norovirus GII QMRA

Working repository for a two-dimensional Monte Carlo assessment of norovirus
GII risk from raw and cooked oysters. The model is under development.

The current workflow includes:

1. initial prevalence and positive-oyster concentration;
2. genome-copy-to-infectious-unit conversion;
3. thermal inactivation;
4. 6- and 12-oyster serving simulations; and
5. dose-response and annual illness estimates.

The analysis is in
[`scripts/oyster_qmra_master.qmd`](scripts/oyster_qmra_master.qmd). Open
`oyster-norovirus-qmra.Rproj` in RStudio and render the QMD to run it.

The current Monte Carlo dimensions are 10,000 variability iterations and 500
uncertainty iterations. These will be checked through convergence testing.
