# Table 1. Model parameters, mathematical formulations, and data sources for quantitative microbial risk assessment of human norovirus GII in raw and cooked oysters

*Arranged in the sequential workflow order of Figure 1 (Exposure Assessment [Stages 1–4] → Hazard Characterization [Stage 5] → Risk Characterization [Stage 6])*

| Parameter & Description | Symbol | Distribution / Mathematical Formulation | Mean / Point Estimate | 95% Uncertainty Interval | Dimension | Reference |
| :--- | :---: | :--- | :---: | :---: | :---: | :--- |
| 1. Harvest & Initial Contamination (Fig. 1: State 1)<br>*Model: $P_{\text{prev}} \sim \text{Beta}(\alpha, \beta)$; $C_+ \sim \text{Gamma}(k, \theta)$, $f(C_+) = \frac{1}{\Gamma(k)\theta^k} C_+^{k-1} e^{-C_+/\theta}$* | | | | | | |
| HuNoV GII harvest prevalence | $P_{\text{prev}}$ | $\text{Beta}(64, 1641)$ conjugate posterior | 0.0375 (3.75%) | [0.0291, 0.0470] | U | Nishida et al. (2003, 2007) |
| Concentration in contaminated oysters (GC/oyster) | $C_+$ | $\text{Gamma}(k, \theta)$ conditional on presence | — | — | U, V | Nishida et al. (2003, 2007) |
| Gamma shape parameter (dimensionless) | $k$ | Joint Bayesian posterior marginal | 0.152 | [0.090, 0.239] | U | Nishida et al. (2003, 2007) |
| Gamma scale parameter (GC/oyster) | $\theta$ | Joint Bayesian posterior marginal | 2,875 | [1,401, 7,488] | U | Nishida et al. (2003, 2007) |
| 2. Genome Copy to Infectious Virion Conversion (Fig. 1: State 2)<br>*Model: $\log_{10}(Q) \sim \text{Uniform}(3.39, 4.22)$; $Q = 10^{\log_{10}(Q)}$; $D_{\text{precook}} = C_{\text{oyster}} / Q$* | | | | | | |
| log10 virion conversion ratio (log10[GC / virion]) | $\log_{10}(Q)$ | $\text{Uniform}(3.39, 4.22)$ | 3.805 | [3.39, 4.22] | U | Mirmahdi et al. (2025) |
| Direct virion conversion ratio (GC / virion) | $Q$ | $10^{\log_{10}(Q)}$ antilog transformation | 6,383 | [2,455, 16,596] | U | Mirmahdi et al. (2025) |
| 3. Cooking Preparation & Dynamic Thermal Inactivation (Fig. 1: State 3)<br>*Models: Geeraerd: $\frac{dN}{dt} = -k_{\text{max}}(T(t)) \cdot [\frac{1}{1 + 10^{-\alpha_0}}] \cdot (N - N_{\text{min}})$; Mafart: $\log_{10}(S(t)) = -[\int \frac{dt}{\delta(T(t))}]^p$* | | | | | | |
| Cooking preparation regimes | — | Baking (230°C / 450°F) and steaming (100°C / 212°F) | Discrete scenarios | — | Scenario | This study |
| Oyster commercial size classes (total weight, g) | — | Small (35.6 ± 4.2 g), Medium (67.4 ± 6.1 g), Large (124.9 ± 11.3 g) | Discrete scenarios | — | Scenario | This study |
| Cooking heating duration (min) | $t$ | Steaming: 0–16 min (2-min steps); Baking: 0–10 min (2-min steps) | Discrete scenarios | — | Scenario | This study |
| Internal thermocouple time–temperature histories (°C) | $T(t)$ | Dynamic core temperature series for 3 sizes × 2 methods | Empirical profiles | — | Empirical | Nalanda et al. (2026) |
| Baseline model: Dynamic Geeraerd with tailing | — | $k_{\text{max}}(T(t)) = \frac{\ln 10}{D_{57^\circ\text{C}} \cdot 10^{-(T(t)-57)/z}}$ | Structural model | — | Model | Dicker et al. (2026) |
| Reference D-value at 57°C (min) | $D_{57^\circ\text{C}}$ | Log-linear inactivation rate at reference temperature | 8.43 | [5.46, 9.94] | U | Dicker et al. (2026) |
| Thermal sensitivity coefficient (°C) | $z$ | Temperature increase required for 1-log reduction in D | 6.71 | [3.98, 9.32] | U | Dicker et al. (2026) |
| Initial physiological state parameter | $\log_{10}(N_0)$ | Initial virus population state (log10 copies/g) | 3.46 | [3.19, 3.70] | U | Dicker et al. (2026) |
| Residual tailing density level | $\log_{10}(N_{\text{min}})$ | Asymptotic resistant tail density (log10 copies/g) | 0.90 | [0.37, 1.36] | U | Dicker et al. (2026) |
| Sensitivity model: Non-isothermal Mafart-Weibull | — | $\delta(T(t)) = D_{60^\circ\text{C}} \cdot 10^{-(T(t)-60)/z}$ | Structural model | — | Model | Wu et al. (2026) |
| Reference scale decimal reduction time at 60°C (min) | $D_{60^\circ\text{C}}$ | Weibull scale parameter at reference temperature | 4.30 | [1.86, 10.79] | U | Wu et al. (2026) |
| Thermal sensitivity coefficient (°C) | $z$ | Temperature change for 10-fold change in $\delta$ | 30.10 | [23.14, 42.04] | U | Wu et al. (2026) |
| Curvature shape parameter (dimensionless) | $p$ | Weibull power exponent ($p < 1$: upward concavity / tailing) | 0.387 | [0.370, 0.404] | U | Wu et al. (2026) |
| 4. Serving Size & Oyster Aggregation (Fig. 1: State 4)<br>*Model: $D_{\text{serving}} = \sum D_{\text{postcook},j}$; $(n_S, n_M, n_L) \sim \text{Multinomial}(n; 1/3, 1/3, 1/3), n \in \{6, 12\}$* | | | | | | |
| Serving size (number of oysters / serving) | $n$ | Fixed consumption clusters of 6 or 12 oysters | Discrete scenarios | — | Scenario | This study |
| Mixed-size serving composition (oysters / class) | $n_S, n_M, n_L$ | $\text{Multinomial}(n; 1/3, 1/3, 1/3), n \in \{6, 12\}$ | Equal ratio (1:1:1) | — | V | This study |
| 5. Hazard Characterization & Host Susceptibility (Fig. 1: State 5)<br>*Model: $P_{\text{inf}}(D|\text{Se}^+) = 1 - {}_1F_1(\alpha, \alpha+\beta, -D)$; $P_{\text{ill}|\text{inf}}(D) = 1 - (1 + \eta D)^{-r}$; $P_{\text{ill}} = P(\text{Se}^+) \cdot P_{\text{inf}} \cdot P_{\text{ill}|\text{inf}}$* | | | | | | |
| Secretor-positive host phenotype frequency | $P(\text{Se}^+)$ | $\text{Beta}(76.0, 29.7)$ [FUT2-functional proportion] | 0.719 (71.9%) | [0.630, 0.800] | U | Messner et al. (2014) |
| Dose–response infection model (confluent hypergeometric) | — | Exact Beta-Poisson ${}_1F_1$ series expansion | Structural model | — | Model | Teunis et al. (2020) |
| Beta shape parameter for infection (dimensionless) | $\alpha$ | Bivariate log-normal MCMC joint posterior draw | 0.233 | [0.012, 4.85] | U | Teunis et al. (2020) |
| Beta scale parameter for infection (dimensionless) | $\beta$ | Bivariate log-normal MCMC joint posterior draw | 5.90 | [0.101, 251] | U | Teunis et al. (2020) |
| Conditional illness-given-infection model | — | Negative binomial conditional dose-dependent illness model | Structural model | — | Model | Teunis et al. (2020) |
| Illness scale parameter (virions^-1) | $\eta$ | Bivariate log-normal MCMC joint posterior draw | 2.02 | [0.009, 279] | U | Teunis et al. (2020) |
| Illness shape parameter (dimensionless) | $r$ | Bivariate log-normal MCMC joint posterior draw | 1.58 | [0.001, 416] | U | Teunis et al. (2020) |
| 6. Risk Characterization & Annual Disease Burden (Fig. 1: State 6)<br>*Model: $E = N_{\text{pop}} \cdot (52 \cdot p_{\text{freq}})$; $I_{\text{ann}} = E \cdot \mathbb{E}_{\text{var}}[P_{\text{ill}}]$; 2D Monte Carlo ($N_{\text{var}} \times N_{\text{unc}}$)* | | | | | | |
| Raw oyster weekly consumption probability | $p_{\text{raw}}$ | Empirical 7-day food recall proportion (FoodNet survey) | 0.007 (0.7%) | — | Empirical | CDC FoodNet (2006–2007) |
| Cooked oyster weekly consumption probability | $p_{\text{cooked}}$ | Empirical 7-day recall proportion (0.026 total - 0.007 raw) | 0.019 (1.9%) | — | Empirical | CDC FoodNet (2006–2007) |
| Total target U.S. consumer population (persons) | $N_{\text{pop}}$ | 2024 U.S. Census Bureau resident population projection | 341,800,000 | — | Empirical | U.S. Census Bureau (2024) |
| Variability iterations (individual servings) | $N_{\text{var}}$ | Monte Carlo sample size for dose & serving variability | 10,000 | — | Simulation | This study |
| Uncertainty iterations (parameter sets) | $N_{\text{unc}}$ | Second-order 2D Monte Carlo uncertainty sample size | 2,000 | — | Simulation | This study |

---
*Notes*:

- Full mathematical models are written as explicit formula equations in section headers and model rows. Each parameter is assigned a dedicated row without bold styling or bullet points.
- Background is pure white; classic APA 3-line table borders are applied.
- Dimensions: U = uncertainty dimension; V = variability dimension; Scenario = discrete predefined regimes; Empirical = surveillance/experimental histories; Simulation = Monte Carlo iterations.
- References: Nishida et al. (2003, 2007); Mirmahdi et al. (2025); Nalanda et al. (2026); Dicker et al. (2026); Wu et al. (2026); Messner et al. (2014); Teunis et al. (2020); CDC FoodNet (2006–2007); U.S. Census Bureau (2024).