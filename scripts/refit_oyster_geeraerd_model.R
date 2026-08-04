# Module 2: Dynamic Geeraerd thermal inactivation scenario analysis

# 0. Setup

library(tidyverse)
library(here)

project_root <- here::here()
data_dir <- file.path(project_root, "data", "module_03_thermal_inactivation_oyster")
output_dir <- file.path(project_root, "results", "module_03_thermal_inactivation_oyster")

source(file.path(project_root, "scripts", "module3_thermal_functions.R"))

# 1. Load Sammi's original data
dyn_path <- file.path(data_dir, "d_oys_pfu.RDS")
comeup_path <- file.path(data_dir, "df_comeup.RDS")
df_dyn <- readRDS(dyn_path) %>% tibble::as_tibble()
df_comeup <- readRDS(comeup_path) %>% tibble::as_tibble()
dyn_names <- colnames(df_dyn)
comeup_names <- colnames(df_comeup)



# 2. Add new dynamic heating replicate data

new_raw <- tibble::tribble(
  ~replicate, ~tissue, ~virus, ~assay, ~target_C, ~time_min, ~temp_C, ~log_value,
  2, "dd", "tv", "pfu", 100, 0, 7.5,  3.228,
  2, "dd", "tv", "pfu", 100, 0.5, 64.7, 2.990,
  2, "dd", "tv", "pfu", 100, 1, 69.1, 1.179,
  2, "dd", "tv", "pfu", 100, 1.25, 71.9, 0.250,
  2, "dd", "tv", "pfu", 100, 1.5, 75.5, 1.329,
  2, "dd", "tv", "pfu", 100, 1.833333333, 82.7, 1.863,
  3, "dd", "tv", "pfu", 100, 0, 8.4,  2.727,
  3, "dd", "tv", "pfu", 100, 0.5, 64.8, 2.505,
  3, "dd", "tv", "pfu", 100, 1, 69.7, 0.551,
  3, "dd", "tv", "pfu", 100, 1.25, 72.8, 0.250,
  3, "dd", "tv", "pfu", 100, 1.5, 76.5, 1.426,
  3, "dd", "tv", "pfu", 100, 1.833333333, 84.4, 0.648
)


# 3. integrate the new data into existing tibble

# Match df_dyn structure
new_dyn <- new_raw %>%
  mutate(
    cond = paste0(assay, "-", target_C, "-", tissue),
    time = time_min,
    temperature = temp_C,
    logN = log_value,
    rep = replicate
  ) %>%
  select(all_of(dyn_names))

df_dyn_updated <- bind_rows(df_dyn, new_dyn) %>%
  prepare_observed_dynamic_profiles()

saveRDS(df_dyn_updated, file.path(data_dir, "df_dyn_updated.rds"))

# Match df_comeup structure
new_comeup <- new_raw %>%
  mutate(rep = replicate,
         matrix = tissue,
         target_temp = target_C,
         time = time_min,
         sample_temp = temp_C) %>%
  select(all_of(comeup_names)) %>%
  distinct()

df_comeup_updated <- bind_rows(df_comeup, new_comeup)

saveRDS(df_comeup_updated, file.path(data_dir, "df_comeup_updated.rds"))

# 4. Observed dynamic heating model helpers are sourced from
# scripts/module3_thermal_functions.R. They preserve replicate-specific
# time-temperature profiles while fitting one shared global Geeraerd parameter
# set across all observed profiles.

# 5. refit global dynamic Geeraerd model

known_ge <- c(temp_ref = 57, C_c0 = 1e-3)
guess_ge <- c(D_R = 10, z = 10, logN0 = 3.6, logNmin = 2)
lower_ge <- c(D_R = 0.1, z = 2, logN0 = 3, logNmin = 0)
upper_ge <- c(D_R = 10, z = 15, logN0 = 6, logNmin = 3)

set.seed(7938)

fit_ge <- fit_global_dynamic_geeraerd(
  d = df_dyn_updated,
  guess = guess_ge,
  known = known_ge,
  lower = lower_ge,
  upper = upper_ge,
  niter = 5000
)

saveRDS(fit_ge, file.path(output_dir, "module3_dynamic_geeraerd_fit.rds"))

fit_ge <- readRDS(file.path(output_dir, "module3_dynamic_geeraerd_fit.rds"))

# 6. Extract posterior parameter draws
post_ge <- as.data.frame(fit_ge$pars) %>%
  transmute(
    D_R = D_R,
    z = z,
    logN0 = logN0,
    logNmin = logNmin,
    N0 = 10^logN0,
    N_min = 10^logNmin
  )

saveRDS(post_ge, file.path(output_dir, "module3_geeraerd_posterior_draws.rds"))

# 7. Posterior model fit check: fitted Geeraerd curves over observed data
geeraerd_fit_overlay <- make_geeraerd_observed_fit_overlay(
  posterior_draws = post_ge,
  observed_data = df_dyn_updated,
  known = known_ge
)

write_csv(
  geeraerd_fit_overlay$prediction_summary,
  file.path(output_dir, "module3_geeraerd_observed_vs_fitted_curve_summary.csv")
)

ggsave(
  file.path(output_dir, "module3_geeraerd_observed_vs_fitted_overlay.png"),
  plot = geeraerd_fit_overlay$plot,
  width = 10,
  height = 12,
  dpi = 300,
  bg = "white"
)
