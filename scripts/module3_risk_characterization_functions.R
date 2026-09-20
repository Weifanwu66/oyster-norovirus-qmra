# Module 4: Hazard & Risk Characterization Function Library
# Contains Teunis et al. (2020) Dose-Response model functions, Host Susceptibility Beta fitting,
# Serving-Level Risk Evaluation with Fully Random Multinomial Size Mixing, Heatmap Visualization,
# Tornado Sensitivity Analysis, and U.S. Population Illness Burden Calculations using FoodNet Data.

library(tidyverse)
library(gsl)
library(MASS)

# 1. Host Susceptibility Beta Parameter Back-Transformation (Messner et al. 2014)

messner_p_mean <- 0.72
messner_p_q025 <- 0.63
messner_p_q975 <- 0.80

#' Fit Beta distribution parameters (alpha, beta) matching a specified mean and 95% CrI
fit_beta_from_mean_ci <- function(mean_value = 0.72, q025 = 0.63, q975 = 0.80) {
  objective <- function(log_par) {
    alpha <- exp(log_par[1])
    beta  <- exp(log_par[2])

    fitted_mean <- alpha / (alpha + beta)
    fitted_q025 <- qbeta(0.025, shape1 = alpha, shape2 = beta)
    fitted_q975 <- qbeta(0.975, shape1 = alpha, shape2 = beta)

    (fitted_mean - mean_value)^2 +
      (fitted_q025 - q025)^2 +
      (fitted_q975 - q975)^2
  }

  fit <- optim(
    par = c(log(75), log(30)),
    fn = objective
  )

  alpha <- exp(fit$par[1])
  beta <- exp(fit$par[2])

  tibble::tibble(
    alpha = alpha,
    beta = beta,
    fitted_mean = alpha / (alpha + beta),
    fitted_q025 = qbeta(0.025, shape1 = alpha, shape2 = beta),
    fitted_q975 = qbeta(0.975, shape1 = alpha, shape2 = beta)
  )
}

# 2. Teunis et al. (2020) Table A14 Dose-Response Parameter Generator & Core Math

#' Sample Teunis et al. (2020) Table A14 Dose-Response MCMC parameter draws
sample_teunis_2020_draws <- function(n_unc = 1000, seed = 1060) {
  set.seed(seed)
  mu_inf <- c(-3.16, 1.82)
  cov_inf <- matrix(c(2.88, -2.01, -2.01, 3.64), ncol = 2)

  mu_ill <- c(0.414, 1.73)
  cov_ill <- matrix(c(7.76, -1.67, -1.67, 6.19), ncol = 2)

  m_inf <- MASS::mvrnorm(n = n_unc, mu = mu_inf, Sigma = cov_inf)
  m_ill <- MASS::mvrnorm(n = n_unc, mu = mu_ill, Sigma = cov_ill)

  tibble::tibble(
    w_eta = m_inf[, 1],
    z_eta = m_inf[, 2],
    w_mu  = m_ill[, 1],
    z_mu  = m_ill[, 2]
  )
}

par1.btr <- function(w, z) {
  exp(z) / (1 + 1 / exp(w))
}

par2.btr <- function(w, z) {
  exp(z) / (exp(w) + 1)
}

#' Teunis 1F1 Kummer Hypergeometric Infection Probability
p.inf.hyp1f1 <- function(a, b, D) {
  res <- numeric(length(D))
  pos_idx <- D > 1e-12
  if (any(pos_idx)) {
    a_pos <- if (length(a) > 1) a[pos_idx] else a
    b_pos <- if (length(b) > 1) b[pos_idx] else b
    val <- gsl::hyperg_1F1(a_pos, a_pos + b_pos, -D[pos_idx])
    val[is.nan(val)] <- 0
    res[pos_idx] <- 1 - val
  }
  return(res)
}

#' Teunis Conditional Illness Probability given Infection
p.ill.inf <- function(n, r, D) {
  1 - (1 + D / n)^(-r)
}

# 3. Comprehensive Teunis Risk Evaluation Functions

evaluate_teunis_risk_matrix <- function(
  dose_matrix,
  teunis_draws = NULL,
  se_alpha = NULL,
  se_beta = NULL
) {
  n_var <- nrow(dose_matrix)
  n_unc <- ncol(dose_matrix)
  
  if (is.null(teunis_draws)) {
    teunis_draws <- sample_teunis_2020_draws(n_unc = n_unc)
  }
  
  if (is.null(se_alpha) || is.null(se_beta)) {
    se_fit <- fit_beta_from_mean_ci()
    se_alpha <- se_fit$alpha
    se_beta  <- se_fit$beta
  }
  
  set.seed(2020)
  p_se_u <- rbeta(n_unc, shape1 = se_alpha, shape2 = se_beta)

  a_u <- par1.btr(teunis_draws$w_eta, teunis_draws$z_eta)
  b_u <- par2.btr(teunis_draws$w_eta, teunis_draws$z_eta)
  r_u <- par1.btr(teunis_draws$w_mu, teunis_draws$z_mu)
  n_u <- par2.btr(teunis_draws$w_mu, teunis_draws$z_mu)
  
  p_inf_matrix <- matrix(0, nrow = n_var, ncol = n_unc)
  p_ill_matrix <- matrix(0, nrow = n_var, ncol = n_unc)
  
  for (u in seq_len(n_unc)) {
    d_vec <- dose_matrix[, u]
    p_inf_se <- p.inf.hyp1f1(a_u[u], b_u[u], d_vec)
    p_ill_given_inf <- p.ill.inf(n_u[u], r_u[u], d_vec)
    
    p_inf_matrix[, u] <- p_se_u[u] * p_inf_se
    p_ill_matrix[, u] <- p_se_u[u] * p_inf_se * p_ill_given_inf
  }
  
  list(
    mean_p_inf = mean(p_inf_matrix),
    mean_p_ill = mean(p_ill_matrix),
    median_p_ill = median(p_ill_matrix),
    p05_p_ill = quantile(p_ill_matrix, 0.05),
    p50_p_ill = median(p_ill_matrix),
    p95_p_ill = quantile(p_ill_matrix, 0.95)
  )
}

#' High-level function to evaluate Module 4 risk across single-size and mixed-size conditions
evaluate_model_risk_summary <- function(
  postcook_nodes_list,
  condition_grid,
  n_var = 10000,
  n_unc = 2000,
  teunis_draws = NULL,
  se_fit = NULL
) {
  if (is.null(teunis_draws)) {
    teunis_draws <- sample_teunis_2020_draws(n_unc = n_unc)
  }
  if (is.null(se_fit)) {
    se_fit <- fit_beta_from_mean_ci()
  }
  
  # 1. Single Size Scenarios (Small, Medium, Large)
  single_size_res <- pmap_dfr(condition_grid, function(method, size, time_min) {
    scen_key <- paste(method, size, paste0(time_min, "min"), sep = "_")
    per_oyster_dose_node <- postcook_nodes_list[[scen_key]]
    node_mat <- unclass(per_oyster_dose_node)[, , 1]
    
    map_dfr(c(6, 12), function(serving_sz) {
      serving_dose_matrix <- matrix(0, nrow = n_var, ncol = n_unc)
      
      for (u in seq_len(n_unc)) {
        col_doses <- node_mat[, u]
        # Fast vectorized matrix sampling for serving size
        sampled_mat <- matrix(sample(col_doses, size = n_var * serving_sz, replace = TRUE), nrow = n_var, ncol = serving_sz)
        serving_dose_matrix[, u] <- rowSums(sampled_mat)
      }
      
      risk_res <- evaluate_teunis_risk_matrix(
        dose_matrix = serving_dose_matrix,
        teunis_draws = teunis_draws,
        se_alpha = se_fit$alpha,
        se_beta = se_fit$beta
      )
      
      res_row <- tibble(
        method = method,
        size = size,
        time_min = time_min,
        serving_size = serving_sz,
        mean_infectious_dose_per_serving = mean(serving_dose_matrix),
        median_infectious_dose_per_serving = median(serving_dose_matrix),
        mean_p_infection = risk_res$mean_p_inf,
        mean_p_illness = risk_res$mean_p_ill,
        median_p_illness = risk_res$median_p_ill,
        p025_p_illness = risk_res$p05_p_ill,
        p975_p_illness = risk_res$p95_p_ill
      )
      
      rm(serving_dose_matrix, risk_res)
      gc(verbose = FALSE)
      
      res_row
    })
  })

  # 2. Mixed Size Scenarios (Fully Random Multinomial Assign: P(S)=1/3, P(M)=1/3, P(L)=1/3)
  unique_methods_times <- condition_grid %>% dplyr::distinct(method, time_min)
  
  mixed_size_res <- pmap_dfr(unique_methods_times, function(method, time_min) {
    key_small  <- paste(method, "small",  paste0(time_min, "min"), sep = "_")
    key_medium <- paste(method, "medium", paste0(time_min, "min"), sep = "_")
    key_large  <- paste(method, "large",  paste0(time_min, "min"), sep = "_")
    
    mat_s <- unclass(postcook_nodes_list[[key_small]])[, , 1]
    mat_m <- unclass(postcook_nodes_list[[key_medium]])[, , 1]
    mat_l <- unclass(postcook_nodes_list[[key_large]])[, , 1]
    
    map_dfr(c(6, 12), function(serving_sz) {
      serving_dose_matrix <- matrix(0, nrow = n_var, ncol = n_unc)
      
      set.seed(3030 + serving_sz)
      for (u in seq_len(n_unc)) {
        col_s <- mat_s[, u]
        col_m <- mat_m[, u]
        col_l <- mat_l[, u]
        
        # Fast vectorized matrix sampling for multinomial mixed sizes per serving
        sz_choices <- matrix(sample.int(3, size = n_var * serving_sz, replace = TRUE), nrow = n_var, ncol = serving_sz)
        
        serving_mat <- matrix(0, nrow = n_var, ncol = serving_sz)
        for (k in seq_len(serving_sz)) {
          c_k <- sz_choices[, k]
          idx1 <- c_k == 1
          idx2 <- c_k == 2
          idx3 <- c_k == 3
          
          d_k <- numeric(n_var)
          if (any(idx1)) d_k[idx1] <- sample(col_s, size = sum(idx1), replace = TRUE)
          if (any(idx2)) d_k[idx2] <- sample(col_m, size = sum(idx2), replace = TRUE)
          if (any(idx3)) d_k[idx3] <- sample(col_l, size = sum(idx3), replace = TRUE)
          
          serving_mat[, k] <- d_k
        }
        
        serving_dose_matrix[, u] <- rowSums(serving_mat)
      }
      
      risk_res <- evaluate_teunis_risk_matrix(
        dose_matrix = serving_dose_matrix,
        teunis_draws = teunis_draws,
        se_alpha = se_fit$alpha,
        se_beta = se_fit$beta
      )
      
      res_row <- tibble(
        method = method,
        size = "mixed",
        time_min = time_min,
        serving_size = serving_sz,
        mean_infectious_dose_per_serving = mean(serving_dose_matrix),
        median_infectious_dose_per_serving = median(serving_dose_matrix),
        mean_p_infection = risk_res$mean_p_inf,
        mean_p_illness = risk_res$mean_p_ill,
        median_p_illness = risk_res$median_p_ill,
        p025_p_illness = risk_res$p05_p_ill,
        p975_p_illness = risk_res$p95_p_ill
      )
      
      rm(serving_dose_matrix, risk_res)
      gc(verbose = FALSE)
      
      res_row
    })
  })

  bind_rows(single_size_res, mixed_size_res)
}

# 4. Risk Formatting & Heatmap Visualization Functions

format_risk_label <- function(x) {
  out <- ifelse(
    x >= 0.001,
    formatC(x, format = "f", digits = 3),
    formatC(x, format = "e", digits = 0)
  )
  out <- sub("e-0", "e-", out)
  sub("e\\+0", "e+", out)
}

format_count_label <- function(x) {
  formatC(x, format = "f", digits = 0, big.mark = ",")
}

format_risk_interval <- function(median_value, lower_value, upper_value) {
  paste0(
    format_risk_label(median_value),
    " [",
    format_risk_label(lower_value),
    ", ",
    format_risk_label(upper_value),
    "]"
  )
}


plot_module4_risk_heatmap <- function(summary_df, output_filepath = NULL, risk_metric = "mean_p_illness", threshold = 1e-4) {
  heatmap_data <- summary_df %>%
    mutate(
      method_label = str_to_title(method),
      size_label = factor(str_to_title(size), levels = c("Small", "Medium", "Large", "Mixed")),
      serving_label = paste0(serving_size, " Oysters / Serving"),
      plot_risk = pmax(.data[[risk_metric]], 1e-6),
      below_threshold = plot_risk < threshold
    )

  metric_title <- if (risk_metric == "mean_p_illness") "Mean Population P(illness)" else "Median P(illness)"

  p <- ggplot(
    heatmap_data,
    aes(
      x = factor(time_min),
      y = size_label,
      fill = plot_risk
    )
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_point(
      data = dplyr::filter(heatmap_data, below_threshold),
      shape = 21,
      color = "grey20",
      fill = NA,
      stroke = 0.55,
      size = 2.2
    ) +
    facet_grid(serving_label ~ method_label, scales = "free_x", space = "free_x") +
    scale_fill_gradientn(
      colors = c("#F7FBFF", "#D6EAF8", "#9ECAE1", "#FDD49E", "#FC8D59", "#B30000"),
      trans = "log10",
      limits = c(1e-6, 1e-1),
      breaks = 10^(-6:-1),
      labels = scales::math_format(10^.x),
      name = metric_title,
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(8, "cm"), barheight = unit(0.4, "cm"))
    ) +
    labs(
      x = "Cooking holding time (minutes)",
      y = "Oyster Size Category"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey55"),
      strip.text = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 12)
    )

  if (!is.null(output_filepath)) ggsave(filename = output_filepath, plot = p, width = 10, height = 7, dpi = 300, compression = "lzw")

  p
}

# 5. Global Spearman Rank Correlation Tornado Plot Function
#' Plot Tornado Sensitivity Diagram for Norovirus Illness Risk P(illness)
plot_tornado_sensitivity <- function(
  postcook_nodes_list,
  output_filepath = NULL,
  teunis_draws = NULL,
  se_alpha = NULL,
  se_beta = NULL,
  n_sample = 5000,
  model_type = c("geeraerd", "meta")
) {
  model_type <- match.arg(model_type)
  set.seed(4040)
  
  if (is.null(teunis_draws)) {
    teunis_draws <- sample_teunis_2020_draws(n_unc = 500)
  }
  
  if (is.null(se_alpha) || is.null(se_beta)) {
    se_fit <- fit_beta_from_mean_ci()
    se_alpha <- se_fit$alpha
    se_beta  <- se_fit$beta
  }
  
  scen_key <- grep("steaming_medium_4min", names(postcook_nodes_list), value = TRUE)[1]
  if (is.na(scen_key)) scen_key <- names(postcook_nodes_list)[1]
  
  post_node <- postcook_nodes_list[[scen_key]]
  post_mat <- unclass(post_node)[, , 1]
  
  n_v <- nrow(post_mat)
  n_u <- ncol(post_mat)
  
  idx_v <- sample.int(n_v, n_sample, replace = TRUE)
  idx_u <- sample.int(n_u, n_sample, replace = TRUE)
  
  postcook_val <- sapply(seq_len(n_sample), function(i) post_mat[idx_v[i], idx_u[i]])
  p_se_val <- rbeta(n_sample, shape1 = se_alpha, shape2 = se_beta)
  
  csv_path <- here::here("results", "module1_3_exposure_dose_calculation", "module1_uncertainty_inputs.csv")
  if (file.exists(csv_path)) {
    module1_u <- read_csv(csv_path, show_col_types = FALSE)
    q_inv_val <- module1_u$q_inv[idx_u]
    prev_mean_val <- module1_u$prevalence_mean[idx_u]
    conc_mean_val <- module1_u$positive_mean_gc[idx_u]
    conc_shape_val <- module1_u$concentration_gamma_shape[idx_u]
  } else {
    stop("module1_uncertainty_inputs.csv not found!")
  }
  
  teunis_sel <- teunis_draws[sample.int(nrow(teunis_draws), n_sample, replace = TRUE), ]
  a_val <- par1.btr(teunis_sel$w_eta, teunis_sel$z_eta)
  b_val <- par2.btr(teunis_sel$w_eta, teunis_sel$z_eta)
  r_val <- par1.btr(teunis_sel$w_mu, teunis_sel$z_mu)
  n_val <- par2.btr(teunis_sel$w_mu, teunis_sel$z_mu)
  
  dose_effective <- postcook_val
  p_inf_se <- p.inf.hyp1f1(a_val, b_val, dose_effective)
  p_ill_given_inf <- p.ill.inf(n_val, r_val, dose_effective)
  p_illness <- p_se_val * p_inf_se * p_ill_given_inf
  
  input_data <- tibble(
    `Harvest Prevalence (P_prev)` = prev_mean_val,
    `Positive Mean Concentration (μ_conc)` = conc_mean_val,
    `Positive Concentration Shape (κ_conc)` = conc_shape_val,
    `RT-qPCR: Infectious Virus Ratio (Q)` = 1 / q_inv_val,
    `Host Secretor Susceptibility P(Se+)` = p_se_val,
    `DR α (Infection)` = a_val,
    `DR β (Infection)` = b_val,
    `DR η (Illness)` = n_val,
    `DR ρ (Illness)` = r_val
  )
  
  if (model_type == "geeraerd") {
    therm_post_raw <- readRDS(here::here("results", "module1_2_thermal_inactivation", "module1_2_geeraerd_thermal_kinetics_posterior.rds"))
    set.seed(1503)
    therm_post <- therm_post_raw %>% dplyr::slice_sample(n = n_u, replace = TRUE)
    input_data <- input_data %>% mutate(
      `Geeraerd Thermal Sensitivity (z-value)` = therm_post$z[idx_u],
      `Geeraerd Reference Inactivation Time (D_r)` = therm_post$D_R[idx_u],
      `Geeraerd Residual Tail (log10 N_min)` = therm_post$logNmin[idx_u],
      `Geeraerd Initial Shoulder (log10 N_0)` = therm_post$logN0[idx_u]
    )
    title_text <- "Dynamic Geeraerd Inactivation Model"
  } else {
    therm_post_raw <- readRDS(here::here("results", "module1_2_thermal_inactivation", "module1_2_meta_analysis_thermal_kinetics_posterior.rds"))
    set.seed(1503)
    therm_post <- therm_post_raw %>% dplyr::slice_sample(n = n_u, replace = TRUE)
    input_data <- input_data %>% mutate(
      `Mafart Thermal Sensitivity (z-value)` = therm_post$z_shellfish[idx_u],
      `Mafart Reference Inactivation Time (D_r)` = therm_post$Dref_shell[idx_u],
      `Mafart Shape Parameter (p)` = therm_post$p_weibull[idx_u]
    )
    title_text <- "Dynamic Meta-Analysis Inactivation Model"
  }
  
  cor_df <- map_dfr(names(input_data), function(var_name) {
    rho <- cor(input_data[[var_name]], p_illness, method = "spearman", use = "pairwise.complete.obs")
    tibble(
      input = var_name,
      spearman_rho = rho,
      abs_rho = abs(rho),
      direction = if_else(rho >= 0, "Positive Association", "Negative Association")
    )
  }) %>%
    arrange(desc(abs_rho)) %>%
    filter(!is.na(spearman_rho)) %>%
    mutate(input = fct_reorder(input, abs_rho))
  
  p <- ggplot(cor_df, aes(x = spearman_rho, y = input, fill = direction)) +
    geom_col(width = 0.65) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("Positive Association" = "#2F6F9F", "Negative Association" = "#B65C45")) +
    scale_x_continuous(limits = c(-0.45, 0.45), breaks = seq(-0.4, 0.4, by = 0.1)) +
    labs(
      title = title_text,
      x = "Spearman Rank Correlation Coefficient (p)",
      y = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    )
  
  if (!is.null(output_filepath)) ggsave(filename = output_filepath, plot = p, width = 8, height = 4.8, dpi = 300, compression = "lzw")
  
  p
}

# 6. Monte Carlo Filtering (MCF) & Regional Sensitivity Analysis (Saltelli et al. 2004, FDA Pouillot et al. 2022)

#' Perform Monte Carlo Filtering (MCF) using Cutoff 1e-4 (0.01% illness risk)
run_monte_carlo_filtering <- function(
  postcook_nodes_list,
  teunis_draws = NULL,
  se_alpha = NULL,
  se_beta = NULL,
  cutoff = 0.0001,
  n_sample = 10000,
  model_type = c("geeraerd", "meta")
) {
  model_type <- match.arg(model_type)
  set.seed(5050)
  
  if (is.null(teunis_draws)) {
    teunis_draws <- sample_teunis_2020_draws(n_unc = 500)
  }
  
  if (is.null(se_alpha) || is.null(se_beta)) {
    se_fit <- fit_beta_from_mean_ci()
    se_alpha <- se_fit$alpha
    se_beta  <- se_fit$beta
  }
  
  scen_key <- grep("steaming_medium_4min", names(postcook_nodes_list), value = TRUE)[1]
  if (is.na(scen_key)) scen_key <- names(postcook_nodes_list)[1]
  
  post_node <- postcook_nodes_list[[scen_key]]
  post_mat <- unclass(post_node)[, , 1]
  
  n_v <- nrow(post_mat)
  n_u <- ncol(post_mat)
  
  idx_v <- sample.int(n_v, n_sample, replace = TRUE)
  idx_u <- sample.int(n_u, n_sample, replace = TRUE)
  
  postcook_val <- sapply(seq_len(n_sample), function(i) post_mat[idx_v[i], idx_u[i]])
  p_se_val <- rbeta(n_sample, shape1 = se_alpha, shape2 = se_beta)
  
  csv_path <- here::here("results", "module1_3_exposure_dose_calculation", "module1_uncertainty_inputs.csv")
  if (file.exists(csv_path)) {
    module1_u <- read_csv(csv_path, show_col_types = FALSE)
    q_inv_val <- module1_u$q_inv[idx_u]
    prev_mean_val <- module1_u$prevalence_mean[idx_u]
    conc_mean_val <- module1_u$positive_mean_gc[idx_u]
    conc_shape_val <- module1_u$concentration_gamma_shape[idx_u]
  } else {
    stop("module1_uncertainty_inputs.csv not found!")
  }
  
  teunis_sel <- teunis_draws[sample.int(nrow(teunis_draws), n_sample, replace = TRUE), ]
  a_val <- par1.btr(teunis_sel$w_eta, teunis_sel$z_eta)
  b_val <- par2.btr(teunis_sel$w_eta, teunis_sel$z_eta)
  r_val <- par1.btr(teunis_sel$w_mu, teunis_sel$z_mu)
  n_val <- par2.btr(teunis_sel$w_mu, teunis_sel$z_mu)
  
  dose_effective <- postcook_val
  p_inf_se <- p.inf.hyp1f1(a_val, b_val, dose_effective)
  p_ill_given_inf <- p.ill.inf(n_val, r_val, dose_effective)
  p_illness <- p_se_val * p_inf_se * p_ill_given_inf
  
  is_unsafe <- p_illness >= cutoff
  
  input_data <- tibble(
    `Harvest Prevalence (P_prev)` = prev_mean_val,
    `Positive Mean Concentration (μ_conc)` = conc_mean_val,
    `Positive Concentration Shape (κ_conc)` = conc_shape_val,
    `RT-qPCR: Infectious Virus Ratio (Q)` = 1 / q_inv_val,
    `Host Secretor Susceptibility P(Se+)` = p_se_val,
    `DR α (Infection)` = a_val,
    `DR β (Infection)` = b_val,
    `DR η (Illness)` = n_val,
    `DR ρ (Illness)` = r_val
  )
  
  if (model_type == "geeraerd") {
    therm_post_raw <- readRDS(here::here("results", "module1_2_thermal_inactivation", "module1_2_geeraerd_thermal_kinetics_posterior.rds"))
    set.seed(1503)
    therm_post <- therm_post_raw %>% dplyr::slice_sample(n = n_u, replace = TRUE)
    input_data <- input_data %>% mutate(
      `Geeraerd Thermal Sensitivity (z-value)` = therm_post$z[idx_u],
      `Geeraerd Reference Inactivation Time (D_r)` = therm_post$D_R[idx_u],
      `Geeraerd Residual Tail (log10 N_min)` = therm_post$logNmin[idx_u],
      `Geeraerd Initial Shoulder (log10 N_0)` = therm_post$logN0[idx_u]
    )
  } else {
    therm_post_raw <- readRDS(here::here("results", "module1_2_thermal_inactivation", "module1_2_meta_analysis_thermal_kinetics_posterior.rds"))
    set.seed(1503)
    therm_post <- therm_post_raw %>% dplyr::slice_sample(n = n_u, replace = TRUE)
    input_data <- input_data %>% mutate(
      `Mafart Thermal Sensitivity (z-value)` = therm_post$z_shellfish[idx_u],
      `Mafart Reference Inactivation Time (D_r)` = therm_post$Dref_shell[idx_u],
      `Mafart Shape Parameter (p)` = therm_post$p_weibull[idx_u]
    )
  }
  
  ks_table <- map_dfr(names(input_data), function(param_name) {
    val_unsafe <- na.omit(input_data[[param_name]][is_unsafe])
    val_safe   <- na.omit(input_data[[param_name]][!is_unsafe])
    
    if (length(val_unsafe) == 0 || length(val_safe) == 0) return(NULL)
    
    ks_res <- ks.test(val_unsafe, val_safe)
    tibble(
      Parameter = param_name,
      `KS_D_statistic` = as.numeric(ks_res$statistic),
      `p_value` = ks_res$p.value,
      n_unsafe = length(val_unsafe),
      n_safe = length(val_safe)
    )
  }) %>%
    mutate(`p_value_BH_FDR` = p.adjust(p_value, method = "BH")) %>%
    arrange(desc(KS_D_statistic))
  
  list(
    input_data = input_data,
    is_unsafe = is_unsafe,
    ks_table = ks_table
  )
}

#' Plot Monte Carlo Filtering Cumulative Distribution Function (CDF) Comparison Curves
plot_mcf_cdf_curves <- function(mcf_result, output_filepath = NULL) {
  top_param <- mcf_result$ks_table$Parameter[1]
  
  unsafe_label <- sprintf("Unsafe / High-Risk (P >= %g)", mcf_result$cutoff)
  safe_label   <- sprintf("Safe / Compliant (P < %g)", mcf_result$cutoff)

  df_plot <- mcf_result$sample_df %>%
    mutate(
      Group = if_else(is_unsafe, unsafe_label, safe_label)
    )
  
  col_scale <- setNames(c("#D95F02", "#1B9E77"), c(unsafe_label, safe_label))

  p <- ggplot(df_plot, aes(x = .data[[top_param]], color = Group)) +
    stat_ecdf(geom = "step", linewidth = 1.1) +
    scale_color_manual(values = col_scale) +
    labs(
      x = top_param,
      y = "Empirical Cumulative Distribution F(x)",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  
  if (!is.null(output_filepath)) ggsave(filename = output_filepath, plot = p, width = 8, height = 4.8, dpi = 300, compression = "lzw")
  
  p
}

# 7. 2D Monte Carlo Convergence Diagnostic Functions (Ballio & Guadagnini 2004, FDA Guidelines)

#' Assess numerical convergence of Mean P(illness) across Variability (N_var) and Uncertainty (N_unc)
assess_simulation_convergence <- function(
  precook_node,
  sizes_var = c(500, 1000, 2500, 5000, 10000),
  sizes_unc = c(50, 100, 250, 500, 1000, 1500, 2000)
) {
  mat <- unclass(precook_node)[, , 1]
  n_v_total <- nrow(mat)
  n_u_total <- ncol(mat)
  
  teunis_draws <- sample_teunis_2020_draws(n_unc = n_u_total)
  se_fit <- fit_beta_from_mean_ci()

  # 1. Variability Convergence (holding N_unc = max)
  var_conv <- map_dfr(sizes_var[sizes_var <= n_v_total], function(sz_v) {
    sub_mat <- mat[1:sz_v, ]
    risk_res <- evaluate_teunis_risk_matrix(
      dose_matrix = sub_mat,
      teunis_draws = teunis_draws,
      se_alpha = se_fit$alpha,
      se_beta = se_fit$beta
    )
    tibble(
      Dimension = "Variability",
      SampleSize = sz_v,
      Mean_P_illness = risk_res$mean_p_ill,
      P95_P_illness = risk_res$p95_p_ill
    )
  })

  # 2. Uncertainty Convergence (holding N_var = max)
  unc_conv <- map_dfr(sizes_unc[sizes_unc <= n_u_total], function(sz_u) {
    sub_mat <- mat[, 1:sz_u]
    sub_teunis <- teunis_draws[1:sz_u, ]
    risk_res <- evaluate_teunis_risk_matrix(
      dose_matrix = sub_mat,
      teunis_draws = sub_teunis,
      se_alpha = se_fit$alpha,
      se_beta = se_fit$beta
    )
    tibble(
      Dimension = "Uncertainty",
      SampleSize = sz_u,
      Mean_P_illness = risk_res$mean_p_ill,
      P95_P_illness = risk_res$p95_p_ill
    )
  })

  bind_rows(var_conv, unc_conv)
}

#' Plot 2D Monte Carlo Numerical Convergence Curves
plot_simulation_convergence <- function(conv_df, output_filepath = NULL) {
  p <- ggplot(conv_df, aes(x = SampleSize, y = Mean_P_illness, color = Dimension)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    facet_wrap(~ Dimension, scales = "free_x") +
    scale_color_manual(values = c("Variability" = "#2980B9", "Uncertainty" = "#C0392B")) +
    labs(
      x = "Iteration",
      y = "Mean HuNoV illness risk P(illness)",
      color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey92", color = "grey55"),
      strip.text = element_text(face = "bold", size = 11)
    )
  if (!is.null(output_filepath)) ggsave(filename = output_filepath, plot = p, width = 8, height = 4, dpi = 300, compression = "lzw")
  p
}

#' Run Monte Carlo Filtering Loop across all cooking scenarios and size classes
run_mcf_scenario_matrix <- function(
  postcook_nodes_list,
  cutoff = 0.0001,
  n_sample = 10000,
  model_type = c("geeraerd", "meta")
) {
  model_type <- match.arg(model_type)
  scen_names <- names(postcook_nodes_list)
  
  col_unsafe <- sprintf("Unsafe Iterations (P >= %g)", cutoff)
  col_safe   <- sprintf("Safe Iterations (P < %g)", cutoff)

  map_dfr(scen_names, function(scen_key) {
    post_node <- postcook_nodes_list[[scen_key]]
    
    mcf_res <- run_monte_carlo_filtering(
      postcook_nodes_list = setNames(list(post_node), scen_key),
      cutoff = cutoff,
      n_sample = n_sample,
      model_type = model_type
    )
    
    parts <- strsplit(scen_key, "_")[[1]]
    method <- str_to_title(parts[1])
    size <- str_to_title(parts[2])
    time_min <- as.numeric(gsub("min", "", parts[3]))
    
    n_unsafe <- sum(mcf_res$is_unsafe)
    n_safe <- sum(!mcf_res$is_unsafe)
    
    top_param <- if (n_unsafe >= 5 && nrow(mcf_res$ks_table) > 0) mcf_res$ks_table$Parameter[1] else "None (Fully Safe)"
    top_d     <- if (n_unsafe >= 5 && nrow(mcf_res$ks_table) > 0) mcf_res$ks_table$KS_D_statistic[1] else NA_real_
    
    res_df <- tibble(
      `Cooking Method` = method,
      `Oyster size` = size,
      `Holding Time (min)` = time_min
    )
    res_df[[col_unsafe]] <- n_unsafe
    res_df[[col_safe]]   <- n_safe
    res_df$`Top Associated Parameter` <- top_param
    res_df$`KS D-Statistic` <- top_d
    res_df
  })
}

# 8. U.S. Annual Illness Burden Calculation Functions (FoodNet Exposure Atlas 2006-2007)

#' Calculate U.S. Annual Norovirus Illness Burden from Raw and Cooked Oyster Consumption
calculate_us_illness_burden <- function(
  precook_node,
  cooked_summary_df,
  us_population = 341800000,
  p_raw_weekly = 0.007,
  p_total_weekly = 0.026,
  n_var = 10000,
  n_unc = 500
) {
  p_cooked_weekly <- p_total_weekly - p_raw_weekly
  annual_raw_events_per_person <- p_raw_weekly * 52
  annual_cooked_events_per_person <- p_cooked_weekly * 52

  total_us_raw_events <- annual_raw_events_per_person * us_population
  total_us_cooked_events <- annual_cooked_events_per_person * us_population

  # Evaluate Raw Oyster Risk (Survival S = 1.0)
  node_mat <- unclass(precook_node)[, , 1]
  
  raw_risk_list <- map_dfr(c(6, 12), function(serving_sz) {
    serving_dose_matrix <- matrix(0, nrow = n_var, ncol = n_unc)
    for (u in seq_len(n_unc)) {
      col_doses <- node_mat[, u]
      s_col <- numeric(n_var)
      for (k in seq_len(serving_sz)) {
        s_col <- s_col + sample(col_doses, size = n_var, replace = TRUE)
      }
      serving_dose_matrix[, u] <- s_col
    }
    
    risk_res <- evaluate_teunis_risk_matrix(dose_matrix = serving_dose_matrix)
    
    tibble(
      food_category = "Raw oysters",
      method = "raw",
      size = "all",
      time_min = 0,
      serving_size = serving_sz,
      mean_p_illness = risk_res$mean_p_ill,
      annual_us_eating_events = total_us_raw_events,
      expected_annual_us_illnesses = risk_res$mean_p_ill * total_us_raw_events
    )
  })

  # Evaluate Cooked Oyster Risk from cooked_summary_df
  cooked_risk_list <- cooked_summary_df %>%
    mutate(
      food_category = "Cooked oysters",
      annual_us_eating_events = total_us_cooked_events,
      expected_annual_us_illnesses = mean_p_illness * total_us_cooked_events
    ) %>%
    dplyr::select(
      food_category,
      method,
      size,
      time_min,
      serving_size,
      mean_p_illness,
      annual_us_eating_events,
      expected_annual_us_illnesses
    )

  combined_summary <- bind_rows(raw_risk_list, cooked_risk_list)

  # Create clean manuscript table
  manuscript_table <- combined_summary %>%
    mutate(
      Cooking_Condition = if_else(
        method == "raw",
        "No cooking (Raw)",
        paste(str_to_title(method), str_to_title(size), paste0(time_min, " min"))
      ),
      `Serving Size` = paste(serving_size, "oysters"),
      `Mean P(illness)` = format_risk_label(mean_p_illness),
      `Annual U.S. Eating Events` = format_count_label(annual_us_eating_events),
      `Expected Annual U.S. Illnesses` = format_count_label(expected_annual_us_illnesses)
    ) %>%
    dplyr::select(
      Food = food_category,
      Cooking_Condition,
      `Serving Size`,
      `Mean P(illness)`,
      `Annual U.S. Eating Events`,
      `Expected Annual U.S. Illnesses`
    )

  list(
    full_summary = combined_summary,
    manuscript_table = manuscript_table
  )
}
