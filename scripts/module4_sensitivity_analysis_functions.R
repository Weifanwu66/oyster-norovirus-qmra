# Module 4: Model Diagnostics & Sensitivity Analysis Function Library
# Contains functions for 2D Monte Carlo Convergence Assessment, Spearman Rank Correlation 
# Tornado Plots, and Monte Carlo Filtering (MCF) Regional Sensitivity Analysis.

# 5. Global Spearman Rank Correlation Tornado Plot Function
#' Plot Tornado Sensitivity Diagram for Norovirus Illness Risk P(illness)
#' Plot Combined Tornado Sensitivity Diagram for Norovirus Illness Risk P(illness)
plot_combined_tornado_sensitivity <- function(
  postcook_nodes_list_geeraerd,
  postcook_nodes_list_meta,
  output_filepath = NULL,
  teunis_draws = NULL,
  se_alpha = NULL,
  se_beta = NULL
) {
  set.seed(4040)
  
  if (is.null(teunis_draws)) {
    teunis_draws <- sample_teunis_2020_draws(n_unc = 500)
  }
  
  if (is.null(se_alpha) || is.null(se_beta)) {
    se_fit <- fit_beta_from_mean_ci()
    se_alpha <- se_fit$alpha
    se_beta  <- se_fit$beta
  }
  
  scen_key <- grep("steaming_medium_4min", names(postcook_nodes_list_geeraerd), value = TRUE)[1]
  if (is.na(scen_key)) scen_key <- names(postcook_nodes_list_geeraerd)[1]
  
  # --- Data for Geeraerd ---
  post_mat_g <- unclass(postcook_nodes_list_geeraerd[[scen_key]])[, , 1]
  # --- Data for Meta ---
  post_mat_m <- unclass(postcook_nodes_list_meta[[scen_key]])[, , 1]
  
  n_u_raw <- min(ncol(post_mat_g), ncol(post_mat_m))
  
  csv_path <- here::here("results", "module1_3_exposure_dose_calculation", "module1_uncertainty_inputs.csv")
  if (file.exists(csv_path)) {
    module1_u <- read_csv(csv_path, show_col_types = FALSE)
    n_u <- min(n_u_raw, nrow(module1_u))
    q_inv_val <- module1_u$q_inv
    prev_mean_val <- module1_u$prevalence_mean
    conc_mean_val <- module1_u$positive_mean_gc
    conc_shape_val <- module1_u$concentration_gamma_shape
    conc_scale_val <- conc_mean_val / conc_shape_val  # theta = mu / k
  } else {
    stop("module1_uncertainty_inputs.csv not found!")
  }
  
  teunis_sel <- teunis_draws[sample.int(nrow(teunis_draws), n_u, replace = TRUE), ]
  a_val <- par1.btr(teunis_sel$w_eta, teunis_sel$z_eta)
  b_val <- par2.btr(teunis_sel$w_eta, teunis_sel$z_eta)
  r_val <- par1.btr(teunis_sel$w_mu, teunis_sel$z_mu)
  n_val <- par2.btr(teunis_sel$w_mu, teunis_sel$z_mu)
  p_se_val <- rbeta(n_u, shape1 = se_alpha, shape2 = se_beta)
  
  # Base Uncertainty Input Table
  base_data <- tibble(
    `Harvest~prevalence~(P[prev])` = prev_mean_val,
    `Positive~concentration~scale~(theta)` = conc_scale_val,
    `Positive~concentration~shape~(kappa)` = conc_shape_val,
    `RT-qPCR~infectious~virus~ratio~(Q)` = 1 / q_inv_val,
    `Host~secretor~susceptibility~(P[Se^"+"])` = p_se_val,
    `DR~alpha~(infection)` = a_val,
    `DR~beta~(infection)` = b_val,
    `DR~eta~(illness)` = n_val,
    `DR~rho~(illness)` = r_val
  )
  
  # Process single model
  process_model <- function(post_mat, model_name, therm_df, therm_cols) {
    mean_p_illness <- numeric(n_u)
    for (u in seq_len(n_u)) {
      dose_u <- post_mat[, u]
      p_inf_se <- p.inf.hyp1f1(a_val[u], b_val[u], dose_u)
      p_ill_given_inf <- p.ill.inf(n_val[u], r_val[u], dose_u)
      p_ill_u <- p_se_val[u] * p_inf_se * p_ill_given_inf
      mean_p_illness[u] <- mean(p_ill_u)
    }
    
    input_data <- bind_cols(base_data, therm_cols)
    
    map_dfr(names(input_data), function(var_name) {
      rho <- cor(input_data[[var_name]], mean_p_illness, method = "spearman", use = "pairwise.complete.obs")
      tibble(
        Model = model_name,
        input = var_name,
        spearman_rho = rho,
        abs_rho = abs(rho),
        direction = if_else(rho >= 0, "Positive Association", "Negative Association")
      )
    }) %>%
      filter(!is.na(spearman_rho))
  }
  
  # Geeraerd
  therm_post_g <- readRDS(here::here("results", "module1_2_thermal_inactivation", "module1_2_geeraerd_thermal_kinetics_posterior.rds"))
  set.seed(1503)
  therm_g <- therm_post_g %>% dplyr::slice_sample(n = n_u, replace = TRUE)
  cols_g <- tibble(
    `Geeraerd~thermal~sensitivity~(z)` = therm_g$z,
    `Geeraerd~reference~inactivation~time~(D[r])` = therm_g$D_R,
    `Geeraerd~residual~tail~(log[10]*N[min])` = therm_g$logNmin,
    `Geeraerd~initial~shoulder~(log[10]*N[0])` = therm_g$logN0
  )
  df_g <- process_model(post_mat_g, "Geeraerd model", therm_g, cols_g)
  
  # Mafart
  therm_post_m <- readRDS(here::here("results", "module1_2_thermal_inactivation", "module1_2_meta_analysis_thermal_kinetics_posterior.rds"))
  set.seed(1503)
  therm_m <- therm_post_m %>% dplyr::slice_sample(n = n_u, replace = TRUE)
  cols_m <- tibble(
    `Mafart-Weibull~thermal~sensitivity~(z)` = therm_m$z_shellfish,
    `Mafart-Weibull~reference~inactivation~time~(D[r])` = therm_m$Dref_shell,
    `Mafart-Weibull~shape~parameter~(p)` = therm_m$p_weibull
  )
  df_m <- process_model(post_mat_m, "Mafart-Weibull model", therm_m, cols_m)
  
  cor_df <- bind_rows(df_g, df_m)
  
  # Ordered levels per facet can be tricky, so we order by abs_rho generally
  # by making input a factor. But each facet has different inputs. We'll order by abs_rho
  # using tidytext::reorder_within
  cor_df <- cor_df %>%
    arrange(abs_rho) %>% mutate(input_ordered = factor(paste(input, Model, sep = "___"), levels = paste(input, Model, sep = "___")))
  
  p <- ggplot(cor_df, aes(x = spearman_rho, y = input_ordered, fill = direction)) +
    geom_col(width = 0.65) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~ Model, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Positive Association" = "#2F6F9F", "Negative Association" = "#B65C45")) +
    scale_y_discrete(labels = function(x) parse(text = sub("___.*$", "", x))) +
    scale_x_continuous(limits = c(-0.45, 0.45), breaks = seq(-0.4, 0.4, by = 0.1)) +
    labs(
      x = expression(paste("Spearman rank correlation coefficient", " (", rho, ")")),
      y = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "#f0f0f0", color = "NA"),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    )
  
  if (!is.null(output_filepath)) ggsave(filename = output_filepath, plot = p, width = 12, height = 5.5, dpi = 300, compression = "lzw")
  
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
  n_u_raw <- ncol(post_mat)
  
  csv_path <- here::here("results", "module1_3_exposure_dose_calculation", "module1_uncertainty_inputs.csv")
  if (file.exists(csv_path)) {
    module1_u <- read_csv(csv_path, show_col_types = FALSE)
    n_u <- min(n_u_raw, nrow(module1_u))
  } else {
    stop("module1_uncertainty_inputs.csv not found!")
  }
  
  idx_v <- sample.int(n_v, n_sample, replace = TRUE)
  idx_u <- sample.int(n_u, n_sample, replace = TRUE)
  
  postcook_val <- sapply(seq_len(n_sample), function(i) post_mat[idx_v[i], idx_u[i]])
  p_se_val <- rbeta(n_sample, shape1 = se_alpha, shape2 = se_beta)
  
    q_inv_val <- module1_u$q_inv[idx_u]
    prev_mean_val <- module1_u$prevalence_mean[idx_u]
    conc_mean_val <- module1_u$positive_mean_gc[idx_u]
    conc_shape_val <- module1_u$concentration_gamma_shape[idx_u]
  
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
    
    ks_res <- suppressWarnings(ks.test(val_unsafe, val_safe))
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
    ks_table = ks_table,
    cutoff = cutoff
  )
}

#' Plot Monte Carlo Filtering Cumulative Distribution Function (CDF) Comparison Curves
#' Plot Monte Carlo Filtering Cumulative Distribution Function (CDF) Comparison Curves
plot_mcf_cdf_curves <- function(mcf_result, output_filepath = NULL) {
  top_param <- mcf_result$ks_table$Parameter[1]
  
  unsafe_label <- sprintf("Unsafe / High-Risk (P >= %g)", mcf_result$cutoff)
  safe_label   <- sprintf("Safe / Compliant (P < %g)", mcf_result$cutoff)

  df_plot <- mcf_result$input_data %>%
    mutate(
      is_unsafe = mcf_result$is_unsafe,
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
  
  if (!is.null(output_filepath)) ggsave(filename = output_filepath, plot = p, width = 6.5, height = 4.5, dpi = 300, compression = "lzw")
  
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
