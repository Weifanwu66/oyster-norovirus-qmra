# Module 3: Risk Characterization Function Library
# Contains Serving-Level Risk Evaluation with Fully Random Multinomial Size Mixing, Heatmap Visualization,
# and U.S. Population Illness Burden Calculations using FoodNet Data.

library(gsl)
library(MASS)


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


plot_module4_risk_heatmap <- function(summary_df, output_filepath = NULL, risk_metric = "mean_p_illness", threshold = 1e-4, include_half_dozen = FALSE) {
  if (!include_half_dozen) {
    summary_df <- summary_df %>% filter(serving_size == 12)
  }
  
  heatmap_data <- summary_df %>%
    mutate(
      method_label = str_to_title(method),
      size_label = factor(str_to_title(size), levels = c("Small", "Medium", "Large", "Mixed")),
      serving_label = paste0(serving_size, " Oysters / Serving"),
      plot_risk = pmax(.data[[risk_metric]], 1e-6),
      below_threshold = plot_risk < threshold
    )

  metric_title <- if (risk_metric == "mean_p_illness") "Mean P(ill)" else "Median P(ill)"

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
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      name = metric_title,
      guide = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(8, "cm"), barheight = unit(0.4, "cm"))
    ) +
    labs(
      x = "Cooking holding time (minutes)",
      y = "Oyster size category"
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
