# Messner et al. fractional-Poisson parameter P
# Reported estimate: P = 0.72, 95% CrI = 0.63–0.80
# Choose alpha and beta so that:
#   mean(P) = alpha / (alpha + beta) = 0.72
#   qbeta(0.025) = 0.63
#   qbeta(0.975) = 0.80

messner_p_mean <- 0.72
messner_p_q025 <- 0.63
messner_p_q975 <- 0.80

# Other multiplicative parameters
p_host <- 0.74
p_ill_given_inf <- 0.30

sample_serving_dose <- function(post_dose_vector, serving_size) {
  sampled_oysters <- matrix(
    sample(
      post_dose_vector,
      size = length(post_dose_vector) * serving_size,
      replace = TRUE
    ),
    ncol = serving_size
  )

  rowSums(sampled_oysters)
}


make_serving_dose_from_one_condition <- function(
    post_cook_table,
    serving_sizes = c(6, 12)
) {

  post_cook_table %>%
    dplyr::group_by(uncertainty_iteration) %>%
    dplyr::group_modify(~{
      purrr::map_dfr(serving_sizes, function(serving_size) {

        dose_vec <- sample_serving_dose(
          post_dose_vector = .x$post_cook_IU_per_oyster,
          serving_size = serving_size
        )

        tibble::tibble(
          target_temperature_C = unique(.x$target_temperature_C),
          time_min = unique(.x$time_min),
          serving_size = serving_size,
          serving_variability_iteration = seq_along(dose_vec),
          dose_IU_per_serving = dose_vec
        )
      })
    }) %>%
    ungroup()
}

fit_beta_from_mean_ci <- function(mean_value, q025, q975) {

  objective <- function(log_par) {

    alpha <- exp(log_par[1])
    beta  <- exp(log_par[2])

    fitted_mean <- alpha / (alpha + beta)
    fitted_q025 <- qbeta(0.025, shape1 = alpha, shape2 = beta)
    fitted_q975 <- qbeta(0.975, shape1 = alpha, shape2 = beta)

    # Minimize mismatch to all three targets:
    # mean, lower 95% CrI, upper 95% CrI
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

messner_p_beta_fit <- fit_beta_from_mean_ci(
  mean_value = messner_p_mean,
  q025 = messner_p_q025,
  q975 = messner_p_q975
)

messner_p_alpha <- messner_p_beta_fit$alpha
messner_p_beta <- messner_p_beta_fit$beta

# Check fit
messner_p_beta_fit

messner_fractional_poisson_risk <- function(
    serving_dose_table,
    dose_response_uncertainty = NULL
) {

  if (is.null(dose_response_uncertainty)) {
    dose_response_uncertainty <- serving_dose_table %>%
      distinct(uncertainty_iteration) %>%
      mutate(
        p_messner = rbeta(
          n = n(),
          shape1 = messner_p_alpha,
          shape2 = messner_p_beta
        )
      )
  } else {
    dose_response_uncertainty <- dose_response_uncertainty %>%
      dplyr::select(uncertainty_iteration, p_messner)
  }

  serving_dose_table %>%
    left_join(dose_response_uncertainty, by = "uncertainty_iteration") %>%
    mutate(
      p_messner = dplyr::coalesce(
        p_messner,
        rbeta(
          n = n(),
          shape1 = messner_p_alpha,
          shape2 = messner_p_beta
        )
      ),
      messner_p_mean = messner_p_mean,
      messner_p_q025 = messner_p_q025,
      messner_p_q975 = messner_p_q975,
      messner_p_alpha = messner_p_alpha,
      messner_p_beta = messner_p_beta,
      p_host = p_host,
      p_ill_given_inf = p_ill_given_inf,
      p_exposure = 1 - exp(-dose_IU_per_serving),
      p_infection = p_host * p_messner * p_exposure,
      p_illness = p_ill_given_inf * p_infection
    )
}

summarize_module4_risk_by_uncertainty <- function(risk_table) {

  risk_table %>%
    ungroup() %>%
    mutate(
      target_temperature_C = parse_number(as.character(target_temperature_C)),
      time_min = parse_number(as.character(time_min)),
      serving_size = parse_number(as.character(serving_size))
    ) %>%
    group_by(target_temperature_C, time_min, serving_size, uncertainty_iteration) %>%
    summarise(
      mean_dose_IU_per_serving = mean(dose_IU_per_serving, na.rm = TRUE),
      mean_p_infection = mean(p_infection, na.rm = TRUE),
      mean_p_illness = mean(p_illness, na.rm = TRUE),
      p05_dose_IU_per_serving = quantile(dose_IU_per_serving, 0.05, na.rm = TRUE),
      p50_dose_IU_per_serving = quantile(dose_IU_per_serving, 0.50, na.rm = TRUE),
      p95_dose_IU_per_serving = quantile(dose_IU_per_serving, 0.95, na.rm = TRUE),
      p05_p_infection = quantile(p_infection, 0.05, na.rm = TRUE),
      p50_p_infection = quantile(p_infection, 0.50, na.rm = TRUE),
      p95_p_infection = quantile(p_infection, 0.95, na.rm = TRUE),
      p05_p_illness = quantile(p_illness, 0.05, na.rm = TRUE),
      p50_p_illness = quantile(p_illness, 0.50, na.rm = TRUE),
      p95_p_illness = quantile(p_illness, 0.95, na.rm = TRUE),
      n_variability = n_distinct(serving_variability_iteration),
      .groups = "drop"
    )
}

summarize_mc2d_outcome <- function(condition_table, outcome_col) {

  condition_table <- condition_table %>%
    arrange(uncertainty_iteration, serving_variability_iteration)

  n_variability <- n_distinct(condition_table$serving_variability_iteration)
  n_uncertainty <- n_distinct(condition_table$uncertainty_iteration)

  risk_summary <- condition_table[[outcome_col]] %>%
    mcdata(type = "VU", nsv = n_variability, nsu = n_uncertainty) %>%
    summary(probs = c(0.05, 0.50, 0.95), lim = c(0.025, 0.975)) %>%
    .[[1]]

  as_tibble(
    risk_summary[, c("mean", "5%", "50%", "95%"), drop = FALSE],
    rownames = "uncertainty_statistic"
  ) %>%
    pivot_longer(
      cols = -uncertainty_statistic,
      names_to = "variability_statistic",
      values_to = "value"
    ) %>%
    mutate(
      outcome = outcome_col,
      n_uncertainty = n_uncertainty,
      n_variability = n_variability,
      variability_statistic = recode(
        variability_statistic,
        "5%" = "p05",
        "50%" = "p50",
        "95%" = "p95"
      ),
      uncertainty_statistic = recode(
        uncertainty_statistic,
        "2.5%" = "p025",
        "97.5%" = "p975"
      )
    ) %>%
    select(
      outcome, variability_statistic, uncertainty_statistic,
      value, n_uncertainty, n_variability
    )
}

summarize_module4_risk_mc2d <- function(risk_table) {

  risk_table %>%
    ungroup() %>%
    mutate(
      target_temperature_C = parse_number(as.character(target_temperature_C)),
      time_min = parse_number(as.character(time_min)),
      serving_size = parse_number(as.character(serving_size))
    ) %>%
    group_by(target_temperature_C, time_min, serving_size) %>%
    group_modify(~ bind_rows(
      summarize_mc2d_outcome(.x, "dose_IU_per_serving"),
      summarize_mc2d_outcome(.x, "p_infection"),
      summarize_mc2d_outcome(.x, "p_illness")
    )) %>%
    ungroup() %>%
    pivot_wider(
      names_from = uncertainty_statistic,
      values_from = value,
      names_prefix = "uncertainty_"
    ) %>%
    arrange(target_temperature_C, time_min, serving_size, outcome, variability_statistic)
}

summarize_condition_risk_compact <- function(
    pre_dose_table,
    condition_survival_table,
    dose_response_uncertainty,
    serving_sizes = c(6, 12)
) {
  base_gc <- pre_dose_table$gii_total_genome_copies_per_oyster
  n_variability <- length(base_gc)

  condition_inputs <- condition_survival_table %>%
    dplyr::left_join(
      dose_response_uncertainty %>%
        dplyr::select(uncertainty_iteration, p_messner),
        by = "uncertainty_iteration"
    )

  purrr::map_dfr(serving_sizes, function(serving_size) {
    serving_gc <- sample_serving_dose(base_gc, serving_size)
    mean_serving_gc <- mean(serving_gc, na.rm = TRUE)
    q05_serving_gc <- quantile(serving_gc, 0.05, na.rm = TRUE)
    q50_serving_gc <- quantile(serving_gc, 0.50, na.rm = TRUE)
    q95_serving_gc <- quantile(serving_gc, 0.95, na.rm = TRUE)

    dose_scale <- condition_inputs$infectious_units_per_genome_copy *
      condition_inputs$survival_fraction
    mean_p_exposure <- colMeans(1 - exp(-outer(serving_gc, dose_scale)))
    mean_dose <- mean_serving_gc * dose_scale
    mean_p_infection <- p_host * condition_inputs$p_messner * mean_p_exposure
    p05_dose <- as.numeric(q05_serving_gc) * dose_scale
    p50_dose <- as.numeric(q50_serving_gc) * dose_scale
    p95_dose <- as.numeric(q95_serving_gc) * dose_scale
    p05_p_infection <- p_host * condition_inputs$p_messner * (1 - exp(-p05_dose))
    p50_p_infection <- p_host * condition_inputs$p_messner * (1 - exp(-p50_dose))
    p95_p_infection <- p_host * condition_inputs$p_messner * (1 - exp(-p95_dose))

    tibble::tibble(
        target_temperature_C = condition_inputs$target_temperature_C,
        time_min = condition_inputs$time_min,
        serving_size = serving_size,
        uncertainty_iteration = condition_inputs$uncertainty_iteration,
        mean_dose_IU_per_serving = mean_dose,
        mean_p_infection = mean_p_infection,
        mean_p_illness = p_ill_given_inf * mean_p_infection,
        p05_dose_IU_per_serving = p05_dose,
        p50_dose_IU_per_serving = p50_dose,
        p95_dose_IU_per_serving = p95_dose,
        p05_p_infection = p05_p_infection,
        p50_p_infection = p50_p_infection,
        p95_p_infection = p95_p_infection,
        p05_p_illness = p_ill_given_inf * p05_p_infection,
        p50_p_illness = p_ill_given_inf * p50_p_infection,
        p95_p_illness = p_ill_given_inf * p95_p_infection,
        n_variability = n_variability
      )
  })
}
