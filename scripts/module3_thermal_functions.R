# Module 3 thermal workflows
# Observed fitting/overlay uses measured replicate-specific profiles.
# Scenario prediction uses averaged/interpolated representative profiles.

prepare_observed_dynamic_profiles <- function(d) {
  d %>%
    dplyr::mutate(
      cond_nominal = sub("-rep[^-]+$", "", .data$cond),
      cond = dplyr::if_else(
        grepl("-rep[^-]+$", .data$cond),
        .data$cond,
        paste0(.data$cond, "-rep", .data$rep)
      )
    ) %>%
    dplyr::select(cond_nominal, cond, rep, time, temperature, logN, dplyr::everything())
}

module3_geeraerd_to_bioinactivation_pars <- function(parameters) {
  all_pars <- as.list(parameters)
  all_pars$N0 <- 10^all_pars$logN0
  all_pars$N_min <- 10^all_pars$logNmin
  all_pars$logN0 <- NULL
  all_pars$logNmin <- NULL
  all_pars
}

get_geeraerd_residuals <- function(p, this_data, known, times, this_temp) {
  pred <- bioinactivation::predict_inactivation(
    "Geeraerd",
    times,
    module3_geeraerd_to_bioinactivation_pars(c(p, known)),
    this_temp
  )$simulation %>%
    as.data.frame()

  merged <- merge(
    this_data[, c("time", "logN")],
    pred[, c("time", "logN")],
    by = "time",
    all.x = TRUE,
    suffixes = c("_obs", "_mod")
  )

  merged$logN_mod - merged$logN_obs
}

get_all_geeraerd_residuals <- function(p, split_data, known, times_list) {
  res <- mapply(
    function(dat, times) {
      this_temp <- dat %>%
        dplyr::select(time, temperature) %>%
        dplyr::distinct(time, .keep_all = TRUE) %>%
        dplyr::arrange(time)

      get_geeraerd_residuals(p, dat, known, times, this_temp)
    },
    split_data,
    times_list,
    SIMPLIFY = FALSE
  )

  sum(unlist(res)^2)
}

fit_global_dynamic_geeraerd <- function(d, guess, known, lower = -Inf, upper = Inf, niter = 5000) {
  d <- d %>%
    dplyr::mutate(time = dplyr::if_else(.data$time < 2e-5, 2e-5, .data$time))

  split_data <- split(d, d$cond)
  times_list <- lapply(split_data, function(x) sort(unique(x$time)))

  FME::modMCMC(
    get_all_geeraerd_residuals,
    p = guess,
    split_data = split_data,
    known = known,
    times_list = times_list,
    lower = lower,
    upper = upper,
    niter = niter
  )
}

make_geeraerd_observed_fit_overlay <- function(
    posterior_draws,
    observed_data,
    known,
    n_posterior_draws = 1000,
    n_time_points = 150,
    seed = 1503
) {
  set.seed(seed)

  observed_data <- observed_data %>%
    dplyr::mutate(
      time = dplyr::if_else(.data$time < 2e-5, 2e-5, .data$time),
      cond = factor(.data$cond, levels = sort(unique(.data$cond))),
      target_temperature_C = as.numeric(sub(".*-([0-9]+(?:\\.[0-9]+)?)-.*", "\\1", .data$cond_nominal))
    ) %>%
    dplyr::arrange(.data$target_temperature_C) %>%
    dplyr::mutate(
      cond_label = factor(
        paste0(.data$target_temperature_C, " °C"),
        levels = paste0(sort(unique(.data$target_temperature_C)), " °C")
      )
    )

  posterior_draws <- posterior_draws %>%
    dplyr::slice_sample(n = min(n_posterior_draws, nrow(posterior_draws)))

  condition_grid <- observed_data %>%
    dplyr::group_by(.data$cond) %>%
    dplyr::summarise(
      cond_nominal = dplyr::first(.data$cond_nominal),
      cond_label = dplyr::first(.data$cond_label),
      min_time = min(.data$time),
      max_time = max(.data$time),
      .groups = "drop"
    )

  prediction_draws <- purrr::map_dfr(seq_len(nrow(condition_grid)), function(i) {
    this_data <- observed_data %>%
      dplyr::filter(.data$cond == condition_grid$cond[i])

    this_temp <- this_data %>%
      dplyr::select(time, temperature) %>%
      dplyr::distinct(time, .keep_all = TRUE) %>%
      dplyr::arrange(time)

    this_times <- seq(condition_grid$min_time[i], condition_grid$max_time[i], length.out = n_time_points)

    purrr::map_dfr(seq_len(nrow(posterior_draws)), function(j) {
      this_pars <- unlist(posterior_draws[j, c("D_R", "z", "logN0", "logNmin")])

      pred <- bioinactivation::predict_inactivation(
        "Geeraerd",
        this_times,
        module3_geeraerd_to_bioinactivation_pars(c(this_pars, known)),
        this_temp
      )$simulation %>%
        as.data.frame()

      tibble::tibble(
        cond_nominal = condition_grid$cond_nominal[i],
        cond_label = condition_grid$cond_label[i],
        cond = condition_grid$cond[i],
        posterior_draw = j,
        time = pred$time,
        fitted_logN = pred$logN
      )
    })
  })

  prediction_summary <- prediction_draws %>%
    dplyr::group_by(.data$cond_nominal, .data$cond_label, .data$time) %>%
    dplyr::summarise(
      fitted_median = median(.data$fitted_logN, na.rm = TRUE),
      fitted_lower_95 = quantile(.data$fitted_logN, 0.025, na.rm = TRUE),
      fitted_upper_95 = quantile(.data$fitted_logN, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  overlay_plot <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = prediction_summary,
      ggplot2::aes(x = time, ymin = fitted_lower_95, ymax = fitted_upper_95),
      fill = "#9ECAE1",
      alpha = 0.45
    ) +
    ggplot2::geom_line(
      data = prediction_summary,
      ggplot2::aes(x = time, y = fitted_median),
      color = "#08519C",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = observed_data,
      ggplot2::aes(x = time, y = logN),
      color = "black",
      fill = "#F03B20",
      shape = 21,
      size = 2.1,
      stroke = 0.6
    ) +
    ggplot2::facet_wrap(~cond_label, ncol = 2) +
    ggplot2::scale_x_continuous(limits = c(0, 20), breaks = 0:20) +
    ggplot2::labs(
      x = "Thermal-process duration (min)",
      y = expression(log[10] * " PFU")
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )

  list(prediction_summary = prediction_summary, plot = overlay_plot)
}

make_interpolated_temp_profile <- function(df_comeup, target_C, times_min) {
  # Scenario prediction only: averages replicate come-up profiles.
  profile_grid <- df_comeup %>%
    dplyr::group_by(target_temp, time) %>%
    dplyr::summarise(sample_temp = mean(sample_temp, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(target_temp) %>%
    dplyr::group_modify(~{
      tibble::tibble(
        time = times_min,
        sample_temp = approx(.x$time, .x$sample_temp, xout = times_min, rule = 2)$y
      )
    }) %>%
    dplyr::ungroup()

  profile_grid %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(
      temperature = approx(target_temp, sample_temp, xout = target_C, rule = 2)$y,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      cond = paste0("target_", target_C, "C_interpolated_representative"),
      profile_type = "interpolated_representative",
      target_temperature_C = target_C
    ) %>%
    dplyr::select(cond, profile_type, target_temperature_C, time, temperature)
}

predict_geeraerd_dynamic_reduction <- function(
    theta,
    temp_profile,
    times_min,
    known = list(temp_ref = 57, C_c0 = 1e-3)
) {
  theta <- as.list(theta)

  pars <- c(
    known,
    list(
      D_R = theta$D_R,
      z = theta$z,
      N0 = 10^theta$logN0,
      N_min = 10^theta$logNmin
    )
  )

  time_pred <- ifelse(times_min == 0, 2e-5, times_min)
  temp_profile <- temp_profile %>%
    dplyr::mutate(time = ifelse(time == 0, 2e-5, time))

  sim <- bioinactivation::predict_inactivation(
    "Geeraerd",
    time_pred,
    pars,
    temp_profile
  )$simulation

  tibble::tibble(
    time_min = times_min,
    pred_logN = sim$logN,
    log_reduction = theta$logN0 - sim$logN,
    survival_fraction = 10^(-(theta$logN0 - sim$logN))
  )
}

make_module3_time_grid <- function(target_time_grid) {
  target_time_grid %>%
    dplyr::transmute(
      target_temperature_C = as.numeric(.data$target_temperature_C),
      time_min = as.numeric(.data$time_min)
    ) %>%
    dplyr::filter(!is.na(.data$target_temperature_C), !is.na(.data$time_min)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data$target_temperature_C, .data$time_min)
}

run_module3_thermal_inactivation <- function(
    geeraerd_post,
    df_comeup,
    n_uncertainty,
    target_time_grid,
    known = list(temp_ref = 57, C_c0 = 1e-3)
) {
  time_grid <- make_module3_time_grid(target_time_grid)

  theta_ge <- geeraerd_post %>%
    dplyr::slice_sample(n = n_uncertainty, replace = TRUE) %>%
    dplyr::mutate(uncertainty_iteration = dplyr::row_number())

  thermal_profiles <- purrr::map_dfr(unique(time_grid$target_temperature_C), function(target_C) {
    this_times <- time_grid %>%
      dplyr::filter(.data$target_temperature_C == target_C) %>%
      dplyr::pull(.data$time_min)

    make_interpolated_temp_profile(df_comeup, target_C, this_times)
  })

  thermal_survival <- tidyr::expand_grid(
    uncertainty_iteration = theta_ge$uncertainty_iteration,
    target_temperature_C = unique(time_grid$target_temperature_C)
  ) %>%
    dplyr::left_join(theta_ge, by = "uncertainty_iteration") %>%
    dplyr::group_by(uncertainty_iteration, target_temperature_C) %>%
    dplyr::group_modify(~{
      this_times <- time_grid %>%
        dplyr::filter(target_temperature_C == .y$target_temperature_C) %>%
        dplyr::pull(.data$time_min)

      temp_profile <- thermal_profiles %>%
        dplyr::filter(target_temperature_C == .y$target_temperature_C) %>%
        dplyr::select(time, temperature)

      predict_geeraerd_dynamic_reduction(.x[1, ], temp_profile, this_times, known)
    }) %>%
    dplyr::ungroup()

  thermal_summary <- thermal_survival %>%
    dplyr::group_by(target_temperature_C, time_min) %>%
    dplyr::summarise(
      median_log_reduction = median(log_reduction, na.rm = TRUE),
      lower_95_log_reduction = quantile(log_reduction, 0.025, na.rm = TRUE),
      upper_95_log_reduction = quantile(log_reduction, 0.975, na.rm = TRUE),
      median_survival_fraction = median(survival_fraction, na.rm = TRUE),
      lower_95_survival_fraction = quantile(survival_fraction, 0.025, na.rm = TRUE),
      upper_95_survival_fraction = quantile(survival_fraction, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    thermal_survival = thermal_survival,
    thermal_summary = thermal_summary,
    thermal_profiles = thermal_profiles,
    time_grid = time_grid,
    theta_ge = theta_ge
  )
}

make_post_cook_dose_for_condition <- function(
    target_C,
    time_value,
    pre_dose_table,
    survival_fraction_table,
    module3_results_dir = NULL,
    save_full = TRUE
) {
  sf_this <- survival_fraction_table %>%
    dplyr::filter(target_temperature_C == target_C, time_min == time_value) %>%
    dplyr::select(
      uncertainty_iteration,
      target_temperature_C,
      time_min,
      log_reduction,
      survival_fraction,
      infectious_units_per_genome_copy
    )

  if (!save_full) {
    post_values <- unlist(lapply(seq_len(nrow(sf_this)), function(i) {
      pre_dose_table$gii_total_genome_copies_per_oyster *
        sf_this$infectious_units_per_genome_copy[i] *
        sf_this$survival_fraction[i]
    }))

    return(tibble::tibble(
      target_temperature_C = target_C,
      time_min = time_value,
      n_variability = dplyr::n_distinct(pre_dose_table$variability_iteration),
      n_uncertainty = dplyr::n_distinct(sf_this$uncertainty_iteration),
      mean_post_cook_IU_per_oyster = mean(post_values, na.rm = TRUE),
      median_post_cook_IU_per_oyster = median(post_values, na.rm = TRUE),
      p90_post_cook_IU_per_oyster = quantile(post_values, 0.90, na.rm = TRUE),
      p95_post_cook_IU_per_oyster = quantile(post_values, 0.95, na.rm = TRUE),
      p99_post_cook_IU_per_oyster = quantile(post_values, 0.99, na.rm = TRUE),
      max_post_cook_IU_per_oyster = max(post_values, na.rm = TRUE),
      median_log_reduction = median(sf_this$log_reduction, na.rm = TRUE),
      median_survival_fraction = median(sf_this$survival_fraction, na.rm = TRUE)
    ))
  }

  if ("uncertainty_iteration" %in% names(pre_dose_table)) {
    post_this <- pre_dose_table %>%
      dplyr::inner_join(sf_this, by = "uncertainty_iteration")
  } else {
    post_this <- tidyr::crossing(pre_dose_table, sf_this)
  }

  post_this <- post_this %>%
    dplyr::mutate(
      pre_cook_IU_per_oyster =
        gii_total_genome_copies_per_oyster *
        infectious_units_per_genome_copy,
      post_cook_IU_per_oyster = pre_cook_IU_per_oyster * survival_fraction
    )

  if (save_full) {
    if (is.null(module3_results_dir)) {
      stop("module3_results_dir is required when save_full = TRUE.", call. = FALSE)
    }
    dir.create(module3_results_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(post_this, file.path(module3_results_dir, paste0("module3_post_cook_IU_target_", gsub("\\.", "p", target_C), "C_time_", gsub("\\.", "p", time_value),"min.rds")
      )
    )
  }

  post_this %>%
    dplyr::summarise(
      target_temperature_C = target_C,
      time_min = time_value,
      n_variability = dplyr::n_distinct(variability_iteration),
      n_uncertainty = dplyr::n_distinct(uncertainty_iteration),
      mean_post_cook_IU_per_oyster = mean(post_cook_IU_per_oyster, na.rm = TRUE),
      median_post_cook_IU_per_oyster = median(post_cook_IU_per_oyster, na.rm = TRUE),
      p90_post_cook_IU_per_oyster = quantile(post_cook_IU_per_oyster, 0.90, na.rm = TRUE),
      p95_post_cook_IU_per_oyster = quantile(post_cook_IU_per_oyster, 0.95, na.rm = TRUE),
      p99_post_cook_IU_per_oyster = quantile(post_cook_IU_per_oyster, 0.99, na.rm = TRUE),
      max_post_cook_IU_per_oyster = max(post_cook_IU_per_oyster, na.rm = TRUE),
      median_log_reduction = median(log_reduction, na.rm = TRUE),
      median_survival_fraction = median(survival_fraction, na.rm = TRUE),
      .groups = "drop"
    )
}
