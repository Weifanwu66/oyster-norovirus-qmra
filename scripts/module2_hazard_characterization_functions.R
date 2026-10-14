# Module 2: Hazard Characterization Function Library
# Contains Host Susceptibility Beta fitting (Messner et al. 2014) and
# Teunis et al. (2020) Dose-Response parameter sampling and mathematics.

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
sample_teunis_2020_draws <- function(n_unc = 1000, seed = 1060L) {
  withr::local_seed(seed,
    .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection")
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


# Same U-level host draws are reused by risk and sensitivity calculations.
sample_host_susceptibility <- function(n_unc, se_alpha, se_beta) {
  withr::with_seed(2020L,
    rbeta(n_unc, shape1 = se_alpha, shape2 = se_beta),
    .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection")
}
