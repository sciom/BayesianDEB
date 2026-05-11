# ===========================================================
# Tests: priors.R  (prior_*, prior_default, prior_to_stan_data*)
# ===========================================================

# --- Constructors: happy path ---

test_that("prior constructors return correct class", {
  p <- prior_lognormal(mu = 1, sigma = 0.5)
  expect_s3_class(p, "bdeb_prior")
  expect_equal(p$family, "lognormal")
  expect_equal(p$mu, 1)
  expect_equal(p$sigma, 0.5)

  p2 <- prior_beta(a = 3, b = 2)
  expect_s3_class(p2, "bdeb_prior")
  expect_equal(p2$family, "beta")
  expect_equal(p2$a, 3)
  expect_equal(p2$b, 2)

  p3 <- prior_halfnormal(sigma = 0.1)
  expect_s3_class(p3, "bdeb_prior")
  expect_equal(p3$family, "halfnormal")

  p4 <- prior_exponential(rate = 2)
  expect_s3_class(p4, "bdeb_prior")
  expect_equal(p4$rate, 2)
})

test_that("prior_default returns complete sets", {
  for (tp in c("individual", "growth_repro", "hierarchical", "debtox")) {
    d <- prior_default(tp)
    expect_type(d, "list")
    # All entries should be bdeb_prior
    for (nm in names(d)) {
      expect_s3_class(d[[nm]], "bdeb_prior")
    }
  }

  # Individual should have base params
  d_ind <- prior_default("individual")
  expect_true(all(c("p_Am", "p_M", "kappa", "v", "E_G", "E0", "L0", "sigma_L")
                  %in% names(d_ind)))

  # Growth_repro adds k_J, E_Hp, k_R, phi_R
  d_gr <- prior_default("growth_repro")
  expect_true(all(c("k_J", "E_Hp", "k_R", "phi_R") %in% names(d_gr)))

  # Hierarchical adds mu_log_p_Am, sigma_log_p_Am
  d_h <- prior_default("hierarchical")
  expect_true(all(c("mu_log_p_Am", "sigma_log_p_Am") %in% names(d_h)))

  # DEBtox adds k_d, z_w, b_w
  d_t <- prior_default("debtox")
  expect_true(all(c("k_d", "z_w", "b_w") %in% names(d_t)))
})


# --- Constructors: edge cases and defaults ---

test_that("prior_lognormal uses correct defaults", {
  p <- prior_lognormal()
  expect_equal(p$mu, 0)
  expect_equal(p$sigma, 1)
})

test_that("prior_normal uses correct defaults", {
  p <- prior_normal()
  expect_equal(p$mu, 0)
  expect_equal(p$sigma, 1)
  expect_equal(p$family, "normal")
})

test_that("prior_beta uses correct defaults", {
  p <- prior_beta()
  expect_equal(p$a, 2)
  expect_equal(p$b, 2)
})

test_that("prior_halfnormal uses correct defaults", {
  p <- prior_halfnormal()
  expect_equal(p$sigma, 1)
})

test_that("prior_halfcauchy returns correct structure", {
  p <- prior_halfcauchy(sigma = 2.5)
  expect_s3_class(p, "bdeb_prior")
  expect_equal(p$family, "halfcauchy")
  expect_equal(p$sigma, 2.5)
})

test_that("prior_exponential uses correct defaults", {
  p <- prior_exponential()
  expect_equal(p$rate, 1)
})

test_that("prior constructors accept extreme values", {
  p1 <- prior_lognormal(mu = -100, sigma = 0.001)
  expect_equal(p1$mu, -100)
  expect_equal(p1$sigma, 0.001)

  p2 <- prior_lognormal(mu = 50, sigma = 100)
  expect_equal(p2$mu, 50)

  p3 <- prior_beta(a = 0.01, b = 0.01)
  expect_equal(p3$a, 0.01)

  p4 <- prior_beta(a = 100, b = 100)
  expect_equal(p4$a, 100)
})


# --- prior_default: fail scenarios ---

test_that("prior_default rejects invalid type", {
  expect_error(prior_default("nonexistent"), "should be one of")
})

test_that("prior_default returns different sets for each type", {
  d_ind <- prior_default("individual")
  d_gr  <- prior_default("growth_repro")
  d_h   <- prior_default("hierarchical")
  d_t   <- prior_default("debtox")

  # growth_repro has more parameters
  expect_true(length(d_gr) > length(d_ind))
  # hierarchical has mu_log_p_Am but not k_J
  expect_true("mu_log_p_Am" %in% names(d_h))
  expect_false("k_J" %in% names(d_h))
  # debtox has k_d but not mu_log_p_Am
  expect_true("k_d" %in% names(d_t))
  expect_false("mu_log_p_Am" %in% names(d_t))
})


# --- prior_to_stan_data: correctness ---

test_that("prior_to_stan_data maps all individual priors", {
  priors <- prior_default("individual")
  sd <- prior_to_stan_data(priors)

  expect_equal(sd$prior_p_Am_mu, priors$p_Am$mu)
  expect_equal(sd$prior_p_Am_sd, priors$p_Am$sigma)
  expect_equal(sd$prior_p_M_mu, priors$p_M$mu)
  expect_equal(sd$prior_p_M_sd, priors$p_M$sigma)
  expect_equal(sd$prior_kappa_a, priors$kappa$a)
  expect_equal(sd$prior_kappa_b, priors$kappa$b)
  expect_equal(sd$prior_v_mu, priors$v$mu)
  expect_equal(sd$prior_v_sd, priors$v$sigma)
  expect_equal(sd$prior_E_G_mu, priors$E_G$mu)
  expect_equal(sd$prior_E_G_sd, priors$E_G$sigma)
  expect_equal(sd$prior_E0_mu, priors$E0$mu)
  expect_equal(sd$prior_E0_sd, priors$E0$sigma)
  expect_equal(sd$prior_L0_mu, priors$L0$mu)
  expect_equal(sd$prior_L0_sd, priors$L0$sigma)
  expect_equal(sd$prior_sigma_L_sd, priors$sigma_L$sigma)
})

test_that("prior_to_stan_data_hierarchical adds hierarchical fields", {
  priors <- prior_default("hierarchical")
  sd <- prior_to_stan_data_hierarchical(priors)

  expect_true("prior_mu_log_p_Am_mu" %in% names(sd))
  expect_true("prior_mu_log_p_Am_sd" %in% names(sd))
  expect_true("prior_sigma_log_p_Am_rate" %in% names(sd))
  expect_equal(sd$prior_mu_log_p_Am_mu, priors$mu_log_p_Am$mu)
  expect_equal(sd$prior_sigma_log_p_Am_rate, priors$sigma_log_p_Am$rate)
  # Hierarchical model does not use individual p_Am prior
  expect_false("prior_p_Am_mu" %in% names(sd))
  # But includes other base priors
  expect_true("prior_p_M_mu" %in% names(sd))
})

test_that("prior_to_stan_data_growth_repro adds repro fields", {
  priors <- prior_default("growth_repro")
  sd <- prior_to_stan_data_growth_repro(priors)

  expect_true("prior_k_J_mu" %in% names(sd))
  expect_true("prior_k_J_sd" %in% names(sd))
  expect_true("prior_E_Hp_mu" %in% names(sd))
  expect_true("prior_phi_R_mu" %in% names(sd))
  expect_true("prior_k_R_sd" %in% names(sd))
})

test_that("prior_to_stan_data_debtox adds tox fields", {
  priors <- prior_default("debtox")
  sd <- prior_to_stan_data_debtox(priors)

  expect_true("prior_k_d_mu" %in% names(sd))
  expect_true("prior_k_d_sd" %in% names(sd))
  expect_true("prior_z_w_mu" %in% names(sd))
  expect_true("prior_b_w_mu" %in% names(sd))
  expect_equal(sd$prior_k_d_mu, priors$k_d$mu)
})

test_that("prior_to_stan_data all values are scalar numerics", {
  priors <- prior_default("individual")
  sd <- prior_to_stan_data(priors)

  for (nm in names(sd)) {
    expect_true(is.numeric(sd[[nm]]), info = paste("field:", nm))
    expect_true(length(sd[[nm]]) == 1, info = paste("field:", nm))
  }
})
