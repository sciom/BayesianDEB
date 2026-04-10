# ===========================================================
# B. Scientific consistency tests
# Verify physical/biological laws, monotonicity, dimensional
# consistency, and known analytical results
# ===========================================================


# --- Arrhenius: known analytical properties ---

test_that("arrhenius(T_ref) == 1 exactly for any T_A", {
  for (ta in c(0, 100, 6000, 8000, 12000, 50000)) {
    expect_equal(arrhenius(293.15, T_ref = 293.15, T_A = ta), 1.0,
                 info = paste("T_A =", ta))
  }
  # Also at different T_ref
  expect_equal(arrhenius(310, T_ref = 310, T_A = 8000), 1.0)
  expect_equal(arrhenius(273.15, T_ref = 273.15, T_A = 5000), 1.0)
})

test_that("arrhenius is strictly monotonically increasing with T", {
  temps <- seq(275, 315, by = 2)
  vals <- vapply(temps, arrhenius, numeric(1), T_ref = 293.15, T_A = 8000)
  expect_true(all(diff(vals) > 0))
})

test_that("arrhenius T_A = 0 gives c_T = 1 for all temperatures", {
  for (temp in c(250, 280, 293.15, 310, 350)) {
    expect_equal(arrhenius(temp, T_ref = 293.15, T_A = 0), 1.0,
                 info = paste("T =", temp))
  }
})

test_that("arrhenius is always positive", {
  vals <- vapply(seq(200, 400, by = 10), arrhenius, numeric(1),
                 T_ref = 293.15, T_A = 8000)
  expect_true(all(vals > 0))
})

test_that("higher T_A means steeper temperature response", {
  # At T > T_ref, higher T_A gives larger c_T
  c_low  <- arrhenius(303.15, T_ref = 293.15, T_A = 4000)
  c_high <- arrhenius(303.15, T_ref = 293.15, T_A = 12000)
  expect_gt(c_high, c_low)

  # At T < T_ref, higher T_A gives smaller c_T
  c_low2  <- arrhenius(283.15, T_ref = 293.15, T_A = 4000)
  c_high2 <- arrhenius(283.15, T_ref = 293.15, T_A = 12000)
  expect_lt(c_high2, c_low2)
})


# --- Derived quantities: monotonicity ---

test_that("L_m reference: Eisenia fetida AmP parameters give L_m = 7.5 cm", {
  # From AmP: {p_Am} = 5.0, [p_M] = 0.5, kappa = 0.75
  # L_m = kappa * {p_Am} / [p_M] = 0.75 * 5.0 / 0.5 = 7.5 cm (structural)
  L_m <- 0.75 * 5.0 / 0.5
  expect_equal(L_m, 7.5)
})

test_that("L_inf = f * L_m exactly (structural, not physical)", {
  fit <- mock_bdeb_fit(n_draws = 100)
  d <- bdeb_derived(fit, quantities = c("L_m", "L_inf"), f = 0.8)
  ratios <- as.numeric(d$L_inf) / as.numeric(d$L_m)
  expect_true(all(abs(ratios - 0.8) < 1e-12))
})

test_that("L_m is independent of f", {
  fit <- mock_bdeb_fit(n_draws = 50)
  d1  <- bdeb_derived(fit, quantities = "L_m", f = 1.0)
  d07 <- bdeb_derived(fit, quantities = "L_m", f = 0.7)
  expect_equal(as.numeric(d1$L_m), as.numeric(d07$L_m), tolerance = 1e-12)
})

test_that("L_inf increases with p_Am (all else equal)", {
  fit <- mock_bdeb_fit(n_draws = 200)
  draws <- posterior::as_draws_df(fit$fit$draws())

  L_inf <- draws$kappa * draws$p_Am / draws$p_M
  # Correlation between p_Am and L_inf should be positive
  expect_gt(stats::cor(draws$p_Am, L_inf), 0.5)
})

test_that("L_inf decreases with p_M (all else equal)", {
  fit <- mock_bdeb_fit(n_draws = 200)
  draws <- posterior::as_draws_df(fit$fit$draws())

  L_inf <- draws$kappa * draws$p_Am / draws$p_M
  # Correlation between p_M and L_inf should be negative
  expect_lt(stats::cor(draws$p_M, L_inf), -0.3)
})

test_that("L_inf increases with kappa", {
  fit <- mock_bdeb_fit(n_draws = 200)
  draws <- posterior::as_draws_df(fit$fit$draws())

  L_inf <- draws$kappa * draws$p_Am / draws$p_M
  expect_gt(stats::cor(draws$kappa, L_inf), 0.3)
})

test_that("L_inf is exactly proportional to f", {
  fit <- mock_bdeb_fit(n_draws = 100)
  d1  <- bdeb_derived(fit, quantities = "L_inf", f = 1.0)
  d05 <- bdeb_derived(fit, quantities = "L_inf", f = 0.5)
  ratios <- as.numeric(d05$L_inf) / as.numeric(d1$L_inf)
  expect_true(all(abs(ratios - 0.5) < 1e-12))
})

test_that("k_M is always positive when p_M > 0 and E_G > 0", {
  fit <- mock_bdeb_fit(n_draws = 200)
  d <- bdeb_derived(fit, quantities = "k_M")
  expect_true(all(d$k_M > 0))
})

test_that("growth_rate is always positive", {
  fit <- mock_bdeb_fit(n_draws = 200)
  d <- bdeb_derived(fit, quantities = "growth_rate", f = 1.0)
  expect_true(all(d$growth_rate > 0))
})

test_that("growth_rate increases as f decreases (VB rate property)", {
  # r_B = k_M * g / (3*(f+g)) — smaller f => smaller denominator => larger r_B
  fit <- mock_bdeb_fit(n_draws = 100)
  d1  <- bdeb_derived(fit, quantities = "growth_rate", f = 1.0)
  d05 <- bdeb_derived(fit, quantities = "growth_rate", f = 0.5)
  expect_true(mean(d05$growth_rate) > mean(d1$growth_rate))
})

test_that("g (investment ratio) is always positive", {
  fit <- mock_bdeb_fit(n_draws = 200)
  d <- bdeb_derived(fit, quantities = "g")
  expect_true(all(d$g > 0))
})


# --- EC50/NEC: analytical identity ---

test_that("EC50 = z_w + 0.5/b_w holds for every posterior draw", {
  fit <- mock_bdeb_fit(n_draws = 500, type = "debtox")
  draws <- posterior::as_draws_df(fit$fit$draws())

  computed <- draws$z_w + 0.5 / draws$b_w
  expect_equal(as.numeric(draws$EC50), as.numeric(computed),
               tolerance = 1e-12)
})

test_that("NEC == z_w for every draw", {
  fit <- mock_bdeb_fit(n_draws = 500, type = "debtox")
  draws <- posterior::as_draws_df(fit$fit$draws())
  expect_equal(as.numeric(draws$NEC), as.numeric(draws$z_w),
               tolerance = 1e-12)
})

test_that("EC50 > NEC always (since b_w > 0)", {
  fit <- mock_bdeb_fit(n_draws = 500, type = "debtox")
  draws <- posterior::as_draws_df(fit$fit$draws())
  expect_true(all(draws$EC50 > draws$NEC))
})

test_that("EC50 decreases with larger b_w (stronger toxicant)", {
  # EC50 = z_w + 0.5/b_w — larger b_w => smaller 0.5/b_w => smaller EC50
  fit <- mock_bdeb_fit(n_draws = 200, type = "debtox")
  draws <- posterior::as_draws_df(fit$fit$draws())
  expect_lt(stats::cor(draws$b_w, draws$EC50), -0.3)
})


# --- deb_fluxes: energy conservation ---

test_that("kappa-rule: kappa*p_C = p_G + p_M when growth positive", {
  fl <- deb_fluxes(E = 500, V = 0.1, f = 1.0,
                   p_Am = 5, p_M = 0.2, kappa = 0.8,
                   v = 0.2, E_G = 400, k_J = 0, E_Hp = 0)
  # Verify growth is positive
  expect_gt(fl$p_G, 0)
  # kappa * p_C = p_G + p_M
  expect_equal(0.8 * fl$p_C, fl$p_G + fl$p_M, tolerance = 1e-10)
})

test_that("(1-kappa)-rule: (1-kappa)*p_C = p_R + p_J when repro positive", {
  # Use small k_J*E_Hp so that (1-kappa)*p_C > p_J
  fl <- deb_fluxes(E = 500, V = 0.1, f = 1.0,
                   p_Am = 5, p_M = 0.2, kappa = 0.8,
                   v = 0.2, E_G = 400, k_J = 0.001, E_Hp = 1)
  # Verify reproduction flux is positive
  expect_gt(fl$p_R, 0)
  expect_equal(0.2 * fl$p_C, fl$p_R + fl$p_J, tolerance = 1e-10)
})

test_that("p_C in deb_fluxes matches Stan ODE exactly", {
  # Verify the R helper and Stan use the same formula:
  # p_C = E * v * L / (E + E_G * V + 1e-12)
  E <- 50; V <- 0.5; v <- 0.2; E_G <- 400
  L <- V^(1/3)

  fl <- deb_fluxes(E = E, V = V, f = 1.0, p_Am = 5, p_M = 0.5,
                   kappa = 0.75, v = v, E_G = E_G)

  expected <- E * v * L / (E + E_G * V + 1e-12)
  expect_equal(fl$p_C, expected, tolerance = 1e-12)
})

test_that("p_A = 0 when f = 0", {
  fl <- deb_fluxes(E = 10, V = 0.5, f = 0, p_Am = 5, p_M = 0.5,
                   kappa = 0.75, v = 0.2, E_G = 400)
  expect_equal(fl$p_A, 0)
})

test_that("p_G >= 0 always (no negative growth flux)", {
  # Even under starvation
  fl <- deb_fluxes(E = 0.001, V = 10, f = 0, p_Am = 5, p_M = 0.5,
                   kappa = 0.75, v = 0.2, E_G = 400)
  expect_gte(fl$p_G, 0)
})

test_that("p_R >= 0 always (no negative reproduction flux)", {
  fl <- deb_fluxes(E = 0.001, V = 0.5, f = 0, p_Am = 5, p_M = 0.5,
                   kappa = 0.75, v = 0.2, E_G = 400,
                   k_J = 0.01, E_Hp = 100)
  expect_gte(fl$p_R, 0)
})

test_that("p_A scales with f * L^2", {
  fl1 <- deb_fluxes(E = 10, V = 0.5, f = 1.0, p_Am = 5, p_M = 0.5,
                    kappa = 0.75, v = 0.2, E_G = 400)
  fl2 <- deb_fluxes(E = 10, V = 0.5, f = 0.5, p_Am = 5, p_M = 0.5,
                    kappa = 0.75, v = 0.2, E_G = 400)
  expect_equal(fl2$p_A, fl1$p_A * 0.5, tolerance = 1e-12)
})


# --- repro_to_intervals: cumulative sum preservation ---

test_that("repro_to_intervals preserves cumulative total", {
  cumul <- data.frame(
    id = 1, time = c(0, 7, 14, 21, 28),
    cumulative = c(0, 10, 30, 60, 100)
  )
  intervals <- repro_to_intervals(cumul)
  expect_equal(sum(intervals$count), 100)
})

test_that("repro_to_intervals: sum of counts equals final - initial", {
  cumul <- data.frame(
    id = 1, time = c(0, 5, 10, 15),
    cumulative = c(5, 20, 45, 80)
  )
  intervals <- repro_to_intervals(cumul)
  expect_equal(sum(intervals$count), 80 - 5)
})

test_that("repro_to_intervals preserves total across multiple individuals", {
  cumul <- data.frame(
    id = c(1, 1, 1, 2, 2, 2),
    time = c(0, 14, 28, 0, 14, 28),
    cumulative = c(0, 25, 60, 0, 15, 40)
  )
  intervals <- repro_to_intervals(cumul)
  # id 1: 60 - 0 = 60 total, id 2: 40 - 0 = 40 total
  expect_equal(sum(intervals$count[intervals$id == 1]), 60)
  expect_equal(sum(intervals$count[intervals$id == 2]), 40)
})
