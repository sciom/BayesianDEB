# ===========================================================
# Tests: plot.R  (plot.bdeb_fit, plot.bdeb_ppc)
# ===========================================================

# --- plot.bdeb_fit ---

test_that("plot.bdeb_fit trace returns ggplot", {
  fit <- mock_bdeb_fit(n_draws = 50)
  p <- plot(fit, type = "trace", pars = c("p_Am", "kappa"))
  expect_s3_class(p, "ggplot")
})

test_that("plot.bdeb_fit posterior returns ggplot", {
  fit <- mock_bdeb_fit(n_draws = 50)
  p <- plot(fit, type = "posterior", pars = c("p_Am", "p_M"))
  expect_s3_class(p, "ggplot")
})

test_that("plot.bdeb_fit trajectory returns ggplot", {
  fit <- mock_bdeb_fit(n_draws = 50)
  p <- plot(fit, type = "trajectory", n_draws = 20)
  expect_s3_class(p, "ggplot")
})

test_that("plot.bdeb_fit uses default pars when NULL", {
  fit <- mock_bdeb_fit(n_draws = 50)
  # Should not error — picks core pars automatically
  p <- plot(fit, type = "trace")
  expect_s3_class(p, "ggplot")
})

test_that("plot.bdeb_fit with debtox uses debtox core pars", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "debtox")
  p <- plot(fit, type = "posterior")
  expect_s3_class(p, "ggplot")
})


# --- plot.bdeb_ppc ---

test_that("plot.bdeb_ppc returns ggplot for growth", {
  fit <- mock_bdeb_fit(n_draws = 50)
  ppc <- bdeb_ppc(fit, type = "growth")
  p <- plot(ppc, n_draws = 20)
  expect_s3_class(p, "ggplot")
})

test_that("plot.bdeb_ppc errors when no data", {
  ppc <- structure(list(type = "growth", model_type = "individual"),
                   class = "bdeb_ppc")
  expect_error(plot(ppc), "No PPC data")
})


# --- get_core_pars (internal) ---

test_that("get_core_pars returns expected parameters for each type", {
  ind  <- BayesianDEB:::get_core_pars("individual")
  gr   <- BayesianDEB:::get_core_pars("growth_repro")
  hier <- BayesianDEB:::get_core_pars("hierarchical")
  tox  <- BayesianDEB:::get_core_pars("debtox")

  expect_true("p_Am" %in% ind)
  expect_true("sigma_L" %in% ind)
  expect_true("k_J" %in% gr)
  expect_true("mu_log_p_Am" %in% hier)
  expect_true("z_w" %in% tox)
  expect_true("b_w" %in% tox)
})
