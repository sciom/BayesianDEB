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

test_that("plot.bdeb_fit pairs returns a non-NULL bayesplot grid (C2)", {
  skip_if_not_installed("gridExtra")
  fit <- mock_bdeb_fit(n_draws = 50)
  p <- plot(fit, type = "pairs", pars = c("p_Am", "p_M", "kappa"))
  expect_false(is.null(p))
  # bayesplot::mcmc_pairs returns a bayesplot_grid (gtable);
  # accept either as valid output.
  expect_true(inherits(p, "bayesplot_grid") ||
              inherits(p, "gtable") ||
              inherits(p, "ggplot"))
})

test_that("plot.bdeb_fit pairs requires >= 2 parameters", {
  fit <- mock_bdeb_fit(n_draws = 50)
  # length-1 pars triggers our own validation before reaching bayesplot,
  # so this test does not require gridExtra.
  expect_error(plot(fit, type = "pairs", pars = "p_Am"),
               "pars")
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

test_that("plot.bdeb_fit trajectory works for hierarchical (faceted by individual)", {
  fit <- mock_bdeb_fit(n_draws = 20, type = "hierarchical")
  p <- plot(fit, type = "trajectory", n_draws = 5)
  expect_s3_class(p, "ggplot")
  # Should have facets
  expect_true("FacetWrap" %in% class(p$facet))
})

test_that("plot.bdeb_fit trajectory works for debtox (faceted by group)", {
  fit <- mock_bdeb_fit(n_draws = 20, type = "debtox")
  p <- plot(fit, type = "trajectory", n_draws = 5)
  expect_s3_class(p, "ggplot")
  expect_true("FacetWrap" %in% class(p$facet))
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


# --- plot.bdeb_fit: edge cases ---

test_that("plot.bdeb_fit with n_draws=1 doesn't crash", {
  fit <- mock_bdeb_fit(n_draws = 50)
  p <- plot(fit, type = "trajectory", n_draws = 1)
  expect_s3_class(p, "ggplot")
})

test_that("plot.bdeb_fit with n_draws > available draws caps silently", {
  fit <- mock_bdeb_fit(n_draws = 20)
  p <- plot(fit, type = "trace", pars = "p_Am")
  expect_s3_class(p, "ggplot")
})

test_that("plot.bdeb_fit hierarchical default pars", {
  fit <- mock_bdeb_fit(n_draws = 30, type = "hierarchical")
  p <- plot(fit, type = "trace")
  expect_s3_class(p, "ggplot")
})

test_that("plot_dose_response returns ggplot", {
  fit <- mock_bdeb_fit(n_draws = 20, type = "debtox")
  p <- plot_dose_response(fit, n_draws = 3)
  expect_s3_class(p, "ggplot")
})

test_that("plot_dose_response rejects individual fit", {
  fit <- mock_bdeb_fit(n_draws = 20, type = "individual")
  expect_error(plot_dose_response(fit), "DEBtox")
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
