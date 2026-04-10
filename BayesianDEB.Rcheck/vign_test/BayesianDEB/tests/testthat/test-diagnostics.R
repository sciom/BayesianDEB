# ===========================================================
# Tests: diagnostics.R  (bdeb_diagnose, bdeb_summary, bdeb_derived)
# ===========================================================

# --- bdeb_diagnose ---

test_that("bdeb_diagnose rejects non-bdeb_fit", {
  expect_error(bdeb_diagnose(list()), "bdeb_fit")
  expect_error(bdeb_diagnose(NULL), "bdeb_fit")
  expect_error(bdeb_diagnose("fit"), "bdeb_fit")
})

test_that("bdeb_diagnose works with mock fit", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "individual")
  result <- bdeb_diagnose(fit)

  expect_type(result, "list")
  expect_true("n_divergent" %in% names(result))
  expect_true("n_max_treedepth" %in% names(result))
  expect_true("ebfmi" %in% names(result))
  expect_true("summary" %in% names(result))

  expect_equal(result$n_divergent, 0)
  expect_equal(result$n_max_treedepth, 0)
})

test_that("bdeb_diagnose returns invisible", {
  fit <- mock_bdeb_fit(n_draws = 50)
  expect_invisible(bdeb_diagnose(fit))
})

test_that("bdeb_diagnose with specific pars", {
  fit <- mock_bdeb_fit(n_draws = 100)
  result <- bdeb_diagnose(fit, pars = c("p_Am", "kappa"))

  # Summary should only have 2 parameters
  expect_equal(nrow(result$summary), 2)
  expect_true(all(result$summary$variable %in% c("p_Am", "kappa")))
})

test_that("bdeb_diagnose summary has correct columns", {
  fit <- mock_bdeb_fit(n_draws = 100)
  result <- bdeb_diagnose(fit, pars = "p_Am")

  s <- result$summary
  expect_true("mean" %in% names(s))
  expect_true("sd" %in% names(s))
  expect_true("rhat" %in% names(s))
  expect_true("ess_bulk" %in% names(s))
  expect_true("ess_tail" %in% names(s))
})


# --- bdeb_loo ---

test_that("bdeb_loo rejects non-bdeb_fit", {
  expect_error(bdeb_loo(list()), "bdeb_fit")
})

test_that("bdeb_loo rejects hierarchical model", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "hierarchical")
  expect_error(bdeb_loo(fit), "not available.*hierarchical")
})

test_that("bdeb_loo rejects debtox model", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "debtox")
  expect_error(bdeb_loo(fit), "not available.*debtox")
})

test_that("bdeb_loo works for individual model", {
  skip_if_not_installed("loo")
  fit <- mock_bdeb_fit(n_draws = 100, type = "individual")
  result <- bdeb_loo(fit)
  expect_s3_class(result, "loo")
})


# --- bdeb_summary ---

test_that("bdeb_summary rejects non-bdeb_fit", {
  expect_error(bdeb_summary(list()), "bdeb_fit")
  expect_error(bdeb_summary(42), "bdeb_fit")
})

test_that("bdeb_summary works with mock fit", {
  fit <- mock_bdeb_fit(n_draws = 100)
  s <- bdeb_summary(fit)

  expect_s3_class(s, "draws_summary")
  expect_true(nrow(s) > 0)
  expect_true("mean" %in% names(s))
  expect_true("rhat" %in% names(s))
})

test_that("bdeb_summary filters by pars", {
  fit <- mock_bdeb_fit(n_draws = 100)
  s <- bdeb_summary(fit, pars = c("p_Am", "p_M", "kappa"))

  expect_equal(nrow(s), 3)
  expect_true(all(s$variable %in% c("p_Am", "p_M", "kappa")))
})

test_that("bdeb_summary respects prob argument", {
  fit <- mock_bdeb_fit(n_draws = 200)

  s90 <- suppressWarnings(bdeb_summary(fit, pars = "p_Am", prob = 0.90))
  s50 <- suppressWarnings(bdeb_summary(fit, pars = "p_Am", prob = 0.50))

  # summarise_draws names quantile columns by percentile (e.g. "5%", "95%")
  s90_df <- as.data.frame(s90)
  s50_df <- as.data.frame(s50)

  width90 <- s90_df[["95%"]] - s90_df[["5%"]]
  width50 <- s50_df[["75%"]] - s50_df[["25%"]]
  expect_gt(width90, width50)
})

test_that("bdeb_summary excludes log_lik and lp__ by default", {
  fit <- mock_bdeb_fit(n_draws = 100)
  s <- bdeb_summary(fit)

  expect_false(any(grepl("^log_lik", s$variable)))
  expect_false(any(grepl("^lp__", s$variable)))
})


# --- bdeb_derived ---

test_that("bdeb_derived rejects non-bdeb_fit", {
  expect_error(bdeb_derived(list()), "bdeb_fit")
  expect_error(bdeb_derived(NULL), "bdeb_fit")
})

test_that("bdeb_derived computes L_inf correctly", {
  fit <- mock_bdeb_fit(n_draws = 100)
  d <- bdeb_derived(fit, quantities = "L_inf", f = 1.0)

  expect_s3_class(d, "draws_df")
  expect_true("L_inf" %in% names(d))
  expect_true(all(d$L_inf > 0))
})

test_that("bdeb_derived computes k_M correctly", {
  fit <- mock_bdeb_fit(n_draws = 100)
  d <- bdeb_derived(fit, quantities = "k_M")

  expect_true("k_M" %in% names(d))
  # k_M = p_M / E_G, both positive => k_M positive
  expect_true(all(d$k_M > 0))
})

test_that("bdeb_derived computes growth_rate correctly", {
  fit <- mock_bdeb_fit(n_draws = 100)
  d <- bdeb_derived(fit, quantities = "growth_rate")

  expect_true("growth_rate" %in% names(d))
  expect_true(all(d$growth_rate > 0))
})

test_that("bdeb_derived computes g (investment ratio)", {
  fit <- mock_bdeb_fit(n_draws = 100)
  d <- bdeb_derived(fit, quantities = "g")

  expect_true("g" %in% names(d))
  expect_true(all(d$g > 0))
})

test_that("bdeb_derived computes all three quantities", {
  fit <- mock_bdeb_fit(n_draws = 100)
  d <- bdeb_derived(fit, quantities = c("L_inf", "k_M", "growth_rate"))

  expect_true(all(c("L_inf", "k_M", "growth_rate") %in% names(d)))
  expect_equal(nrow(d), 100)
})

test_that("bdeb_derived L_inf scales with f", {
  fit <- mock_bdeb_fit(n_draws = 100)
  d1 <- bdeb_derived(fit, quantities = "L_inf", f = 1.0)
  d07 <- bdeb_derived(fit, quantities = "L_inf", f = 0.7)

  # L_inf at f=0.7 should be ~70% of L_inf at f=1.0
  expect_true(mean(d07$L_inf) < mean(d1$L_inf))
  ratio <- mean(d07$L_inf) / mean(d1$L_inf)
  expect_true(abs(ratio - 0.7) < 0.01)
})

test_that("bdeb_derived handles hierarchical fit (mu_log_p_Am)", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "hierarchical")
  d <- bdeb_derived(fit, quantities = "L_inf")

  expect_true("L_inf" %in% names(d))
  expect_true(all(d$L_inf > 0))
})

test_that("bdeb_derived with empty quantities returns only .draw", {
  fit <- mock_bdeb_fit(n_draws = 50)
  d <- bdeb_derived(fit, quantities = character(0))

  expect_s3_class(d, "draws_df")
})
