# ===========================================================
# Tests: fit.R  (bdeb_fit input validation, print, summary S3)
# ===========================================================

test_that("bdeb_fit rejects non-bdeb_model input", {
  expect_error(bdeb_fit(list(a = 1)), "bdeb_model")
  expect_error(bdeb_fit(NULL), "bdeb_model")
  expect_error(bdeb_fit(data.frame(x = 1)), "bdeb_model")
})

test_that("bdeb_fit requires cmdstanr", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  mod <- bdeb_model(dat, type = "individual")

  # We can't easily test this without unloading cmdstanr,
  # but we verify the function signature accepts all arguments
  expect_true(is.function(bdeb_fit))
  args <- formals(bdeb_fit)
  expect_true("chains" %in% names(args))
  expect_true("iter_warmup" %in% names(args))
  expect_true("iter_sampling" %in% names(args))
  expect_true("adapt_delta" %in% names(args))
  expect_true("max_treedepth" %in% names(args))
  expect_true("seed" %in% names(args))
  expect_true("parallel_chains" %in% names(args))
  expect_true("refresh" %in% names(args))
})

test_that("print.bdeb_fit works with mock", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  expect_invisible(print(fit))
})

test_that("coef.bdeb_fit returns named numeric vector", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  cc <- coef(fit)
  expect_type(cc, "double")
  expect_true(length(cc) > 0)
  expect_true("p_Am" %in% names(cc))
  expect_true("kappa" %in% names(cc))
  # Should not contain generated quantities
  expect_false(any(grepl("^L_hat|^L_rep|^log_lik|^lp__", names(cc))))
})

test_that("coef.bdeb_fit type='mean' differs from 'median'", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "individual")
  cc_med  <- coef(fit, type = "median")
  cc_mean <- coef(fit, type = "mean")
  expect_equal(length(cc_med), length(cc_mean))
  # Not identical (unless symmetric, which mock data isn't exactly)
  expect_false(identical(cc_med, cc_mean))
})

test_that("summary.bdeb_fit returns a draws_summary", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  s <- summary(fit)
  expect_s3_class(s, "draws_summary")
})


# --- confint.bdeb_fit ---

test_that("confint.bdeb_fit returns a matrix with lower/upper columns", {
  fit <- mock_bdeb_fit(n_draws = 100)
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2)
  expect_true(nrow(ci) > 0)
  expect_true(all(ci[, 1] <= ci[, 2]))
  # Default level = 0.95 -> "2.5%" / "97.5%"
  expect_equal(colnames(ci), c("2.5%", "97.5%"))
})

test_that("confint.bdeb_fit respects level", {
  fit <- mock_bdeb_fit(n_draws = 200)
  ci95 <- confint(fit, level = 0.95)
  ci80 <- confint(fit, level = 0.80)
  width95 <- ci95[, 2] - ci95[, 1]
  width80 <- ci80[, 2] - ci80[, 1]
  expect_true(all(width95 >= width80 - 1e-12))
})

test_that("confint.bdeb_fit filters by parm", {
  fit <- mock_bdeb_fit(n_draws = 100)
  ci <- confint(fit, parm = c("p_Am", "kappa"))
  expect_equal(nrow(ci), 2)
  expect_setequal(rownames(ci), c("p_Am", "kappa"))
})

test_that("confint.bdeb_fit errors on unknown parm", {
  fit <- mock_bdeb_fit(n_draws = 50)
  expect_error(confint(fit, parm = "nonsense"), "Unknown")
})

test_that("confint.bdeb_fit validates level", {
  fit <- mock_bdeb_fit(n_draws = 50)
  expect_error(confint(fit, level = 0),    "level")
  expect_error(confint(fit, level = 1),    "level")
  expect_error(confint(fit, level = c(0.5, 0.9)), "level")
})

test_that("confint.bdeb_fit rejects non-bdeb_fit", {
  expect_error(confint.bdeb_fit(list()), "bdeb_fit")
})


# --- nobs.bdeb_fit ---

test_that("nobs.bdeb_fit returns positive integer", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "individual")
  n <- nobs(fit)
  expect_type(n, "integer")
  expect_gt(n, 0L)
})

test_that("nobs.bdeb_fit rejects non-bdeb_fit", {
  expect_error(nobs.bdeb_fit(list()), "bdeb_fit")
})


# --- fitted.bdeb_fit ---

test_that("fitted.bdeb_fit returns named numeric for individual model", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "individual")
  f <- fitted(fit)
  expect_type(f, "double")
  expect_true(length(f) > 0)
  expect_true(all(grepl("^L_hat\\[", names(f))))
})

test_that("fitted.bdeb_fit median != mean (in general)", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "individual")
  fm <- fitted(fit, type = "median")
  fa <- fitted(fit, type = "mean")
  expect_false(identical(fm, fa))
})

test_that("fitted.bdeb_fit errors for hierarchical model", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "hierarchical")
  expect_error(fitted(fit), "not available")
})

test_that("fitted.bdeb_fit rejects non-bdeb_fit", {
  expect_error(fitted.bdeb_fit(list()), "bdeb_fit")
})


# --- residuals.bdeb_fit ---

test_that("residuals.bdeb_fit returns observed minus fitted", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "individual")
  r <- residuals(fit)
  expect_type(r, "double")
  expect_equal(length(r), length(fitted(fit)))
})

test_that("residuals.bdeb_fit rejects non-bdeb_fit", {
  expect_error(residuals.bdeb_fit(list()), "bdeb_fit")
})


# --- vcov.bdeb_fit ---

test_that("vcov.bdeb_fit returns symmetric matrix on parameters", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "individual")
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), ncol(V))
  expect_equal(rownames(V), colnames(V))
  expect_true(all(diag(V) >= 0))
  # Should not include log_lik / L_hat / lp__
  expect_false(any(grepl("^(log_lik|L_hat|lp__)", rownames(V))))
})

test_that("vcov.bdeb_fit rejects non-bdeb_fit", {
  expect_error(vcov.bdeb_fit(list()), "bdeb_fit")
})


# --- logLik.bdeb_fit ---

test_that("logLik.bdeb_fit returns logLik object with df and nobs", {
  fit <- mock_bdeb_fit(n_draws = 100, type = "individual")
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_type(as.numeric(ll), "double")
  expect_false(is.null(attr(ll, "df")))
  expect_false(is.null(attr(ll, "nobs")))
  expect_gt(attr(ll, "df"), 0L)
  expect_gt(attr(ll, "nobs"), 0L)
})

test_that("logLik.bdeb_fit errors for hierarchical model", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "hierarchical")
  expect_error(logLik(fit), "not available")
})

test_that("logLik.bdeb_fit errors for debtox model", {
  fit <- mock_bdeb_fit(n_draws = 50, type = "debtox")
  expect_error(logLik(fit), "not available")
})

test_that("logLik.bdeb_fit rejects non-bdeb_fit", {
  expect_error(logLik.bdeb_fit(list()), "bdeb_fit")
})
