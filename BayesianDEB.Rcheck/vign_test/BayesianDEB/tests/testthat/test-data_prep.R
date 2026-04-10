# ===========================================================
# Tests: data_prep.R  (bdeb_data, repro_to_intervals, build_stan_data_*)
# ===========================================================

# --- bdeb_data: basic happy path ---

test_that("bdeb_data validates growth data", {
  df_good <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df_good)

  expect_s3_class(dat, "bdeb_data")
  expect_equal(dat$n_ind, 1)
  expect_equal(dat$endpoints, "growth")

  # Missing columns
  df_bad <- data.frame(id = 1, time = 0:5)
  expect_error(bdeb_data(growth = df_bad), "missing columns")

  # No data at all
  expect_error(bdeb_data(), "At least one")
})

test_that("bdeb_data handles multiple individuals", {
  df <- data.frame(
    id = rep(1:5, each = 10),
    time = rep(0:9, 5),
    length = runif(50, 0.1, 0.5)
  )
  dat <- bdeb_data(growth = df)
  expect_equal(dat$n_ind, 5)
  expect_equal(length(dat$ids), 5)
})

test_that("bdeb_data validates reproduction data", {
  df <- data.frame(id = 1, t_start = 0, t_end = 28, count = 50)
  dat <- bdeb_data(reproduction = df)
  expect_equal(dat$endpoints, "reproduction")

  # Bad interval
  df_bad <- data.frame(id = 1, t_start = 28, t_end = 0, count = 50)
  expect_error(bdeb_data(reproduction = df_bad), "t_end")
})

test_that("repro_to_intervals works", {
  cumul <- data.frame(
    id = rep(1, 5),
    time = c(0, 7, 14, 21, 28),
    cumulative = c(0, 10, 30, 60, 100)
  )
  intervals <- repro_to_intervals(cumul)

  expect_equal(nrow(intervals), 4)
  expect_equal(intervals$count, c(10, 20, 30, 40))
  expect_equal(intervals$t_start, c(0, 7, 14, 21))
  expect_equal(intervals$t_end, c(7, 14, 21, 28))
})

test_that("print.bdeb_data works", {
  df <- data.frame(id = 1, time = 0:5, length = seq(0.1, 0.6, by = 0.1))
  dat <- bdeb_data(growth = df)
  # cli output goes to stderr; just check it runs without error
  expect_invisible(print(dat))
})


# --- bdeb_data: edge cases ---

test_that("bdeb_data rejects negative times", {
  df <- data.frame(id = 1, time = c(-1, 0, 1), length = c(0.1, 0.2, 0.3))
  expect_error(bdeb_data(growth = df), "non-negative")
})

test_that("bdeb_data rejects negative lengths", {
  df <- data.frame(id = 1, time = 0:2, length = c(0.1, -0.2, 0.3))
  expect_error(bdeb_data(growth = df), "non-negative")
})

test_that("bdeb_data sorts by id and time", {
  df <- data.frame(
    id     = c(2, 1, 2, 1),
    time   = c(5, 3, 0, 0),
    length = c(0.4, 0.3, 0.1, 0.1)
  )
  dat <- bdeb_data(growth = df)
  expect_equal(dat$growth$id, c(1, 1, 2, 2))
  expect_equal(dat$growth$time, c(0, 3, 0, 5))
})

test_that("bdeb_data handles single observation per individual", {
  df <- data.frame(id = 1, time = 0, length = 0.1)
  dat <- bdeb_data(growth = df)
  expect_equal(dat$n_ind, 1)
  expect_equal(nrow(dat$growth), 1)
})

test_that("bdeb_data works with character IDs", {
  df <- data.frame(
    id     = c("A", "A", "B", "B"),
    time   = c(0, 7, 0, 7),
    length = c(0.1, 0.2, 0.1, 0.3)
  )
  dat <- bdeb_data(growth = df)
  expect_equal(dat$n_ind, 2)
  expect_true(all(c("A", "B") %in% dat$ids))
})

test_that("bdeb_data stores f_food correctly", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  dat <- bdeb_data(growth = df, f_food = 0.7)
  expect_equal(dat$f_food, 0.7)
})

test_that("bdeb_data stores concentration", {
  df <- data.frame(id = 1, time = 0:3, length = 1:4 * 0.1)
  conc <- c("1" = 0)
  dat <- bdeb_data(growth = df, concentration = conc)
  expect_equal(dat$concentration, conc)
})

test_that("bdeb_data unifies IDs across endpoints", {
  g <- data.frame(id = c(1, 1), time = c(0, 7), length = c(0.1, 0.2))
  r <- data.frame(id = c(2), t_start = 0, t_end = 28, count = 50)
  dat <- bdeb_data(growth = g, reproduction = r)
  expect_equal(dat$n_ind, 2)
  expect_true(all(c("1", "2") %in% dat$ids))
  expect_equal(sort(dat$endpoints), c("growth", "reproduction"))
})

test_that("bdeb_data rejects non-data.frame growth input gracefully", {
  expect_error(bdeb_data(growth = list(id = 1, time = 0, length = 0.1)))
})

test_that("bdeb_data accepts zero-length observations", {
  df <- data.frame(id = 1, time = 0:2, length = c(0.0, 0.1, 0.2))
  dat <- bdeb_data(growth = df)
  expect_equal(nrow(dat$growth), 3)
})


# --- Reproduction validation: edge cases ---

test_that("bdeb_data rejects equal t_start and t_end", {
  df <- data.frame(id = 1, t_start = 5, t_end = 5, count = 10)
  expect_error(bdeb_data(reproduction = df), "t_end")
})

test_that("bdeb_data rejects negative reproduction counts", {
  df <- data.frame(id = 1, t_start = 0, t_end = 7, count = -5)
  expect_error(bdeb_data(reproduction = df), "non-negative")
})

test_that("bdeb_data reproduction missing columns gives informative error", {
  df <- data.frame(id = 1, t_start = 0, count = 10)
  expect_error(bdeb_data(reproduction = df), "missing columns")
})




# --- repro_to_intervals: edge cases ---

test_that("repro_to_intervals rejects missing columns", {
  df <- data.frame(id = 1, time = 0:3)
  expect_error(repro_to_intervals(df), "missing columns")
})

test_that("repro_to_intervals returns NULL for single-row individual", {
  df <- data.frame(id = 1, time = 0, cumulative = 0)
  result <- repro_to_intervals(df)
  expect_null(result)
})

test_that("repro_to_intervals handles multiple individuals", {
  df <- data.frame(
    id         = c(1, 1, 1, 2, 2, 2),
    time       = c(0, 7, 14, 0, 7, 14),
    cumulative = c(0, 10, 25, 0, 5, 20)
  )
  result <- repro_to_intervals(df)
  expect_equal(nrow(result), 4)
  expect_true(all(c(1, 2) %in% result$id))
  # Individual 1: 10, 15
  r1 <- result[result$id == 1, ]
  expect_equal(r1$count, c(10, 15))
  # Individual 2: 5, 15
  r2 <- result[result$id == 2, ]
  expect_equal(r2$count, c(5, 15))
})

test_that("repro_to_intervals handles unsorted input", {
  df <- data.frame(
    id         = c(1, 1, 1),
    time       = c(14, 0, 7),
    cumulative = c(30, 0, 10)
  )
  result <- repro_to_intervals(df)
  expect_equal(result$count, c(10, 20))
  expect_equal(result$t_start, c(0, 7))
  expect_equal(result$t_end, c(7, 14))
})

test_that("repro_to_intervals handles zero differences", {
  df <- data.frame(
    id = 1, time = c(0, 7, 14),
    cumulative = c(0, 0, 10)
  )
  result <- repro_to_intervals(df)
  expect_equal(result$count, c(0, 10))
})


# --- build_stan_data_individual: edge cases ---

test_that("build_stan_data_individual rejects multiple individuals", {
  df <- data.frame(
    id = rep(1:2, each = 5), time = rep(0:4, 2),
    length = runif(10, 0.1, 0.5)
  )
  dat <- bdeb_data(growth = df)
  priors <- prior_default("individual")
  expect_error(
    build_stan_data_individual(dat, priors),
    "single individual"
  )
})

test_that("build_stan_data_individual produces correct structure", {
  df <- data.frame(id = 1, time = c(0, 7, 14), length = c(0.1, 0.2, 0.3))
  dat <- bdeb_data(growth = df, f_food = 0.8)
  priors <- prior_default("individual")
  sd <- build_stan_data_individual(dat, priors)

  expect_equal(sd$N_obs, 3)
  expect_equal(sd$t_obs, c(1e-3, 7, 14))  # t=0 replaced with epsilon
  expect_equal(sd$L_obs, c(0.1, 0.2, 0.3))
  expect_equal(sd$f_food, 0.8)
  # Prior hyperparameters present
  expect_true("prior_p_Am_mu" %in% names(sd))
  expect_true("prior_kappa_a" %in% names(sd))
  expect_true("prior_sigma_L_sd" %in% names(sd))
})


# --- build_stan_data_hierarchical: edge cases ---

test_that("build_stan_data_hierarchical pads matrices correctly", {
  # Unequal observation counts
  df <- data.frame(
    id     = c(1, 1, 1, 2, 2),
    time   = c(0, 7, 14, 0, 7),
    length = c(0.1, 0.2, 0.3, 0.1, 0.2)
  )
  dat <- bdeb_data(growth = df)
  priors <- prior_default("hierarchical")
  sd <- build_stan_data_hierarchical(dat, priors)

  expect_equal(sd$N_ind, 2)
  expect_equal(sd$max_N_obs, 3)
  expect_equal(as.integer(sd$N_obs), c(3L, 2L))
  expect_true(is.matrix(sd$t_obs))
  expect_equal(nrow(sd$t_obs), 2)
  expect_equal(ncol(sd$t_obs), 3)
  # Padded entry should be 0 for time
  expect_equal(sd$t_obs[2, 3], 0)
  # Padded entry should be NaN for length
  expect_true(is.nan(sd$L_obs[2, 3]))
})


# --- build_stan_data_debtox: edge cases ---

test_that("build_stan_data_debtox maps reproduction to concentration groups", {
  # 2 conc groups, growth + reproduction
  g <- data.frame(
    id = rep(1:4, each = 3),
    time = rep(c(0, 14, 28), 4),
    length = runif(12, 0.1, 0.5)
  )
  r <- data.frame(
    id      = 1:4,
    t_start = rep(0, 4),
    t_end   = rep(28, 4),
    count   = c(50, 60, 30, 10)
  )
  conc <- c("1" = 0, "2" = 0, "3" = 100, "4" = 100)
  dat <- bdeb_data(growth = g, reproduction = r, concentration = conc)
  priors <- prior_default("debtox")
  sd <- build_stan_data_debtox(dat, priors)

  expect_equal(sd$has_repro, 1L)
  expect_equal(sd$N_groups, 2)
  # Group 1 (conc=0): ids 1,2 => 2 repro records
  # Group 2 (conc=100): ids 3,4 => 2 repro records
  expect_equal(as.integer(sd$N_R), c(2L, 2L))
  expect_equal(sd$max_N_R, 2L)
  # R_counts matrix should have correct values
  expect_true(all(sd$R_counts[1, 1:2] %in% c(50, 60)))
  expect_true(all(sd$R_counts[2, 1:2] %in% c(30, 10)))
  # idx_R_end should point to valid time indices
  expect_true(all(sd$idx_R_end[, 1:2] >= 1))
  expect_true(all(sd$idx_R_end[, 1:2] <= sd$max_N_obs))
})

test_that("build_stan_data_debtox without reproduction sets has_repro=0", {
  g <- data.frame(
    id = rep(1:2, each = 3),
    time = rep(c(0, 7, 14), 2),
    length = runif(6, 0.1, 0.5)
  )
  conc <- c("1" = 0, "2" = 50)
  dat <- bdeb_data(growth = g, concentration = conc)
  priors <- prior_default("debtox")
  sd <- build_stan_data_debtox(dat, priors)

  expect_equal(sd$has_repro, 0L)
  expect_equal(as.integer(sd$N_R), c(0L, 0L))
})

test_that("build_stan_data_debtox handles repro with unmatched ids gracefully", {
  g <- data.frame(id = rep(1:2, each = 3), time = rep(0:2, 2),
                  length = runif(6))
  # Repro for id=99 which doesn't have a concentration mapping
  r <- data.frame(id = 99, t_start = 0, t_end = 2, count = 50)
  conc <- c("1" = 0, "2" = 50)
  dat <- bdeb_data(growth = g, reproduction = r, concentration = conc)
  priors <- prior_default("debtox")
  sd <- build_stan_data_debtox(dat, priors)

  # Unmatched repro id should be dropped, has_repro stays 0
  expect_equal(sd$has_repro, 0L)
})

test_that("build_stan_data_debtox aggregates multi-individual groups with warning", {
  g <- data.frame(
    id = rep(1:4, each = 3),
    time = rep(c(0, 14, 28), 4),
    concentration = rep(c(0, 0, 100, 100), each = 3),
    length = c(0.1, 0.2, 0.3, 0.11, 0.21, 0.31,
               0.1, 0.15, 0.18, 0.09, 0.14, 0.17)
  )
  conc <- c("1" = 0, "2" = 0, "3" = 100, "4" = 100)
  dat <- bdeb_data(growth = g, concentration = conc)
  priors <- prior_default("debtox")

  expect_warning(
    sd <- build_stan_data_debtox(dat, priors),
    "multiple individuals"
  )
  # After aggregation: 2 groups, 3 time points each
  expect_equal(sd$N_groups, 2)
  expect_equal(as.integer(sd$N_obs), c(3L, 3L))
})

test_that("build_stan_data_debtox single-individual groups produce no warning", {
  g <- data.frame(
    id = rep(1:2, each = 3),
    time = rep(c(0, 7, 14), 2),
    concentration = rep(c(0, 50), each = 3),
    length = runif(6, 0.1, 0.5)
  )
  conc <- c("1" = 0, "2" = 50)
  dat <- bdeb_data(growth = g, concentration = conc)
  priors <- prior_default("debtox")

  expect_silent(build_stan_data_debtox(dat, priors))
})

test_that("build_stan_data_debtox rejects missing concentration", {
  df <- data.frame(id = 1, time = 0:3, length = runif(4))
  dat <- bdeb_data(growth = df)
  priors <- prior_default("debtox")
  expect_error(build_stan_data_debtox(dat, priors), "concentration")
})
