# ===========================================================
# Tests: simulate.R (deb_simulate, debtox_simulate)
# ===========================================================

test_that("deb_simulate returns correct structure", {
  traj <- deb_simulate(t_max = 84, p_Am = 5, p_M = 0.5,
    kappa = 0.75, v = 0.2, E_G = 400, E0 = 1, L0 = 0.1)

  expect_s3_class(traj, "data.frame")
  expect_named(traj, c("time", "E", "V", "L"))
  expect_equal(traj$time[1], 0)
  expect_true(max(traj$time) >= 84)
  expect_true(all(traj$L > 0))
  expect_true(all(traj$E > 0))
  expect_true(all(traj$V > 0))
})

test_that("deb_simulate L grows over time under ad libitum", {
  traj <- deb_simulate(84, 5, 0.5, 0.75, 0.2, 400, 1, 0.1, f = 1, dt = 0.1)
  # Final > initial (overall growth)
  expect_gt(tail(traj$L, 1), traj$L[1])
})

test_that("deb_simulate final L increases with longer simulation", {
  t1 <- deb_simulate(84, 5, 0.5, 0.75, 0.2, 400, 1, 0.1, dt = 0.1)
  t2 <- deb_simulate(840, 5, 0.5, 0.75, 0.2, 400, 1, 0.1, dt = 0.1)
  expect_gt(tail(t2$L, 1), tail(t1$L, 1))
})

test_that("deb_simulate f = 0 means no growth", {
  traj <- deb_simulate(84, 5, 0.5, 0.75, 0.2, 400, 1, 0.1, f = 0)
  # Length should not increase much (only initial transient)
  expect_lt(max(traj$L) - traj$L[1], 0.01)
})

test_that("deb_simulate validates inputs", {
  expect_error(deb_simulate(-1, 5, 0.5, 0.75, 0.2, 400, 1, 0.1), "positive")
  expect_error(deb_simulate(84, 0, 0.5, 0.75, 0.2, 400, 1, 0.1), "positive")
})

test_that("debtox_simulate returns correct structure", {
  traj <- debtox_simulate(42, 5, 0.5, 0.75, 0.2, 400, 1, 0.1,
    k_d = 0.3, z_w = 15, b_w = 0.003, C_w = 80)

  expect_named(traj, c("time", "E", "V", "L", "R", "Dw"))
  expect_true(all(traj$L > 0))
  expect_true(all(traj$Dw >= 0))
})

test_that("debtox_simulate C_w = 0 equals deb_simulate", {
  t1 <- deb_simulate(84, 5, 0.5, 0.75, 0.2, 400, 1, 0.1)
  t2 <- debtox_simulate(84, 5, 0.5, 0.75, 0.2, 400, 1, 0.1,
    k_d = 0.3, z_w = 15, b_w = 0.003, C_w = 0)
  expect_equal(t1$L, t2$L, tolerance = 1e-10)
})

test_that("debtox_simulate higher C_w gives shorter final length", {
  t_ctrl <- debtox_simulate(42, 5, 0.5, 0.75, 0.2, 400, 1, 0.1,
    k_d = 0.3, z_w = 15, b_w = 0.003, C_w = 0)
  t_high <- debtox_simulate(42, 5, 0.5, 0.75, 0.2, 400, 1, 0.1,
    k_d = 0.3, z_w = 15, b_w = 0.003, C_w = 200)
  expect_gt(tail(t_ctrl$L, 1), tail(t_high$L, 1))
})

test_that("debtox_simulate C_w below NEC has no effect", {
  t_ctrl <- debtox_simulate(42, 5, 0.5, 0.75, 0.2, 400, 1, 0.1,
    k_d = 0.3, z_w = 50, b_w = 0.003, C_w = 0)
  t_below <- debtox_simulate(42, 5, 0.5, 0.75, 0.2, 400, 1, 0.1,
    k_d = 0.3, z_w = 50, b_w = 0.003, C_w = 30)
  # C_w=30 < z_w=50 => no damage => same trajectory
  expect_equal(t_ctrl$L, t_below$L, tolerance = 1e-10)
})
