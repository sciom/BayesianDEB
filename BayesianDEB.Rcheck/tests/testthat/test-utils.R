# ===========================================================
# Tests: utils.R  (arrhenius, deb_fluxes, check_cmdstanr, assert_positive)
# ===========================================================

# --- arrhenius: happy path ---

test_that("arrhenius correction works", {
  # At reference temperature, correction should be 1
  expect_equal(arrhenius(293.15, T_ref = 293.15, T_A = 8000), 1.0)

  # Higher temperature => factor > 1
  expect_gt(arrhenius(303.15, T_ref = 293.15, T_A = 8000), 1.0)

  # Lower temperature => factor < 1
  expect_lt(arrhenius(283.15, T_ref = 293.15, T_A = 8000), 1.0)
})


# --- arrhenius: numerical properties ---

test_that("arrhenius is monotonically increasing with temperature", {
  temps <- seq(273.15, 313.15, by = 5)
  factors <- vapply(temps, arrhenius, numeric(1),
                    T_ref = 293.15, T_A = 8000)
  # Each subsequent factor should be larger
  expect_true(all(diff(factors) > 0))
})

test_that("arrhenius scales with T_A", {
  # Higher T_A means steeper temperature dependence
  low  <- arrhenius(303.15, T_ref = 293.15, T_A = 4000)
  high <- arrhenius(303.15, T_ref = 293.15, T_A = 12000)
  expect_gt(high, low)
})

test_that("arrhenius handles extreme temperatures", {
  # Very cold but above 0 K
  val_cold <- arrhenius(200, T_ref = 293.15, T_A = 8000)
  expect_true(is.finite(val_cold))
  expect_gt(val_cold, 0)

  # Very hot
  val_hot <- arrhenius(400, T_ref = 293.15, T_A = 8000)
  expect_true(is.finite(val_hot))
  expect_gt(val_hot, 1)
})

test_that("arrhenius with T_A = 0 returns 1 for any temperature", {
  expect_equal(arrhenius(250, T_ref = 293.15, T_A = 0), 1.0)
  expect_equal(arrhenius(350, T_ref = 293.15, T_A = 0), 1.0)
})

test_that("arrhenius with different T_ref shifts reference point", {
  # Factor should be 1 at whatever T_ref is used
  expect_equal(arrhenius(310, T_ref = 310, T_A = 8000), 1.0)
})

test_that("arrhenius is symmetric around T_ref", {
  # c_T(T_ref + delta) * c_T(T_ref - delta) is NOT 1 (exponential, not linear)
  # but c_T(T_ref + d) > 1 and c_T(T_ref - d) < 1
  up   <- arrhenius(298.15, T_ref = 293.15, T_A = 8000)
  down <- arrhenius(288.15, T_ref = 293.15, T_A = 8000)
  expect_gt(up, 1)
  expect_lt(down, 1)
})


# --- deb_fluxes: happy path ---

test_that("deb_fluxes returns correct structure", {
  fl <- deb_fluxes(E = 10, V = 0.5, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)

  expect_type(fl, "list")
  expect_named(fl, c("p_A", "p_C", "p_M", "p_G", "p_J", "p_R", "L", "e"))

  # All fluxes should be numeric
  for (nm in names(fl)) {
    expect_type(fl[[nm]], "double")
  }

  # L should be V^(1/3)
  expect_equal(fl$L, 0.5^(1/3))

  # p_A should be f * p_Am * L^2
  expect_equal(fl$p_A, 1.0 * 5 * (0.5^(1/3))^2)
})


# --- deb_fluxes: edge cases and numerical stability ---

test_that("deb_fluxes with zero reserve does not crash", {
  fl <- deb_fluxes(E = 0, V = 0.5, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  expect_true(is.finite(fl$p_C))
  expect_true(is.finite(fl$L))
  # Zero reserve => zero mobilisation
  expect_equal(fl$p_C, 0)
})

test_that("deb_fluxes with zero volume does not crash", {
  fl <- deb_fluxes(E = 10, V = 0, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  expect_true(is.finite(fl$L))
  expect_equal(fl$L, 0)
  expect_equal(fl$p_A, 0)   # L^2 = 0
  expect_equal(fl$p_M, 0)   # p_M * V = 0
})

test_that("deb_fluxes with very small volume is stable", {
  fl <- deb_fluxes(E = 1e-10, V = 1e-15, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  expect_true(all(is.finite(unlist(fl))))
})

test_that("deb_fluxes with very large values is stable", {
  fl <- deb_fluxes(E = 1e8, V = 1e6, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  expect_true(all(is.finite(unlist(fl))))
})

test_that("deb_fluxes growth is non-negative", {
  fl <- deb_fluxes(E = 10, V = 0.5, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  expect_gte(fl$p_G, 0)
})

test_that("deb_fluxes reproduction is non-negative", {
  fl <- deb_fluxes(E = 50, V = 2.0, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400, k_J = 0.002, E_Hp = 50)
  expect_gte(fl$p_R, 0)
})

test_that("deb_fluxes f=0 means zero assimilation", {
  fl <- deb_fluxes(E = 10, V = 0.5, f = 0.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  expect_equal(fl$p_A, 0)
})

test_that("deb_fluxes p_J depends on k_J and E_Hp", {
  fl <- deb_fluxes(E = 10, V = 0.5, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400,
                   k_J = 0.01, E_Hp = 100)
  expect_equal(fl$p_J, 0.01 * 100)
})

test_that("deb_fluxes with k_J=0 and E_Hp=0 gives p_J=0 and p_R>=0", {
  fl <- deb_fluxes(E = 10, V = 0.5, f = 1.0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400, k_J = 0, E_Hp = 0)
  expect_equal(fl$p_J, 0)
  expect_gte(fl$p_R, 0)
})

test_that("deb_fluxes starvation: maintenance exceeds mobilisation", {
  # Very low E relative to V => kappa*p_C < p_M
  fl <- deb_fluxes(E = 0.001, V = 10, f = 0,
                   p_Am = 5, p_M = 0.5, kappa = 0.75,
                   v = 0.2, E_G = 400)
  # Growth should be zero (clamped), not negative
  expect_equal(fl$p_G, 0)
})


# --- assert_positive ---

test_that("assert_positive accepts valid input", {
  expect_silent(assert_positive(1.0, "test"))
  expect_silent(assert_positive(0.001, "test"))
  expect_silent(assert_positive(1e10, "test"))
})

test_that("assert_positive rejects zero", {
  expect_error(assert_positive(0, "x"), "positive")
})

test_that("assert_positive rejects negative", {
  expect_error(assert_positive(-1, "x"), "positive")
})

test_that("assert_positive rejects non-numeric", {
  expect_error(assert_positive("a", "x"), "positive")
})

test_that("assert_positive rejects vector", {
  expect_error(assert_positive(c(1, 2), "x"), "positive")
})

test_that("assert_positive rejects NA", {
  expect_error(assert_positive(NA_real_, "x"), "positive")
})

test_that("assert_positive rejects NULL", {
  expect_error(assert_positive(NULL, "x"), "positive")
})

test_that("assert_positive rejects Inf", {
  expect_error(assert_positive(Inf, "x"), "positive")
})
