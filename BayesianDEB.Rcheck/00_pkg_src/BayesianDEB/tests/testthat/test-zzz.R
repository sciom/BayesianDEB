# ===========================================================
# Tests: zzz.R  (stan_file, package hooks)
# ===========================================================

test_that("stan_file returns valid path for all model names", {
  models <- c("bdeb_individual_growth", "bdeb_growth_repro",
              "bdeb_hierarchical_growth", "bdeb_debtox")

  for (m in models) {
    path <- stan_file(m)
    expect_true(file.exists(path), info = paste("model:", m))
    expect_true(grepl("\\.stan$", path), info = paste("model:", m))
  }
})

test_that("stan_file rejects nonexistent model", {
  expect_error(stan_file("nonexistent_model"), "not found")
})

test_that("stan_file rejects empty string", {
  expect_error(stan_file(""), "not found")
})

test_that("bundled Stan files contain expected blocks", {
  models <- c("bdeb_individual_growth", "bdeb_growth_repro",
              "bdeb_hierarchical_growth", "bdeb_debtox")

  for (m in models) {
    path <- stan_file(m)
    content <- readLines(path)
    text <- paste(content, collapse = "\n")

    expect_true(grepl("functions\\s*\\{", text),
                info = paste(m, "missing functions block"))
    expect_true(grepl("data\\s*\\{", text),
                info = paste(m, "missing data block"))
    expect_true(grepl("parameters\\s*\\{", text),
                info = paste(m, "missing parameters block"))
    expect_true(grepl("model\\s*\\{", text),
                info = paste(m, "missing model block"))
  }
})

test_that("individual and hierarchical Stan models contain ode_bdf", {
  for (m in c("bdeb_individual_growth", "bdeb_hierarchical_growth")) {
    path <- stan_file(m)
    content <- paste(readLines(path), collapse = "\n")
    expect_true(grepl("ode_bdf", content), info = m)
  }
})

test_that("debtox Stan model contains EC50 in generated quantities", {
  path <- stan_file("bdeb_debtox")
  content <- paste(readLines(path), collapse = "\n")
  expect_true(grepl("EC50", content))
  expect_true(grepl("NEC", content))
  expect_true(grepl("generated quantities", content))
})

test_that("growth_repro Stan model contains neg_binomial", {
  path <- stan_file("bdeb_growth_repro")
  content <- paste(readLines(path), collapse = "\n")
  expect_true(grepl("neg_binomial_2", content))
})

test_that("hierarchical Stan model contains non-centred parameterisation", {
  path <- stan_file("bdeb_hierarchical_growth")
  content <- paste(readLines(path), collapse = "\n")
  expect_true(grepl("z_log_p_Am", content))
  expect_true(grepl("mu_log_p_Am", content))
  expect_true(grepl("sigma_log_p_Am", content))
})
