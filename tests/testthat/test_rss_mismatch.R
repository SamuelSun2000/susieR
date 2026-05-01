context("RSS R-reference mismatch (R_bias correction)")

# ---- API surface guards ----

test_that("R_bias = 'map_qc' is accepted by match.arg but errors as NYI", {
  set.seed(11)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
              R_bias = "map_qc", max_iter = 2, verbose = FALSE),
    "not yet implemented"
  )
})

test_that("R_bias = 'map_qc' with lambda > 0 errors at the rss_lambda constructor", {
  set.seed(13)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  # lambda > 0 routes through rss_lambda_constructor; "map_qc" is rejected
  # there with a clear message before reaching the susie.R NYI guard.
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
              lambda = 0.1, R_bias = "map_qc",
              max_iter = 2, verbose = FALSE),
    "map_qc"
  )
})

test_that("Optional artifact args validate ranges", {
  set.seed(17)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
              R_bias = "map", artifact_threshold = -0.1,
              max_iter = 2, verbose = FALSE),
    "artifact_threshold"
  )
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
              R_bias = "map", artifact_threshold = 1.1,
              max_iter = 2, verbose = FALSE),
    "artifact_threshold"
  )
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
              R_bias = "map", eig_delta_rel = -1,
              max_iter = 2, verbose = FALSE),
    "eig_delta_rel"
  )
  expect_error(
    susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
              R_bias = "map", B_artifact = 0,
              max_iter = 2, verbose = FALSE),
    "B_artifact"
  )
})

# ---- Region-level scalar lambda_bias on the SS path ----

test_that("SS path stores scalar lambda_bias / B_corrected (not per-slot)", {
  set.seed(101)
  p <- 25
  n <- 1500
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
                   R_bias = "map", max_iter = 5, verbose = FALSE)
  expect_length(fit$finite_R_diagnostics$lambda_bias, 1)
  expect_length(fit$finite_R_diagnostics$B_corrected, 1)
  expect_true(fit$finite_R_diagnostics$lambda_bias >= 0)
  expect_equal(fit$finite_R_diagnostics$B_corrected,
               1 / (1 / fit$finite_R_diagnostics$B +
                      fit$finite_R_diagnostics$lambda_bias),
               tolerance = 1e-12)
})

test_that("rss_lambda path keeps per-slot lambda_bias (out of scope)", {
  set.seed(107)
  p <- 25
  n <- 1500
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
                   lambda = 0.1, R_bias = "map",
                   max_iter = 5, verbose = FALSE)
  expect_length(fit$finite_R_diagnostics$lambda_bias, 3)
  expect_length(fit$finite_R_diagnostics$B_corrected, 3)
})
