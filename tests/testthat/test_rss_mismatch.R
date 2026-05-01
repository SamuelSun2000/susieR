context("RSS R-reference mismatch (R_bias correction)")

# ---- API surface guards ----

test_that("R_bias = 'map_qc' runs and returns Q_art diagnostics", {
  set.seed(11)
  p <- 20
  n <- 1000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
                   R_bias = "map_qc", max_iter = 2, verbose = FALSE)
  d <- fit$finite_R_diagnostics
  expect_true(!is.null(d$Q_art))
  expect_true(d$Q_art >= 0 && d$Q_art <= 1)
  expect_true(is.logical(d$artifact_flag))
  expect_true(d$mode_label %in% c("normal", "warning", "conservative"))
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

# ---- Q_art unit tests ----

test_that("compute_Q_art recovers Q ~ 1 when r_fit lies in low-eigen direction", {
  # Diagonal R with eigenvalues (2, 1, 1e-6). Default eig_delta_rel=1e-3
  # selects only the third eigenvalue.
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  r_fit <- c(0, 0, 1)  # purely in the low-eigen direction
  out <- susieR:::compute_Q_art(eig, r_fit)
  expect_equal(out$Q_art, 1, tolerance = 1e-12)
  expect_true(out$evaluable)
  expect_equal(out$low_eigen_count, 1L)
})

test_that("compute_Q_art returns Q ~ 0 when r_fit avoids low-eigen directions", {
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  r_fit <- c(1, 0.5, 0)  # fully in top two eigen directions
  out <- susieR:::compute_Q_art(eig, r_fit)
  expect_equal(out$Q_art, 0, tolerance = 1e-12)
})

test_that("compute_Q_art is non-evaluable when r_fit has negligible energy", {
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  out <- susieR:::compute_Q_art(eig, rep(0, 3))
  expect_equal(out$Q_art, 0)
  expect_false(out$evaluable)
})

test_that("compute_Q_art is non-evaluable when no low-eigenvalues exist", {
  V <- diag(3)
  d <- c(2, 1, 0.5)  # all > 1e-3 * 2 = 2e-3
  eig <- list(values = d, vectors = V)
  out <- susieR:::compute_Q_art(eig, c(1, 0, 0))
  expect_equal(out$low_eigen_count, 0L)
  expect_false(out$evaluable)
})

test_that("compute_Q_art is in [0, 1] for typical inputs", {
  V <- diag(3)
  d <- c(2, 1, 1e-6)
  eig <- list(values = d, vectors = V)
  for (r_fit in list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1),
                     c(0.5, 0.5, 0.5), c(-1, 1, -1))) {
    out <- susieR:::compute_Q_art(eig, r_fit)
    expect_true(out$Q_art >= 0 && out$Q_art <= 1)
  }
})

# ---- map_qc end-to-end smoke ----

test_that("map_qc on well-behaved data yields Q_art near 0 and no flag", {
  set.seed(11)
  p <- 25
  n <- 2000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)
  z[3] <- 4

  fit <- susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
                   R_bias = "map_qc", max_iter = 5, verbose = FALSE)
  d <- fit$finite_R_diagnostics
  expect_lt(d$Q_art, 0.1)
  expect_false(d$artifact_flag)
  expect_equal(d$mode_label, "normal")
})

test_that("map_qc emits a true R warning when artifact_flag triggers", {
  rho <- 0.9999
  z <- c(-8, -8)
  R <- matrix(c(1, -rho, -rho, 1), 2, 2)
  expect_warning(
    fit <- susie_rss(z = z, R = R, n = 5000, L = 1, finite_R = 1e6,
                     R_bias = "map_qc", max_iter = 5,
                     prior_variance = 50, estimate_prior_variance = FALSE,
                     estimate_residual_variance = FALSE, verbose = FALSE,
                     check_R = FALSE),
    "Residual R-bias artifact detected"
  )
  expect_true(fit$finite_R_diagnostics$artifact_flag)
  expect_equal(fit$finite_R_diagnostics$Q_art, 1, tolerance = 1e-6)
})

test_that("map_qc surfaces Q_art and mode_label diagnostics", {
  set.seed(12)
  p <- 25; n <- 2000
  X <- matrix(rnorm(n * p), n, p)
  R <- cor(X)
  z <- rnorm(p)

  fit <- susie_rss(z = z, R = R, n = n, L = 3, finite_R = 5000,
                   R_bias = "map_qc", max_iter = 3, verbose = FALSE)
  d <- fit$finite_R_diagnostics
  for (fld in c("Q_art", "artifact_flag", "artifact_evaluable",
                "low_eigen_count", "low_eigen_fraction", "eig_delta",
                "mode_label", "lambda_bias", "B_corrected"))
    expect_true(!is.null(d[[fld]]), info = paste("missing diagnostic:", fld))
})

test_that("map_qc with X-input runs and surfaces Q_art", {
  set.seed(15)
  p <- 25; n <- 2000
  X <- matrix(rnorm(n * p), n, p)
  X <- scale(X, center = TRUE, scale = TRUE)
  beta <- rep(0, p); beta[5] <- 0.4
  y <- drop(X %*% beta + rnorm(n))
  z <- as.numeric(crossprod(X, y) / sqrt(diag(crossprod(X))))
  fit <- susie_rss(z = z, X = X, n = n, L = 3, finite_R = 5000,
                   R_bias = "map_qc", max_iter = 3, verbose = FALSE)
  d <- fit$finite_R_diagnostics
  expect_true(!is.null(d$Q_art))
  expect_true(d$Q_art >= 0 && d$Q_art <= 1)
})

test_that("map_qc works on lambda=0 multi-panel SS path", {
  set.seed(19)
  n <- 80
  p <- 12
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p), n, p)
  z <- rnorm(p)

  fit <- susie_rss(z = z, X = list(X1, X2), n = 1000, L = 3,
                   finite_R = TRUE, R_bias = "map_qc", max_iter = 3,
                   verbose = FALSE)
  d <- fit$finite_R_diagnostics
  expect_true(!is.null(d$Q_art))
  expect_true(d$Q_art >= 0 && d$Q_art <= 1)
  expect_length(d$lambda_bias, 1)
})
