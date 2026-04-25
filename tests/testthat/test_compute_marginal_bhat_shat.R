# compute_marginal_bhat_shat
#
# Per-position marginal OLS regression helper. Used by susieR's own
# T = 1 SER path (cosmetic refactor candidate), mvsusieR's OLS
# branch, and mfsusieR's prior init / SER. Three contracts:
#   1. Vector-Y input is treated as a single-column matrix.
#   2. predictor_weights override matches recompute from colSums(X^2).
#   3. sigma2 supplied gives single-effect-residual Shat shape.

set.seed(42)
N <- 50
J <- 5

test_that("vector Y is treated as a single-column matrix", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  y <- rnorm(N)

  out_vec <- compute_marginal_bhat_shat(X, y)
  out_mat <- compute_marginal_bhat_shat(X, matrix(y, ncol = 1))

  expect_equal(dim(out_vec$Bhat), c(J, 1L))
  expect_equal(dim(out_vec$Shat), c(J, 1L))
  expect_equal(out_vec$Bhat, out_mat$Bhat, tolerance = 0)
  expect_equal(out_vec$Shat, out_mat$Shat, tolerance = 0)
})

test_that("matrix Y returns J x T Bhat / Shat", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- matrix(rnorm(N * 3), N, 3)

  out <- compute_marginal_bhat_shat(X, Y)

  expect_equal(dim(out$Bhat), c(J, 3L))
  expect_equal(dim(out$Shat), c(J, 3L))
})

test_that("predictor_weights override matches recompute from colSums(X^2)", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- matrix(rnorm(N * 4), N, 4)

  pw <- colSums(X^2)

  out_default <- compute_marginal_bhat_shat(X, Y)
  out_override <- compute_marginal_bhat_shat(X, Y, predictor_weights = pw)

  expect_equal(out_default$Bhat, out_override$Bhat, tolerance = 0)
  expect_equal(out_default$Shat, out_override$Shat, tolerance = 0)
})

test_that("sigma2 supplied gives single-effect-residual Shat (sqrt(sigma2 / pw))", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- matrix(rnorm(N * 2), N, 2)

  out <- compute_marginal_bhat_shat(X, Y, sigma2 = 0.5)

  pw <- colSums(X^2)
  expected_shat <- matrix(sqrt(0.5 / pw), nrow = J, ncol = 2)
  expect_equal(out$Shat, expected_shat, tolerance = 0)
})

test_that("Bhat = X'Y / colSums(X^2) for centred X", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- matrix(rnorm(N * 3), N, 3)

  out <- compute_marginal_bhat_shat(X, Y)
  expected_bhat <- crossprod(X, Y) / colSums(X^2)

  expect_equal(out$Bhat, expected_bhat, tolerance = 0)
})

test_that("Shat (no sigma2) matches per-column residual SD / sqrt(n-1)", {
  X <- matrix(rnorm(N * J), N, J)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- matrix(rnorm(N * 2), N, 2)

  out <- compute_marginal_bhat_shat(X, Y)

  # Manual recompute, no Rfast.
  Bhat <- crossprod(X, Y) / colSums(X^2)
  expected_shat <- matrix(0, nrow = J, ncol = 2)
  for (t in 1:2) {
    for (j in 1:J) {
      r <- Y[, t] - X[, j] * Bhat[j, t]
      expected_shat[j, t] <- sqrt(var(r))
    }
  }
  expected_shat <- expected_shat / sqrt(N - 1)

  expect_equal(out$Shat, expected_shat, tolerance = 1e-12)
})
