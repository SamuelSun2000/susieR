# Greedy-L outer loop in susie_workhorse.
# Contracts: (1) L_greedy = NULL is bit-identical to fixed-L susie.
# (2) L_greedy != NULL grows L until min(lbf) < lbf_min or L reaches
# params$L.

set.seed(42)
N <- 200
J <- 100
X <- matrix(rnorm(N * J), N, J)
true_idx <- c(10, 30)                   # K = 2 real effects
beta <- numeric(J)
beta[true_idx] <- c(2.5, -1.8)
y <- X %*% beta + rnorm(N, sd = 0.3)

test_that("L_greedy = NULL is bit-identical to fixed-L susie", {
  fit_fixed <- susie(X, y, L = 5)

  obj <- susie(X, y, L = 5, init_only = TRUE)
  obj$params$L_greedy <- NULL
  fit_direct <- susie_workhorse(obj$data, obj$params)

  expect_equal(fit_direct$alpha, fit_fixed$alpha, tolerance = 0)
  expect_equal(fit_direct$lbf,   fit_fixed$lbf,   tolerance = 0)
  expect_equal(fit_direct$elbo,  fit_fixed$elbo,  tolerance = 0)
})

test_that("L_greedy grows L in steps of L_greedy, capped at params$L", {
  obj <- susie(X, y, L = 10, init_only = TRUE)
  obj$params$L_greedy <- 3
  obj$params$lbf_min  <- 0.1
  fit <- susie_workhorse(obj$data, obj$params)

  final_L <- nrow(fit$alpha)
  expect_true(final_L %in% c(3, 6, 9, 10))     # multiple of 3, capped at 10
  expect_lte(final_L, 9)                        # K = 2 real, saturates early
})

test_that("L_greedy >= K_true saturates in a single round", {
  # L_greedy = 6, K = 2 real. Round 1 at L = 6 has empty slots so
  # min(lbf) < lbf_min fires immediately, no wasted second round.
  obj <- susie(X, y, L = 12, init_only = TRUE)
  obj$params$L_greedy <- 6
  obj$params$lbf_min  <- 0.1
  fit <- susie_workhorse(obj$data, obj$params)

  expect_identical(nrow(fit$alpha), 6L)
})

test_that("K_true > L_greedy keeps growing past the first round", {
  set.seed(7)
  Xh <- matrix(rnorm(N * J), N, J)
  bh <- numeric(J)
  bh[c(5, 20, 45, 70)] <- c(2.5, -2.0, 1.8, -1.5)   # K = 4 real effects
  yh <- Xh %*% bh + rnorm(N, sd = 0.3)

  obj <- susie(Xh, yh, L = 10, init_only = TRUE)
  obj$params$L_greedy <- 3
  obj$params$lbf_min  <- 0.1
  fit <- susie_workhorse(obj$data, obj$params)

  expect_gte(nrow(fit$alpha), 6)
  expect_lte(nrow(fit$alpha), 10)
})

test_that("L_greedy = L stops after one round at L", {
  obj <- susie(X, y, L = 3, init_only = TRUE)
  obj$params$L_greedy <- 3
  fit <- susie_workhorse(obj$data, obj$params)

  expect_identical(nrow(fit$alpha), 3L)
})
