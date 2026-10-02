# =============================================================================
# Multi-panel two-stage fit (control$multi_panel_refit)
# =============================================================================

test_that("multi-panel R list: two-stage fit equals a single fit on sum_k omega_k R_k", {
  set.seed(20)
  p <- 15; Bn <- 80
  X1 <- matrix(rnorm(Bn * p), Bn, p)
  X2 <- matrix(rnorm(Bn * p), Bn, p)
  R1 <- cor(X1); R2 <- cor(X2)
  z <- rnorm(p); z[4] <- 5
  fit <- suppressWarnings(susie_rss(z = z, R = list(R1, R2), n = 1000,
                                    L = 3, max_iter = 20,
                                    R_finite = c(80, 120),
                                    R_mismatch = "eb_mix"))
  joint <- suppressWarnings(susie_rss(z = z, R = list(R1, R2), n = 1000,
                                      L = 3, max_iter = 20,
                                      R_finite = c(80, 120),
                                      R_mismatch = "eb_mix",
                                      control = list(multi_panel_refit = FALSE)))
  om <- fit$omega_weights
  expect_equal(om, joint$omega_weights)
  expect_equal(sum(om), 1, tolerance = 1e-12)
  expect_false(is.null(fit$single_panel_fits))

  direct <- suppressWarnings(susie_rss(z = z, R = om[1] * R1 + om[2] * R2,
                                       n = 1000, L = 3, max_iter = 20,
                                       R_finite = 1 / sum(om^2 / c(80, 120)),
                                       R_mismatch = "eb_mix"))
  expect_equal(fit$pip, direct$pip)
  expect_equal(fit$alpha, direct$alpha)
  expect_equal(fit$R_finite_diagnostics$B, 1 / sum(om^2 / c(80, 120)))
  expect_equal(fit$R_finite_diagnostics$lambda_bias,
               direct$R_finite_diagnostics$lambda_bias)
})

test_that("multi-panel X list: weighted sketch reproduces sum_k omega_k R_k", {
  set.seed(21)
  p <- 40; B1 <- 12; B2 <- 9
  X1 <- matrix(rnorm(B1 * p, mean = 2, sd = 3), B1, p)
  X2 <- matrix(rnorm(B2 * p), B2, p)
  om <- c(0.3, 0.7)
  ref <- form_weighted_reference(list(X1, X2), "X", om)
  expect_equal(nrow(ref$X), B1 + B2)
  expect_equal(crossprod(ref$X), om[1] * cor(X1) + om[2] * cor(X2))
  # A zero weight drops that panel's rows.
  expect_equal(nrow(form_weighted_reference(list(X1, X2), "X", c(0, 1))$X), B2)

  z <- rnorm(p); z[7] <- 4.5
  fit <- suppressWarnings(susie_rss(z = z, X = list(X1, X2), n = 500, L = 2,
                                    max_iter = 20, R_finite = TRUE,
                                    R_mismatch = "eb"))
  om <- fit$omega_weights
  ref <- form_weighted_reference(list(X1, X2), "X", om)
  direct <- suppressWarnings(susie_rss(z = z, X = ref$X, n = 500, L = 2,
                                       max_iter = 20,
                                       R_finite = 1 / sum(om^2 / c(B1, B2)),
                                       R_mismatch = "eb"))
  expect_equal(fit$pip, direct$pip)
  expect_equal(fit$R_finite_diagnostics$B, 1 / sum(om^2 / c(B1, B2)))
})
