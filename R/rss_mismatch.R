# RSS R-reference mismatch handling.
#
# Single home for code that targets the discrepancy between the
# supplied R reference and the target population:
#
#   * 1-D MAP optimizer for the variance component lambda_bias
#     (estimate_lambda_bias)
#   * per-variable inflation factor used inside the SER step
#     (compute_shat2_inflation, compute_shat2_inflation_rss)
#   * model-state storage helper for per-slot inflation diagnostics
#     (apply_inflation_state)
#   * per-sweep region-level fit (fit_R_bias) -- the new piece
#   * residual R-mismatch QC diagnostic Q_art (later commit)
#
# Storage convention on the model:
#   model$lambda_bias    SS / ss_mixture: scalar (set once per sweep
#                        by fit_R_bias).
#                        rss_lambda:      length-L vector (legacy
#                        per-slot path).
#   model$B_corrected    same shape as model$lambda_bias.
#   model$shat2_inflation per-variable inflation vector of length p,
#                        consumed by the SER step.
#
# The dispatch files (sufficient_stats_methods.R, ss_mixture_methods.R,
# rss_lambda_methods.R, iterative_bayesian_stepwise_selection.R) only
# call into this module; the bodies live here.

# =============================================================================
# 1-D MAP OPTIMIZER FOR lambda_bias
# =============================================================================

# Estimate extra R-bias variance beyond finite-reference uncertainty.
# Likelihood on the z-score residual scale,
#   tau_j^2 = sigma2 + (1/finite_R_B + lambda_bias) * s_j,
# with a half-Cauchy(prior_scale) prior on u = sqrt(lambda_bias).
# The Fisher-information boundary SE,
#   SE_0 = sqrt(2) * sigma2 / sqrt(sum(s^2)),
# defines a data-driven floor: estimates below 0.1 * SE_0 are zeroed.
# This both suppresses Brent boundary noise and replaces ad-hoc display
# thresholds with one rule; "none" short-circuits before optimization.
#' @keywords internal
estimate_lambda_bias <- function(r, s, sigma2, finite_R_B, method) {
  if (is.null(method) || method == "none")
    return(0)
  keep <- is.finite(r) & is.finite(s) & s > .Machine$double.eps
  if (!any(keep) || !is.finite(sigma2) || sigma2 <= .Machine$double.eps)
    return(0)

  cache <- list(r2 = r[keep]^2, s = s[keep])
  cache$base <- sigma2 + cache$s / finite_R_B
  pos <- (cache$r2 - cache$base) / cache$s
  pos <- pos[is.finite(pos) & pos > 0]
  prior_scale <- sqrt(max(1 / finite_R_B, 1 / 10000))
  upper_lambda <- max(c(1, 100 / finite_R_B, 100 * prior_scale^2,
                        10 * pos), na.rm = TRUE)
  upper_u <- sqrt(upper_lambda)

  nll <- function(u) {
    lambda_bias <- u^2
    tau <- cache$base + lambda_bias * cache$s
    0.5 * sum(log(tau) + cache$r2 / tau) + log1p((u / prior_scale)^2)
  }
  lambda_hat <- optimize(nll, interval = c(0, upper_u))$minimum^2

  ss2 <- sum(cache$s^2)
  if (ss2 <= 0) return(0)
  se_boundary <- sqrt(2) * sigma2 / sqrt(ss2)
  if (lambda_hat < 0.1 * se_boundary) 0 else lambda_hat
}

# =============================================================================
# PER-VARIABLE INFLATION
# =============================================================================

# SS-path per-variable inflation factor tau_j^2 / sigma2 with
#   tau_j^2 = sigma2 + (1/finite_R_B + lambda_bias) * (eta_j^2 + v_g),
#   eta_j^2 = XtXr_without_l[j]^2 / (n-1)   (z-score scale)
#   v_g     = sum(b_minus_l * XtXr_without_l).
# Reads the region-level scalar lambda_bias from model (set once per
# sweep by fit_R_bias). Per-slot lambda_bias re-fitting was removed:
# the previous design re-estimated lambda_bias inside every SER step
# from the leave-one-effect residual, which intentionally contains the
# lth sparse signal and so confounded signal with R-bias. The fix is
# the per-sweep fit_R_bias hook in ibss_fit; this function only
# applies the scalar to the slot-specific xi_l.
# Returns NULL when no inflation applies, otherwise a list with the
# inflation vector and lambda_bias / B_corrected = NULL so that
# apply_inflation_state does not write per-slot diagnostics on the SS
# path (those are scalars on the model, set by fit_R_bias).
#' @keywords internal
compute_shat2_inflation <- function(data, model, XtXr_without_l, b_minus_l, r) {
  finite_R_B <- data$finite_R_B
  if (is.null(finite_R_B) ||
      model$sigma2 <= .Machine$double.eps) {
    return(NULL)
  }
  v_g     <- max(sum(b_minus_l * XtXr_without_l), 0)
  eta2    <- XtXr_without_l^2 / (data$n - 1)
  s <- eta2 + v_g
  lambda_bias <- if (!is.null(model$lambda_bias)) model$lambda_bias[1] else 0
  infl <- 1 + (1 / finite_R_B + lambda_bias) * s / model$sigma2
  list(infl = infl, lambda_bias = NULL, B_corrected = NULL)
}

# rss_lambda counterpart of compute_shat2_inflation: z-score scale,
#   eta_j^2 = Rz_without_l[j]^2 (no n-1 division on z-scale),
#   v_g = max(b_minus_l' Rz_without_l, 0).
# In multi-panel the model-level finite_R_B is the omega-weighted
# (sum_k omega_k^2 / B_k)^{-1}.
# This function preserves the legacy per-slot lambda_bias re-fit (using
# the post-sweep fitted residual r_full = z - Rz_full) because the
# rss_lambda dispatch is out of scope of this redesign and will be
# retired by a separate change.
#' @keywords internal
compute_shat2_inflation_rss <- function(data, model, Rz_without_l, b_minus_l) {
  # Use model-level finite_R_B (updated by omega) if available, else data-level.
  finite_R_B <- if (!is.null(model$finite_R_B)) model$finite_R_B else data$finite_R_B
  if (is.null(finite_R_B) || model$sigma2 <= .Machine$double.eps) return(NULL)
  v_g  <- max(sum(b_minus_l * Rz_without_l), 0)
  eta2 <- Rz_without_l^2   # z-score scale: no (n-1) division needed
  s <- eta2 + v_g
  R_bias <- if (!is.null(data$R_bias)) data$R_bias else "none"
  if (R_bias == "none") {
    lambda_bias <- 0
  } else {
    # Generative-model target: lambda_bias is local R-mismatch variance
    # estimated from the residual after all currently active effects
    # have been explained. Do not estimate it from the leave-one-effect
    # residual, which intentionally contains the lth sparse signal
    # during the SER update.
    b_full <- if (!is.null(model$zbar)) model$zbar else {
      sw <- if (!is.null(model$slot_weights)) model$slot_weights else
              rep(1, nrow(model$alpha))
      colSums(sw * model$alpha * model$mu)
    }
    Rz_full <- if (!is.null(model$Rz))
                 model$Rz
               else as.vector(compute_Rv(data, b_full, model$X_meta))
    r_full <- data$z - Rz_full
    v_g_full <- max(sum(b_full * Rz_full), 0)
    s_full <- Rz_full^2 + v_g_full
    lambda_bias <- estimate_lambda_bias(r_full, s_full, model$sigma2,
                                        finite_R_B, R_bias)
  }
  infl <- 1 + (1 / finite_R_B + lambda_bias) * s / model$sigma2
  if (R_bias == "none") {
    list(infl = infl, lambda_bias = NULL, B_corrected = NULL)
  } else {
    list(infl = infl,
         lambda_bias = lambda_bias,
         B_corrected    = 1 / (1 / finite_R_B + lambda_bias))
  }
}

# =============================================================================
# MODEL-STATE STORAGE FOR PER-SLOT INFLATION DIAGNOSTICS
# =============================================================================

# Unpack the inflation list from compute_shat2_inflation* into the model:
# bare per-variant vector at model$shat2_inflation, plus per-slot
# diagnostics model$lambda_bias[l] and model$B_corrected[l] when the
# caller (rss_lambda dispatch) provides them. SS-path callers pass
# infl_state$lambda_bias = NULL because the SS path stores lambda_bias
# as a scalar set once per sweep by fit_R_bias.
#' @keywords internal
apply_inflation_state <- function(model, infl_state, l) {
  if (is.null(infl_state)) {
    model$shat2_inflation <- NULL
    return(model)
  }
  model$shat2_inflation <- infl_state$infl
  L <- nrow(model$alpha)
  if (!is.null(infl_state$lambda_bias)) {
    if (is.null(model$lambda_bias) || length(model$lambda_bias) != L)
      model$lambda_bias <- rep(0, L)
    model$lambda_bias[l] <- infl_state$lambda_bias
  }
  if (!is.null(infl_state$B_corrected)) {
    if (is.null(model$B_corrected) || length(model$B_corrected) != L)
      model$B_corrected <- rep(NA_real_, L)
    model$B_corrected[l] <- infl_state$B_corrected
  }
  model
}

# =============================================================================
# PER-SWEEP REGION-LEVEL fit_R_bias
# =============================================================================

#' Fit the region-level lambda_bias from the post-sweep fitted residual.
#'
#' Math (see archive/ld_mismatch_generativemodel.tex):
#'   beta_bar    = colSums(slot_weight * alpha * mu)        (full posterior mean, betahat scale)
#'   XtXr_full   = X'X * beta_bar = (n-1) * R * beta_bar
#'   r_fit       = (data$Xty - XtXr_full) / sqrt(n-1)       (z-scale fitted residual)
#'   eta_fit_j^2 = XtXr_full[j]^2 / (n-1)                   (z-scale per-variant signal)
#'   v_g_fit     = sum(beta_bar * XtXr_full)                (= beta_bar_z' R beta_bar_z)
#'   xi_fit_j    = eta_fit_j^2 + v_g_fit
#' MAP for lambda_bias on the working likelihood
#'   r_fit_j ~ N(0, sigma2 + (1/B + lambda) * xi_fit_j)
#' with half-Cauchy(scale = sqrt(max(1/B, 1e-4))) prior on sqrt(lambda).
#' Fisher SE zero-mask applied (see estimate_lambda_bias).
#'
#' Replaces the per-slot re-fit that used to live inside
#' compute_shat2_inflation, which estimated lambda_bias from the
#' intra-sweep r_full_z that drifts as the slot loop progresses.
#'
#' @keywords internal
#' @noRd
fit_R_bias <- function(data, params, model) {
  R_bias <- if (!is.null(params$R_bias)) params$R_bias else "none"
  if (R_bias == "none") return(model)
  finite_R_B <- data$finite_R_B
  if (is.null(finite_R_B) || !is.finite(model$sigma2) ||
      model$sigma2 <= .Machine$double.eps)
    return(model)
  if (!inherits(data, c("ss", "ss_mixture"))) return(model)

  sw <- if (!is.null(model$slot_weights)) model$slot_weights else
          rep(1, nrow(model$alpha))
  b_full    <- colSums(sw * model$alpha * model$mu)
  XtXr_full <- if (!is.null(model$XtXr))
                 model$XtXr
               else compute_Rv(data, b_full)
  nm1 <- if (!is.null(data$nm1)) data$nm1 else (data$n - 1)
  if (!is.finite(nm1) || nm1 <= 0) return(model)

  r_fit_z  <- (data$Xty - XtXr_full) / sqrt(nm1)
  v_g_full <- max(sum(b_full * XtXr_full), 0)
  s_full   <- XtXr_full^2 / nm1 + v_g_full

  model$lambda_bias <- estimate_lambda_bias(r_fit_z, s_full, model$sigma2,
                                            finite_R_B, R_bias)
  model$B_corrected <- 1 / (1 / finite_R_B + model$lambda_bias)
  model
}
