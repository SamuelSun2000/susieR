# RSS R-reference mismatch handling.
#
# Single home for code that targets the discrepancy between the
# supplied R reference and the target population. Active on the SS
# / ss_mixture dispatches; the rss_lambda dispatch (lambda > 0) does
# NOT use any of this (entry-level errors block lambda > 0 with
# R_finite or R_mismatch != "none").
#
#   * 1-D MAP optimizer for the variance component lambda_bias
#     (estimate_lambda_bias)
#   * per-variable inflation factor used inside the SER step
#     (compute_shat2_inflation)
#   * model-state storage helper for per-slot inflation diagnostics
#     (apply_inflation_state)
#   * SER-protected initialization for the recommended EB path
#     (initialize_R_mismatch)
#   * per-sweep region-level fit (fit_R_mismatch)
#   * residual R-mismatch QC diagnostic Q_art (always with R_mismatch)
#
# Storage convention on the model:
#   model$lambda_bias    scalar set once per sweep by fit_R_mismatch
#   model$B_corrected    1 / (1/B + lambda_bias)
#   model$shat2_inflation per-variable inflation vector of length p,
#                        consumed by the SER step.

# =============================================================================
# FINITE-REFERENCE SETUP AND DIAGNOSTICS
# =============================================================================

# Resolve R_finite into an explicit reference sample size B. FALSE is an
# explicit "off" setting; NULL means unspecified. R_finite = TRUE is only
# meaningful when the reference factor X is available; for precomputed R,
# the caller must provide B explicitly.
#' @keywords internal
resolve_R_finite <- function(R_finite, X = NULL, is_multi_panel = FALSE) {
  if (is.null(R_finite))
    return(NULL)
  if (identical(R_finite, FALSE))
    return(NULL)
  if (isTRUE(R_finite)) {
    if (is.null(X))
      stop("R_finite = TRUE requires X input. ",
           "When using a precomputed R matrix, provide a positive number ",
           "specifying the reference sample size B instead.")
    if (is_multi_panel)
      return(vapply(X, nrow, integer(1)))
    return(nrow(X))
  }
  if (!is.numeric(R_finite) || any(!is.finite(R_finite)) ||
      any(R_finite <= 0)) {
    stop("R_finite must be NULL, FALSE, TRUE, or positive numeric value(s).")
  }
  if (is_multi_panel) {
    K <- if (is.null(X)) length(R_finite) else length(X)
    if (length(R_finite) == 1)
      return(rep(as.numeric(R_finite), K))
    if (length(R_finite) == K)
      return(as.numeric(R_finite))
    stop("For multi-panel input, R_finite must be FALSE, TRUE, a single ",
         "positive number, or one positive number per panel.")
  }
  if (length(R_finite) == 1)
    return(as.numeric(R_finite))
  stop("R_finite must be NULL, FALSE, TRUE, or a single positive number.")
}

# Compute finite-reference R diagnostics (debiased Frobenius norm,
# effective rank, r/B ratio, per-variant diagonal deviation from 1).
# Used by both summary_stats_constructor and rss_lambda_constructor.
#
# @param X Factor matrix (B x p), or NULL.
# @param R Precomputed R matrix (p x p), or NULL.
# @param B Reference panel sample size.
# @param p Number of variants.
# @param x_is_standardized If TRUE, X has been standardized so X'X = R_hat
#   directly (no normalization). If FALSE, R_hat = X'X/B so the Frobenius
#   norm needs a /B^2 correction.
# @return List with B, p, R_frob_sq_debiased, effective_rank, r_over_B,
#   Rhat_diag_deviation.
#' @keywords internal
compute_R_finite_diagnostics <- function(X = NULL, R = NULL, B, p,
                                         x_is_standardized = FALSE) {
  if (!is.null(X)) {
    A <- tcrossprod(X)           # B x B Gram matrix
    R_frob_sq <- sum(A * A)      # ||XX'||_F^2 = ||X'X||_F^2
    if (!x_is_standardized)
      R_frob_sq <- R_frob_sq / nrow(X)^2
    Rhat_diag <- colSums(X^2)
    if (!x_is_standardized)
      Rhat_diag <- Rhat_diag / nrow(X)
  } else if (!is.null(R)) {
    R_frob_sq <- sum(R * R)
    Rhat_diag <- diag(R)
  } else {
    R_frob_sq <- p               # identity fallback
    Rhat_diag <- rep(1, p)
  }

  # Debiased Frobenius norm (Ledoit-Wolf unbiased estimator). In the
  # B = Inf limit there is no finite-reference debiasing term.
  R_frob_sq_db <- if (is.infinite(B)) R_frob_sq else
                    (B * R_frob_sq - p^2) / (B + 1)
  eff_rank <- p^2 / max(R_frob_sq_db, 1)

  list(
    B = B,
    p = p,
    R_frob_sq_debiased = R_frob_sq_db,
    effective_rank = eff_rank,
    r_over_B = eff_rank / B,
    Rhat_diag_deviation = abs(Rhat_diag - 1)
  )
}

# =============================================================================
# 1-D MAP OPTIMIZER FOR lambda_bias
# =============================================================================

# Estimate R-bias variance beyond any supplied finite-reference uncertainty.
# With R_finite_B = Inf this is the B^{-1} = 0 limit, so lambda_bias is the
# total continuous R-mismatch variance component.
# Likelihood on the z-score residual scale,
#   tau_j^2 = sigma2 + (1/R_finite_B + lambda_bias) * s_j,
# with a half-Cauchy(prior_scale) prior on u = sqrt(lambda_bias).
# The Fisher-information boundary SE,
#   SE_0 = sqrt(2) * sigma2 / sqrt(sum(s^2)),
# defines a data-driven floor: estimates below 0.1 * SE_0 are zeroed.
# This both suppresses Brent boundary noise and replaces ad-hoc display
# thresholds with one rule; "none" short-circuits before optimization.
#' @keywords internal
estimate_lambda_bias <- function(r, s, sigma2, R_finite_B, method) {
  if (is.null(method) || method == "none")
    return(0)
  keep <- is.finite(r) & is.finite(s) & s > .Machine$double.eps
  if (!any(keep) || !is.finite(sigma2) || sigma2 <= .Machine$double.eps)
    return(0)

  cache <- list(r2 = r[keep]^2, s = s[keep])
  cache$base <- sigma2 + cache$s / R_finite_B
  pos <- (cache$r2 - cache$base) / cache$s
  pos <- pos[is.finite(pos) & pos > 0]
  prior_scale <- sqrt(max(1 / R_finite_B, 1 / 10000))
  upper_lambda <- max(c(1, 100 / R_finite_B, 100 * prior_scale^2,
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
#   tau_j^2 = sigma2 + (1/R_finite_B + lambda_bias) * (eta_j^2 + v_g),
#   eta_j^2 = XtXr_without_l[j]^2 / (n-1)   (z-score scale)
#   v_g     = sum(b_minus_l * XtXr_without_l).
# Reads the region-level scalar lambda_bias from model (set once per
# sweep by fit_R_mismatch). Per-slot lambda_bias re-fitting was removed:
# the previous design re-estimated lambda_bias inside every SER step
# from the leave-one-effect residual, which intentionally contains the
# lth sparse signal and so confounded signal with R-bias. The fix is
# the per-sweep fit_R_mismatch hook in ibss_fit; this function only
# applies the scalar to the slot-specific xi_l.
# Returns NULL when no inflation applies, otherwise a list with the
# inflation vector and lambda_bias / B_corrected = NULL so that
# apply_inflation_state does not write per-slot diagnostics on the SS
# path (those are scalars on the model, set by fit_R_mismatch).
#' @keywords internal
compute_shat2_inflation <- function(data, model, XtXr_without_l, b_minus_l, r) {
  R_finite_B <- if (!is.null(model$R_finite_B)) model$R_finite_B else data$R_finite_B
  if (is.null(R_finite_B) ||
      model$sigma2 <= .Machine$double.eps) {
    return(NULL)
  }
  v_g     <- max(sum(b_minus_l * XtXr_without_l), 0)
  eta2    <- XtXr_without_l^2 / (data$n - 1)
  s <- eta2 + v_g
  lambda_bias <- if (is.null(model$lambda_bias)) 0 else model$lambda_bias
  infl <- 1 + (1 / R_finite_B + lambda_bias) * s / model$sigma2
  list(infl = infl, lambda_bias = NULL, B_corrected = NULL)
}

# =============================================================================
# MODEL-STATE STORAGE FOR PER-SLOT INFLATION DIAGNOSTICS
# =============================================================================

# Unpack the inflation list from compute_shat2_inflation into the model.
# Sets model$shat2_inflation to the per-variant inflation vector. The
# per-slot writes to model$lambda_bias[l] / model$B_corrected[l] gated
# below are dormant in the current code: SS / ss_mixture callers always
# pass infl_state$lambda_bias = NULL (the scalar lambda_bias is set
# once per sweep by fit_R_mismatch), and the rss_lambda path no longer
# calls this function. The per-slot machinery is retained as inert
# back-compat shim and will be removed when the next constructor pass
# converges on a single storage shape.
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
# PER-SWEEP REGION-LEVEL fit_R_mismatch
# =============================================================================

#' Fit the region-level lambda_bias from a fitted residual.
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
#' The same lambda_bias fit is followed by the Q_art residual artifact
#' diagnostic; see compute_Q_art. The artifact diagnostic emits an R
#' warning when flagged; it does not change lambda_bias or the SER
#' likelihood.
#'
#' @keywords internal
#' @noRd
compute_R_mismatch_state <- function(data, params, model, phase = "sweep") {
  R_mismatch <- if (!is.null(params$R_mismatch)) params$R_mismatch else "none"
  if (R_mismatch == "none") return(model)
  R_finite_B <- if (!is.null(model$R_finite_B)) model$R_finite_B else data$R_finite_B
  if (is.null(R_finite_B) || !is.finite(model$sigma2) ||
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
                                            R_finite_B, R_mismatch)
  model$B_corrected <- 1 / (1 / R_finite_B + model$lambda_bias)

  eigen_R <- get_R_mismatch_eigen(data, model)
  if (is.null(eigen_R))
    stop("R_mismatch requires data$eigen_R; ",
         "summary_stats_constructor should have cached it.")
  eig_delta_rel <- if (!is.null(params$eig_delta_rel))
                     params$eig_delta_rel else 1e-3
  eig_delta_abs <- if (!is.null(params$eig_delta_abs))
                     params$eig_delta_abs else 0
  art <- compute_Q_art(eigen_R, r_fit_z, eig_delta_rel, eig_delta_abs)
  threshold <- if (!is.null(params$artifact_threshold))
                 params$artifact_threshold else 0.1
  flagged <- isTRUE(art$evaluable) && isTRUE(art$Q_art > threshold)

  model$Q_art              <- art$Q_art
  model$artifact_evaluable <- art$evaluable
  model$artifact_flag      <- flagged
  model$low_eigen_count    <- art$low_eigen_count
  model$low_eigen_fraction <- art$low_eigen_count /
                              length(eigen_R$values)
  model$eig_delta          <- art$eig_delta

  if (flagged) {
    msg <- paste0("Residual R-bias artifact detected (Q_art = ",
                  sprintf("%.3g", art$Q_art),
                  " > threshold ", sprintf("%.3g", threshold),
                  "). Fine-mapping results may be unreliable with ",
                  "this R reference. Consider allele/QC review, ",
                  "multi-reference analysis, or conservative fallback.")
    model$mode_label <- "warning"
    warning_message(msg)
    warning(msg, call. = FALSE)
  } else {
    model$mode_label <- "normal"
  }

  if (isTRUE(params$track_fit)) {
    keep <- is.finite(r_fit_z) & is.finite(s_full) &
            s_full > .Machine$double.eps
    r2 <- r_fit_z[keep]^2
    s_keep <- s_full[keep]
    trace_row <- list(
      sweep = if (is.null(model$R_mismatch_trace)) 1L else
                length(model$R_mismatch_trace) + 1L,
      phase = phase,
      R_mismatch = R_mismatch,
      lambda_bias = model$lambda_bias,
      B_corrected = model$B_corrected,
      B = R_finite_B,
      sigma2 = model$sigma2,
      n_effects = nrow(model$alpha),
      n_nonzero_lbf = if (!is.null(model$lbf))
                        sum(is.finite(model$lbf) & model$lbf > 0)
                      else NA_integer_,
      mean_r2 = if (length(r2) > 0) mean(r2) else NA_real_,
      median_r2 = if (length(r2) > 0) stats::median(r2) else NA_real_,
      max_r2 = if (length(r2) > 0) max(r2) else NA_real_,
      mean_s = if (length(s_keep) > 0) mean(s_keep) else NA_real_,
      median_s = if (length(s_keep) > 0) stats::median(s_keep) else NA_real_,
      max_s = if (length(s_keep) > 0) max(s_keep) else NA_real_,
      cor_r2_s = if (length(r2) > 1 && stats::sd(r2) > 0 &&
                     stats::sd(s_keep) > 0)
                   suppressWarnings(stats::cor(r2, s_keep,
                                               method = "spearman"))
                 else NA_real_,
      Q_art = if (!is.null(model$Q_art)) model$Q_art else NA_real_,
      artifact_flag = if (!is.null(model$artifact_flag))
                        model$artifact_flag else NA,
      low_eigen_count = if (!is.null(model$low_eigen_count))
                          model$low_eigen_count else NA_integer_
    )
    model$R_mismatch_trace[[trace_row$sweep]] <- trace_row
  }

  model
}

#' SER-protected initialization for the R_mismatch EB path.
#'
#' The joint EB/sparse objective is path dependent. Starting lambda_bias at
#' zero can let secondary R-mismatch patterns enter as sparse effects before
#' the variance component is estimated. For R_mismatch = "eb", initialize only
#' in the B = Inf limit, where no finite-reference component is available for
#' early protection. R_mismatch = "eb_force_init" always initializes this way;
#' R_mismatch = "eb_no_init" always skips it.
#'
#' @keywords internal
#' @noRd
initialize_R_mismatch <- function(data, params, model) {
  R_mismatch <- if (!is.null(params$R_mismatch)) params$R_mismatch else "none"
  R_finite_B <- if (!is.null(model$R_finite_B)) model$R_finite_B else data$R_finite_B
  should_init <- R_mismatch == "eb_force_init" ||
                 (R_mismatch == "eb" && is.infinite(R_finite_B))
  if (!should_init || !inherits(data, c("ss", "ss_mixture")) ||
      nrow(model$alpha) < 1)
    return(model)

  model <- single_effect_update(data, params, model, 1L)
  model <- compute_R_mismatch_state(data, params, model, phase = "init_ser")
  model$R_mismatch_init <- list(
    method = "ser",
    lambda_bias = model$lambda_bias,
    B_corrected = model$B_corrected,
    Q_art = if (!is.null(model$Q_art)) model$Q_art else NA_real_,
    artifact_flag = if (!is.null(model$artifact_flag))
                      model$artifact_flag else NA
  )
  model
}

#' @keywords internal
#' @noRd
fit_R_mismatch <- function(data, params, model) {
  compute_R_mismatch_state(data, params, model, phase = "sweep")
}

# Eigen accessor for R-mismatch QC. The ordinary SS path stores data$eigen_R.
# The ss_mixture path can change R through omega, so recover the current
# mixture spectrum from panel_R when omega is available; otherwise fall
# back to the initialized X_meta crossproduct.
#' @keywords internal
get_R_mismatch_eigen <- function(data, model) {
  if (!is.null(model$eigen_R))
    return(model$eigen_R)
  if (!is.null(data$eigen_R) && !inherits(data, "ss_mixture"))
    return(data$eigen_R)
  if (inherits(data, "ss_mixture")) {
    if (!is.null(model$omega) && !is.null(data$omega_cache)) {
      eig <- eigen_from_reduced(data$omega_cache, model$omega,
                                data$K, data$p)
      eig$values <- pmax(eig$values, 0)
      return(eig)
    }
    if (!is.null(model$omega) && !is.null(data$panel_R)) {
      R_mix <- Reduce("+", Map(function(w, R) w * R, model$omega, data$panel_R))
      R_mix <- 0.5 * (R_mix + t(R_mix))
      eig <- eigen(R_mix, symmetric = TRUE)
      eig$values <- pmax(eig$values, 0)
      return(eig)
    }
    if (!is.null(data$X)) {
      R_init <- crossprod(data$X) / data$nm1
      R_init <- 0.5 * (R_init + t(R_init))
      eig <- eigen(R_init, symmetric = TRUE)
      eig$values <- pmax(eig$values, 0)
      return(eig)
    }
  }
  data$eigen_R
}

# =============================================================================
# Q_art residual R-bias artifact diagnostic
# =============================================================================

# Fraction of the fitted residual projected onto low-eigenvalue directions of R.
#   delta = max(eig_delta_abs, eig_delta_rel * max(d))
#   A_delta = {k : d_k <= delta}
#   Q_art = sum_{k in A_delta} (v_k' r_fit)^2 / sum(r_fit^2)
# This extends the column-space check used for z or Xty in the Zou et al.
# (2022) RSS likelihood: the original check asks whether the input summary
# vector lies in the non-zero eigenspace of R; Q_art asks the same question
# of the fitted residual after R_mismatch correction. A large Q_art means
# the residual still has projection in directions where the supplied R is
# nearly singular, so fine-mapping results should be treated with caution.
#
# Returns a list with Q_art (in [0, 1]), evaluable (FALSE when no
# low-eigenvalues exist or r_fit has negligible norm),
# low_eigen_count, eig_delta. Q_art is a heuristic proportion, not a
# calibrated test statistic; see archive/ld_mismatch_generativemodel.tex
# Sec. "Detecting residual R-bias artifacts".
#' @keywords internal
compute_Q_art <- function(eigen_R, r_fit, eig_delta_rel = 1e-3,
                          eig_delta_abs = 0,
                          residual_norm_floor = 1e-12) {
  d <- eigen_R$values
  delta  <- max(eig_delta_abs, eig_delta_rel * max(d))
  proj <- low_eigen_projection_fraction(eigen_R, r_fit, delta,
                                        residual_norm_floor)
  if (!proj$evaluable) {
    return(list(Q_art = 0, evaluable = FALSE,
                low_eigen_count = proj$low_eigen_count, eig_delta = delta))
  }
  list(Q_art = proj$fraction, evaluable = TRUE,
       low_eigen_count = proj$low_eigen_count, eig_delta = delta)
}
