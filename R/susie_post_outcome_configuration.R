# Post-hoc causal-configuration probabilities for one or more SuSiE-class fits.
#
# Two algorithms live here, exposed through one entry point:
#
#   * SuSiEx (Nature Genetics, 2024): N-trait 2^N enumeration. Per CS tuple
#     (one CS chosen from each trait), report posterior probabilities over
#     all 2^N "which traits share the causal" patterns plus per-trait
#     marginals. Legacy reference implementation:
#     `mvf.susie.alpha::posthoc_multfsusie`.
#
#   * Coloc pairwise ABF (Wallace, 2020 / `coloc::coloc.bf_bf`): pairwise
#     H0/H1/H2/H3/H4 posteriors for every (trait, trait') pair across every
#     (CS in trait, CS in trait') pair. Implemented inline here as a
#     verbatim port of `coloc:::combine.abf` so susieR has no soft
#     dependency on coloc.
#
# This file deliberately does NOT introduce S3 generics. The public function
# normalises any supported input shape (single fit, list of fits, or a single
# multi-output fit treated outcome-wise) to a flat list of "trait views",
# then runs the requested algorithms against that list. Class-aware branches
# use `inherits()` and are confined to one helper.

#' Post-hoc causal-configuration probabilities for one or more SuSiE-class fits
#'
#' Combines two complementary post-hoc analyses into a single call:
#' (1) the SuSiEx \eqn{2^N} combinatorial enumeration, reporting the
#' posterior probability of every binary causality pattern across the
#' \eqn{N} input traits; and (2) the coloc pairwise ABF, reporting the
#' five colocalisation hypothesis posteriors (H0/H1/H2/H3/H4) for every
#' pair of traits. Both run by default; pick a subset via \code{methods}.
#'
#' Two grouping modes are supported through the \code{by} argument:
#' \describe{
#'   \item{\code{"fit"}}{Each input fit contributes a single trait view.
#'     Multi-output fits (\code{mvsusie}, \code{mfsusie}) are kept whole: the
#'     trait's per-(CS, SNP) log Bayes factors are the joint composite
#'     stored on the fit as \code{lbf_variable}. Configuration enumeration
#'     loops over the cross-product \eqn{L_1 \times \dots \times L_N} of CS
#'     indices.}
#'   \item{\code{"outcome"}}{Multi-output fits fan out into per-outcome views,
#'     each with its own per-(CS, SNP) log Bayes factors read from
#'     \code{fit$lbf_outcome} (an \eqn{L \times J \times R} or
#'     \eqn{L \times J \times M} array). All per-outcome views share the
#'     joint fit's PIP matrix and CS list, so the configuration enumeration
#'     reduces to a single index \eqn{l \in 1..L}. Single-output \code{susie}
#'     fits pass through unchanged. Requires \code{$lbf_outcome} on the
#'     fit (set \code{attach_lbf_outcome = TRUE} when fitting).}
#' }
#'
#' \subsection{SuSiEx algorithm}{
#' For each credible-set tuple \eqn{(l_1, \dots, l_N)}:
#' \enumerate{
#'   \item Per-trait CS-level log BF (alpha-weighted SNP average):
#'     \deqn{\log\mathrm{BF}^{(n)}_{l_n} = \sum_j \alpha_{n,l_n,j}\,
#'       \log\mathrm{BF}_{n,l_n,j}.}
#'   \item Enumerate the \eqn{2^N} binary configurations
#'     \eqn{c \in \{0,1\}^N}.
#'   \item Configuration log BF:
#'     \deqn{\log\mathrm{BF}^{(c)} = \sum_{n: c_n = 1} \log\mathrm{BF}^{(n)}_{l_n}.}
#'   \item Normalise under a uniform prior over the \eqn{2^N} configurations.
#'   \item Per-trait marginal: \eqn{P(\mathrm{trait}\,n\,\mathrm{causal}) =
#'     \sum_{c: c_n = 1} P(c \mid \mathrm{tuple})}.
#' }
#' }
#'
#' \subsection{Coloc pairwise algorithm}{
#' For each unordered trait pair \eqn{(n, n')} and each CS pair
#' \eqn{(l_n, l_{n'})}, with per-SNP log BFs
#' \eqn{\ell_1 = \log\mathrm{BF}_{n,l_n,\cdot}} and
#' \eqn{\ell_2 = \log\mathrm{BF}_{n',l_{n'},\cdot}} (length \eqn{J}), the
#' five hypothesis log-BFs are
#' \deqn{\log\mathrm{BF}_{H_0} = 0,\quad
#'       \log\mathrm{BF}_{H_1} = \log p_1 + \mathrm{LSE}(\ell_1),\quad
#'       \log\mathrm{BF}_{H_2} = \log p_2 + \mathrm{LSE}(\ell_2),}
#' \deqn{\log\mathrm{BF}_{H_3} = \log p_1 + \log p_2 +
#'       \mathrm{logdiff}(\mathrm{LSE}(\ell_1) + \mathrm{LSE}(\ell_2),\;
#'                        \mathrm{LSE}(\ell_1 + \ell_2)),}
#' \deqn{\log\mathrm{BF}_{H_4} = \log p_{12} + \mathrm{LSE}(\ell_1 + \ell_2),}
#' and the corresponding posteriors are
#' \eqn{\mathrm{PP.H}_h = \exp(\log\mathrm{BF}_{H_h} -
#'       \mathrm{LSE}(\log\mathrm{BF}_{H_0:H_4}))}, where
#' \eqn{\mathrm{LSE}} is the log-sum-exp.
#' \itemize{
#'   \item H0: no causal variant in either CS.
#'   \item H1: causal in trait \eqn{n} only.
#'   \item H2: causal in trait \eqn{n'} only.
#'   \item H3: distinct causals in the two traits.
#'   \item H4: a single shared causal variant.
#' }
#' }
#'
#' @param input A single fit of class \code{susie}, \code{mvsusie}, or
#'   \code{mfsusie}, OR a list of such fits.
#' @param by Either \code{"fit"} (one trait per input fit; default) or
#'   \code{"outcome"} (multi-output fits expand into per-outcome traits).
#' @param methods Character vector. Subset of \code{c("susiex",
#'   "coloc_pairwise")}. By default both run.
#' @param prob_thresh Per-trait marginal threshold for the convenience
#'   \code{$active} flags in the SuSiEx output. Default \code{0.8}.
#' @param cs_only Logical. If \code{TRUE} (default) only enumerate over CSs
#'   present in each fit's \code{$sets$cs}; if \code{FALSE} loop over all L
#'   rows of \code{$alpha}. Either way, effects whose entire alpha row is
#'   zero are skipped. When \code{TRUE}, every fit must carry a non-null
#'   \code{$sets$cs} or the function errors.
#' @param p1,p2,p12 Coloc per-SNP causal priors: \code{p1} for trait 1
#'   alone, \code{p2} for trait 2 alone, \code{p12} for shared causal.
#'   Defaults match \code{coloc::coloc.bf_bf}: \code{p1 = p2 = 1e-4},
#'   \code{p12 = 5e-6}. Only used when \code{"coloc_pairwise"} is in
#'   \code{methods}.
#' @param ... Currently ignored.
#'
#' @return A list with one element per requested method:
#' \describe{
#'   \item{\code{$susiex}}{A list of length equal to the number of CS tuples
#'     considered. Each element has components \code{cs_indices} (length-N
#'     integer tuple), \code{logBF_trait} (length N), \code{configs}
#'     (\eqn{2^N \times N} binary matrix), \code{config_prob} (length
#'     \eqn{2^N}), \code{posthoc} (length-N marginals), and \code{active}
#'     (logical at \code{prob_thresh}).}
#'   \item{\code{$coloc_pairwise}}{A data.frame with one row per
#'     (trait1, trait2, l1, l2) combination, columns \code{trait1, trait2,
#'     l1, l2, hit1, hit2, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4}.}
#' }
#'
#' @references
#' SuSiEx, Nature Genetics 2024 (combinatorial \eqn{2^N} enumeration).
#' Wallace, PLoS Genetics 2020 (coloc pairwise H0/H1/H2/H3/H4 ABF).
#'
#' @export
susie_post_outcome_configuration <- function(input,
                                            by          = c("fit", "outcome"),
                                            methods     = c("susiex",
                                                            "coloc_pairwise"),
                                            prob_thresh = 0.8,
                                            cs_only     = TRUE,
                                            p1          = 1e-4,
                                            p2          = 1e-4,
                                            p12         = 5e-6,
                                            ...) {
  by      <- match.arg(by)
  methods <- match.arg(methods, c("susiex", "coloc_pairwise"),
                       several.ok = TRUE)

  if (!is.numeric(prob_thresh) || length(prob_thresh) != 1L ||
      !is.finite(prob_thresh) || prob_thresh < 0 || prob_thresh > 1) {
    stop("`prob_thresh` must be a single numeric in [0, 1].")
  }
  if (!is.logical(cs_only) || length(cs_only) != 1L || is.na(cs_only)) {
    stop("`cs_only` must be a single logical (TRUE or FALSE).")
  }
  for (nm in c("p1", "p2", "p12")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
        v <= 0 || v >= 1) {
      stop("`", nm, "` must be a single numeric in (0, 1).")
    }
  }

  views <- normalise_to_views(input, by = by, cs_only = cs_only)

  out <- list()
  if ("susiex" %in% methods) {
    out$susiex <- susiex_configurations(views,
                                        by          = by,
                                        prob_thresh = prob_thresh)
  }
  if ("coloc_pairwise" %in% methods) {
    out$coloc_pairwise <- coloc_pairwise_abf(views,
                                             p1  = p1,
                                             p2  = p2,
                                             p12 = p12)
  }
  out
}

# -----------------------------------------------------------------------------
# Input normalisation
# -----------------------------------------------------------------------------

is_susie_fit <- function(x) {
  inherits(x, "susie") || inherits(x, "mvsusie") || inherits(x, "mfsusie")
}

normalise_to_views <- function(input, by, cs_only) {
  fits <- if (is_susie_fit(input)) list(input) else as.list(input)

  if (length(fits) == 0L) {
    stop("`input` must be a SuSiE-class fit or a non-empty list of fits.")
  }
  for (k in seq_along(fits)) {
    if (!is_susie_fit(fits[[k]])) {
      stop("Element ", k,
           " of `input` is not a SuSiE-class fit (`susie`, `mvsusie`, or ",
           "`mfsusie`).")
    }
    if (cs_only && is.null(fits[[k]]$sets$cs)) {
      stop("Fit ", k, ": `cs_only = TRUE` requires `$sets$cs` to be present. ",
           "Either pass `cs_only = FALSE` or attach a credible-set list via ",
           "susie_get_cs() before calling.")
    }
  }

  raw_names <- names(fits)
  if (is.null(raw_names)) raw_names <- character(length(fits))
  default_names <- paste0("trait_", seq_along(fits))
  raw_names[!nzchar(raw_names)] <- default_names[!nzchar(raw_names)]

  views <- vector("list", 0)
  for (k in seq_along(fits)) {
    views <- c(views, expand_one_fit(fits[[k]], raw_names[k], by = by))
  }
  views
}

expand_one_fit <- function(fit, base_name, by) {
  if (by == "fit") {
    return(list(make_view(
      name    = base_name,
      alpha   = fit$alpha,
      lbf     = fit$lbf_variable,
      sets_cs = fit$sets$cs
    )))
  }

  # by = "outcome": multi-output fits fan out; single-output fits pass
  # through as one view.
  if (inherits(fit, "mvsusie") || inherits(fit, "mfsusie")) {
    if (is.null(fit$lbf_outcome)) {
      stop("Fit '", base_name, "': `by = \"outcome\"` requires `$lbf_outcome` ",
           "(an L x J x R or L x J x M array) on the fit. ",
           "Refit with `attach_lbf_outcome = TRUE` (the default in mfsusie / ",
           "mvsusie), or pass `by = \"fit\"` to use the joint composite log ",
           "BF instead.")
    }
    R <- dim(fit$lbf_outcome)[3L]
    out_names <- dimnames(fit$lbf_outcome)[[3L]]
    if (is.null(out_names)) out_names <- paste0("outcome_", seq_len(R))
    views <- vector("list", R)
    for (r in seq_len(R)) {
      views[[r]] <- make_view(
        name    = paste0(base_name, "_", out_names[r]),
        alpha   = fit$alpha,
        lbf     = fit$lbf_outcome[, , r, drop = TRUE],
        sets_cs = fit$sets$cs
      )
    }
    return(views)
  }

  # Single-output `susie` under by = "outcome": same as by = "fit".
  list(make_view(
    name    = base_name,
    alpha   = fit$alpha,
    lbf     = fit$lbf_variable,
    sets_cs = fit$sets$cs
  ))
}

make_view <- function(name, alpha, lbf, sets_cs) {
  if (is.null(alpha) || is.null(lbf)) {
    stop("Trait '", name, "': both `$alpha` and `$lbf_variable` (or per-",
         "outcome lbf row) must be non-null.")
  }
  if (!is.matrix(alpha)) alpha <- as.matrix(alpha)
  if (!is.matrix(lbf))   lbf   <- as.matrix(lbf)
  if (!identical(dim(alpha), dim(lbf))) {
    stop("Trait '", name, "': `alpha` and `lbf` must have identical shape; ",
         "got ", paste(dim(alpha), collapse = "x"), " vs ",
         paste(dim(lbf), collapse = "x"), ".")
  }
  list(name = name, alpha = alpha, lbf = lbf, sets_cs = sets_cs)
}

# -----------------------------------------------------------------------------
# CS-tuple enumeration shared by both algorithms.
# -----------------------------------------------------------------------------

# Per-view CS index set, restricted to $sets$cs when cs_only = TRUE.
view_cs_indices <- function(view, cs_only) {
  L_n <- nrow(view$alpha)
  if (!cs_only) return(seq_len(L_n))

  idx <- attr(view$sets_cs, "cs_idx")
  if (is.null(idx)) {
    # Fall back to the names of $sets$cs ("L1", "L2", ... in susieR's format).
    if (length(view$sets_cs) > 0L && !is.null(names(view$sets_cs))) {
      idx <- as.integer(sub("^L", "", names(view$sets_cs)))
    } else {
      idx <- seq_len(L_n)
    }
  }
  idx[idx >= 1L & idx <= L_n]
}

# Returns a list of length-N integer tuples (one CS index per view).
# Under by = "outcome" all views share CSs and we use the diagonal.
# Under by = "fit" we use the cross-product.
enumerate_cs_tuples <- function(views, by, cs_only) {
  per_view <- lapply(views, view_cs_indices, cs_only = cs_only)
  if (any(vapply(per_view, length, integer(1)) == 0L)) return(list())

  if (by == "outcome") {
    common <- Reduce(intersect, per_view)
    lapply(common, function(l) rep.int(l, length(views)))
  } else {
    grid <- expand.grid(per_view, KEEP.OUT.ATTRS = FALSE)
    lapply(seq_len(nrow(grid)), function(i) as.integer(grid[i, , drop = TRUE]))
  }
}

# -----------------------------------------------------------------------------
# SuSiEx 2^N configuration enumeration.
# -----------------------------------------------------------------------------

susiex_configurations <- function(views, by, prob_thresh,
                                  max_traits = 20L) {
  N <- length(views)
  if (N > max_traits) {
    stop("susiex: N = ", N, " traits exceeds the safety ceiling (",
         max_traits, "); 2^N enumeration would be too large.")
  }

  cs_tuples <- enumerate_cs_tuples(views, by = by, cs_only = TRUE)
  if (length(cs_tuples) == 0L) return(list())

  configs <- as.matrix(expand.grid(rep(list(c(0L, 1L)), N)))
  colnames(configs) <- paste0("trait_", seq_len(N))
  trait_names <- vapply(views, function(v) v$name, character(1))

  out <- vector("list", length(cs_tuples))
  for (ti in seq_along(cs_tuples)) {
    tuple <- cs_tuples[[ti]]

    logBF_trait <- numeric(N)
    skip <- FALSE
    for (n in seq_len(N)) {
      l_n       <- tuple[n]
      alpha_row <- views[[n]]$alpha[l_n, ]
      lbf_row   <- views[[n]]$lbf  [l_n, ]
      if (all(alpha_row == 0)) { skip <- TRUE; break }
      logBF_trait[n] <- sum(alpha_row * lbf_row)   # alpha-weighted SNP avg
    }
    if (skip) {
      out[[ti]] <- NULL
      next
    }

    logBF_conf <- as.vector(configs %*% logBF_trait)
    maxlog     <- max(logBF_conf)
    prob_conf  <- exp(logBF_conf - maxlog)
    prob_conf  <- prob_conf / sum(prob_conf)
    posthoc    <- as.vector(crossprod(configs, prob_conf))

    out[[ti]] <- list(
      cs_indices  = setNames(as.integer(tuple), trait_names),
      logBF_trait = setNames(logBF_trait,        trait_names),
      configs     = configs,
      config_prob = prob_conf,
      posthoc     = setNames(posthoc,            trait_names),
      active      = setNames(posthoc >= prob_thresh, trait_names)
    )
  }

  out[!vapply(out, is.null, logical(1))]
}

# -----------------------------------------------------------------------------
# Coloc pairwise ABF (verbatim port of coloc:::combine.abf).
# -----------------------------------------------------------------------------

# Numerically stable log(sum(exp(x))).
.logsum <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Numerically stable log(exp(a) - exp(b)) for a > b.
.logdiff <- function(a, b) {
  m <- max(a, b)
  m + log(exp(a - m) - exp(b - m))
}

# Compute (PP.H0, PP.H1, PP.H2, PP.H3, PP.H4) from per-SNP log-BF vectors,
# matching coloc:::combine.abf line-for-line.
combine_abf_pair <- function(l1, l2, p1, p2, p12) {
  stopifnot(length(l1) == length(l2))
  lsum    <- l1 + l2
  lH0     <- 0
  lH1     <- log(p1)  + .logsum(l1)
  lH2     <- log(p2)  + .logsum(l2)
  lH3     <- log(p1)  + log(p2) +
             .logdiff(.logsum(l1) + .logsum(l2), .logsum(lsum))
  lH4     <- log(p12) + .logsum(lsum)
  all_lH  <- c(lH0, lH1, lH2, lH3, lH4)
  pp      <- exp(all_lH - .logsum(all_lH))
  names(pp) <- paste0("PP.H", 0:4)
  pp
}

coloc_pairwise_abf <- function(views, p1, p2, p12) {
  N <- length(views)
  empty <- data.frame(trait1 = character(0), trait2 = character(0),
                      l1 = integer(0), l2 = integer(0),
                      hit1 = character(0), hit2 = character(0),
                      PP.H0 = numeric(0), PP.H1 = numeric(0),
                      PP.H2 = numeric(0), PP.H3 = numeric(0),
                      PP.H4 = numeric(0),
                      stringsAsFactors = FALSE)
  if (N < 2L) return(empty)

  trait_names <- vapply(views, function(v) v$name, character(1))
  rows <- list()

  for (a in seq_len(N - 1L)) {
    for (b in (a + 1L):N) {
      L1 <- view_cs_indices(views[[a]], cs_only = TRUE)
      L2 <- view_cs_indices(views[[b]], cs_only = TRUE)
      if (length(L1) == 0L || length(L2) == 0L) next

      var_names_a <- colnames(views[[a]]$lbf)
      var_names_b <- colnames(views[[b]]$lbf)

      for (i in L1) {
        if (all(views[[a]]$alpha[i, ] == 0)) next
        l1_row <- views[[a]]$lbf[i, ]
        for (j in L2) {
          if (all(views[[b]]$alpha[j, ] == 0)) next
          l2_row <- views[[b]]$lbf[j, ]

          pp <- combine_abf_pair(l1_row, l2_row, p1 = p1, p2 = p2, p12 = p12)

          hit1 <- if (!is.null(var_names_a)) {
                    var_names_a[which.max(l1_row)]
                  } else {
                    paste0("snp_", which.max(l1_row))
                  }
          hit2 <- if (!is.null(var_names_b)) {
                    var_names_b[which.max(l2_row)]
                  } else {
                    paste0("snp_", which.max(l2_row))
                  }

          rows[[length(rows) + 1L]] <- data.frame(
            trait1 = trait_names[a], trait2 = trait_names[b],
            l1     = i,              l2     = j,
            hit1   = hit1,           hit2   = hit2,
            PP.H0  = pp["PP.H0"],    PP.H1  = pp["PP.H1"],
            PP.H2  = pp["PP.H2"],    PP.H3  = pp["PP.H3"],
            PP.H4  = pp["PP.H4"],
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        }
      }
    }
  }

  if (length(rows) == 0L) return(empty)
  do.call(rbind, rows)
}
