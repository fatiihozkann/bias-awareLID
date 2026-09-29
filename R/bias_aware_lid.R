#' Bias-aware detection of group-specific local item dependence
#'
#' Centered bias-aware pipeline with three options added in version 0.2.0:
#' \code{centering = "model"} centers each residual correlation by its
#' model-implied null baseline (recommended; exact under local independence for
#' any loading pattern), \code{null = "parametric"} draws the reference
#' distribution from the fitted measurement model instead of column
#' permutation, and \code{heterogeneity_test = TRUE} tests each qualified
#' pair's group disparity against a group-label permutation null with BH control.
#'
#' @inheritParams cfa_residualize
#' @param group factor of group memberships. @param B permutations. @param q BH level.
#' @param tau magnitude floor. @param floor retention proportion. @param seed seed.
#' @param centering "model" (default) or "mean".
#' @param null "permutation" (default) or "parametric".
#' @param heterogeneity_test logical; test disparity of qualified pairs.
#' @param refit logical; re-estimate the measurement model in every parametric replicate.
#' @param direct_test logical; with \code{null = "parametric"}, test each final pair's direct
#'   residual correlation |Q3*| against the model-based reference with BH control.
#' @return object of class \code{biasaware_lid}.
#' @export
bias_aware_lid <- function(data, group, model = NULL, B = 1000, q = .05,
                           tau = .10, floor = .90, seed = 42,
                           centering = c("model", "mean"),
                           null = c("permutation", "parametric"),
                           heterogeneity_test = FALSE, refit = FALSE, direct_test = FALSE) {
  centering <- match.arg(centering); null <- match.arg(null)
  X <- as.matrix(data); items <- colnames(X); group <- factor(group)
  stopifnot(length(group) == nrow(X), nlevels(group) >= 2)
  rz <- cfa_residualize(X, model); Xres <- rz$residuals
  base <- if (centering == "model") model_baseline(rz$fit) else NULL
  obs <- .centered_W(Xres, base); W <- obs$W
  pQ <- NULL
  if (null == "permutation") pv <- .permute_pvalues(Xres, group, W, B, seed) else {
    pr <- .parametric_pvalues(rz$fit, nrow(X), W, B, seed, base, refit = refit,
                              model = if (is.null(model)) paste0("F1 =~ ", paste(items, collapse = " + ")) else model,
                              Q_obs = if (direct_test) obs$Q else NULL, items = items)
    pv <- pr$pW; pQ <- pr$pQ }
  ut <- upper.tri(W); S <- .bh_select(pv, ut, q); E0 <- S & (W > tau)
  gidx <- split(seq_len(nrow(X)), group)
  cwg <- lapply(gidx, function(ix) .centered_W(Xres[ix, , drop = FALSE], base))
  Wg <- lapply(cwg, `[[`, "W"); Qg <- lapply(cwg, `[[`, "Q")
  P <- Reduce(pmax, lapply(Wg, function(w) abs(w - W)))
  sel <- .beta_select(W, P, E0, tau, floor); final <- sel$final
  rev_ <- (Reduce(pmin, Qg) < 0) & (Reduce(pmax, Qg) > 0)
  het <- NULL
  if (heterogeneity_test && sum(E0) > 0) {
    h <- .heterogeneity_pvalues(Xres, group, W, E0, base, B, seed + 7L)
    m <- nrow(h); ps <- sort(h$p_het); k <- which(ps <= seq_len(m) * q / m)
    h$het_BH <- h$p_het <= (if (length(k)) ps[max(k)] else 0)
    h$item_i <- items[h$i]; h$item_j <- items[h$j]; het <- h[, c("item_i", "item_j", "P", "p_het", "het_BH")]
  }
  direct <- rep(NA, sum(ut)); pdirect <- rep(NA_real_, sum(ut))
  if (direct_test && !is.null(pQ) && sum(final) > 0) {
    pdirect <- pQ[ut]; f <- final[ut]; pf <- pdirect[f]; mf <- length(pf)
    ps_ <- sort(pf); k <- which(ps_ <= seq_len(mf) * q / mf); cut <- if (length(k)) ps_[max(k)] else 0
    direct <- f & (pdirect <= cut) }
  ij <- which(ut, arr.ind = TRUE)
  edges <- data.frame(item_i = items[ij[, 1]], item_j = items[ij[, 2]],
    wTO = W[ij], Q3_star_signed = obs$Q[ij], p_value = pv[ij], group_disparity = P[ij],
    BH = S[ut], E0 = E0[ut], final = final[ut], sign_reversal = rev_[ut],
    p_direct = pdirect, direct_BH = direct)
  structure(list(edges = edges, beta_star = sel$beta, n_bh = sum(S), n_E0 = sum(E0),
    n_final = sum(final), removed = edges[edges$E0 & !edges$final, c("item_i", "item_j")],
    reversals_E0 = sum(rev_ & E0), reversals_final = sum(rev_ & final),
    heterogeneity = het, W = W, P = P, W_groups = Wg, Q_pooled = obs$Q, Q_groups = Qg, fit = rz$fit,
    settings = list(B = B, q = q, tau = tau, floor = floor, seed = seed, centering = centering,
                    null = null, refit = refit, direct_test = direct_test, groups = table(group))), class = "biasaware_lid")
}
