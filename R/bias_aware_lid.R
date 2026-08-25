#' Bias-aware detection of group-specific local item dependence
#'
#' Runs the full centered bias-aware pipeline: CFA residualization, signed
#' centering of residual correlations, weighted topological overlap, within-group
#' permutation inference with BH-FDR control, the qualified set E0, the
#' retention-constrained group-disparity penalty, and signed reversal
#' diagnostics. Implements the pipeline of Ozkan (2026).
#'
#' @param data item response data (rows = persons, columns = items).
#' @param group factor or vector of group memberships (length nrow(data)).
#' @param model optional lavaan measurement model; \code{NULL} = one factor.
#' @param B number of permutations (default 1000).
#' @param q BH false discovery level (default .05).
#' @param tau unpenalized magnitude floor (default .10).
#' @param floor retention proportion for beta selection (default .90).
#' @param seed permutation seed (default 42).
#' @return object of class \code{biasaware_lid}.
#' @export
bias_aware_lid <- function(data, group, model = NULL, B = 1000, q = .05,
                           tau = .10, floor = .90, seed = 42) {
  X <- as.matrix(data); items <- colnames(X)
  group <- factor(group)
  stopifnot(length(group) == nrow(X), nlevels(group) >= 2)
  rz <- cfa_residualize(X, model)
  Xres <- rz$residuals
  obs <- .centered_W(Xres); W <- obs$W
  pv <- .permute_pvalues(Xres, group, W, B, seed)
  ut <- upper.tri(W)
  S <- .bh_select(pv, ut, q)
  E0 <- S & (W > tau)
  gidx <- split(seq_len(nrow(X)), group)
  cwg <- lapply(gidx, function(ix) .centered_W(Xres[ix, , drop = FALSE]))
  Wg <- lapply(cwg, `[[`, "W"); Qg <- lapply(cwg, `[[`, "Q")
  P <- Reduce(pmax, lapply(Wg, function(w) abs(w - W)))
  sel <- .beta_select(W, P, E0, tau, floor)
  final <- sel$final
  rev_ <- (Reduce(pmin, Qg) < 0) & (Reduce(pmax, Qg) > 0)
  ij <- which(ut, arr.ind = TRUE)
  edges <- data.frame(item_i = items[ij[, 1]], item_j = items[ij[, 2]],
    wTO = W[ij], Q3_star_signed = obs$Q[ij], p_value = pv[ij],
    group_disparity = P[ij], BH = S[ut], E0 = E0[ut],
    final = final[ut], sign_reversal = rev_[ut])
  structure(list(edges = edges, beta_star = sel$beta,
    n_bh = sum(S), n_E0 = sum(E0), n_final = sum(final),
    removed = edges[edges$E0 & !edges$final, c("item_i", "item_j")],
    reversals_E0 = sum(rev_ & E0), reversals_final = sum(rev_ & final),
    W = W, P = P, Q_pooled = obs$Q, Q_groups = Qg, fit = rz$fit,
    settings = list(B = B, q = q, tau = tau, floor = floor, seed = seed,
                    groups = table(group))), class = "biasaware_lid")
}
