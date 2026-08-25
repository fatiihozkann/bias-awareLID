#' Stratified full-pipeline bootstrap stability
#'
#' Repeats the complete bias-aware pipeline (CFA refit, permutation inference,
#' E0 construction, beta reselection) in stratified person-level resamples that
#' preserve the original group sizes.
#'
#' @param data,group,model,B,q,tau,floor as in \code{\link{bias_aware_lid}}.
#' @param reps number of bootstrap replications (default 200).
#' @param seed base seed; replication r uses seed + r.
#' @param criteria stability criteria (default c(.80, .90)).
#' @return list: \code{edge_stability} (per-pair final-selection frequency and
#'   stability flags) and \code{replications} (per-replication summary).
#' @export
lid_bootstrap <- function(data, group, model = NULL, reps = 200, B = 1000,
                          q = .05, tau = .10, floor = .90, seed = 880000,
                          criteria = c(.80, .90)) {
  X <- as.matrix(data); group <- factor(group)
  gidx <- split(seq_len(nrow(X)), group)
  freq <- NULL; summ <- list()
  for (r in seq_len(reps)) {
    set.seed(seed + r)
    idx <- unlist(lapply(gidx, function(ix) sample(ix, length(ix), TRUE)))
    res <- tryCatch(bias_aware_lid(X[idx, , drop = FALSE], group[idx], model,
              B = B, q = q, tau = tau, floor = floor, seed = seed + 100000 + r),
              error = function(e) NULL)
    if (is.null(res)) { summ[[r]] <- data.frame(rep = r, ok = FALSE); next }
    f <- as.integer(res$edges$final)
    freq <- if (is.null(freq)) f else freq + f
    summ[[r]] <- data.frame(rep = r, ok = TRUE, beta = res$beta_star,
      n_bh = res$n_bh, n_E0 = res$n_E0, n_final = res$n_final)
  }
  base <- bias_aware_lid(X, group, model, B = B, q = q, tau = tau,
                         floor = floor, seed = seed)
  ok_n <- sum(vapply(summ, function(s) isTRUE(s$ok), NA))
  es <- data.frame(base$edges[, c("item_i", "item_j")],
                   final_frequency = freq / ok_n)
  for (cr in criteria) es[[paste0("stable_", cr * 100)]] <- es$final_frequency >= cr
  list(edge_stability = es, replications = do.call(rbind, summ),
       baseline = base)
}
