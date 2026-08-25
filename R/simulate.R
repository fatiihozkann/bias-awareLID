#' Simulate data with designed local-dependence truth classes
#'
#' Generates multi-group item responses with planted dependent pairs: half
#' consistent (dependence in every group) and half group-specific (dependence
#' in one randomly chosen affected group only), following the design of the
#' validation simulation in Ozkan (2026).
#'
#' @param N total sample size. @param p number of items.
#' @param rho residual dependence strength. @param n_groups number of groups.
#' @param loading common factor loading (default .6). @param seed RNG seed.
#' @return list: \code{data}, \code{group}, \code{truth} (data frame of
#'   planted pairs with class and affected group).
#' @export
simulate_lid_data <- function(N = 500, p = 20, rho = .3, n_groups = 2,
                              loading = .6, seed = 1) {
  set.seed(seed)
  k_pairs <- max(1, floor(p / 5))
  perm <- sample.int(p, p)
  pairs <- split(perm[1:(2 * k_pairs)], rep(1:k_pairs, each = 2))
  pairs <- lapply(pairs, function(v) sort(as.integer(v)))
  n_gs <- floor(k_pairs / 2)
  gs_idx <- if (n_gs > 0) sample(seq_len(k_pairs), n_gs) else integer(0)
  affected <- sample.int(n_groups, 1)
  sizes <- rep(floor(N / n_groups), n_groups)
  sizes[1] <- sizes[1] + N - sum(sizes)
  sim_group <- function(n, active) {
    eta <- stats::rnorm(n)
    Psi <- diag(1 - loading^2, p)
    for (pr in active) Psi[pr[1], pr[2]] <- Psi[pr[2], pr[1]] <-
        rho * (1 - loading^2)
    ev <- eigen(Psi, symmetric = TRUE)
    Ehalf <- ev$vectors %*% diag(sqrt(pmax(ev$values, 1e-8))) %*% t(ev$vectors)
    E <- matrix(stats::rnorm(n * p), n, p) %*% Ehalf
    outer(eta, rep(loading, p)) + E
  }
  dfl <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    active <- c(pairs[setdiff(seq_len(k_pairs), gs_idx)],
                if (g == affected) pairs[gs_idx] else list())
    dfl[[g]] <- sim_group(sizes[g], active)
  }
  X <- do.call(rbind, dfl); colnames(X) <- paste0("X", seq_len(p))
  truth <- do.call(rbind, lapply(seq_len(k_pairs), function(i) data.frame(
    item_i = paste0("X", pairs[[i]][1]), item_j = paste0("X", pairs[[i]][2]),
    class = if (i %in% gs_idx) "group_specific" else "consistent",
    affected_group = if (i %in% gs_idx) affected else NA_integer_)))
  list(data = X, group = factor(rep(seq_len(n_groups), sizes)), truth = truth)
}
