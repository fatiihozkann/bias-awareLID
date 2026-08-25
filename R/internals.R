#' @keywords internal
.centered_W <- function(M) {
  C <- suppressWarnings(stats::cor(M, use = "pairwise.complete.obs"))
  C <- as.matrix(Matrix::nearPD(C, corr = TRUE)$mat)
  cb <- mean(C[upper.tri(C)])
  Q <- C; Q[upper.tri(Q)] <- Q[upper.tri(Q)] - cb
  Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
  A <- abs(Q); diag(A) <- 0
  k <- rowSums(A); L <- A %*% A
  W <- (L + A) / (outer(k, k, pmin) + 1 - A); diag(W) <- 0
  list(W = W, A = A, Q = Q)
}

#' @keywords internal
.permute_pvalues <- function(Xres, groups, W_obs, B, seed) {
  set.seed(seed)
  p <- ncol(Xres); cnt <- matrix(0, p, p)
  gidx <- split(seq_len(nrow(Xres)), groups)
  for (b in seq_len(B)) {
    Xp <- Xres
    for (rows in gidx) for (j in seq_len(p)) Xp[rows, j] <- Xres[sample(rows), j]
    cnt <- cnt + (.centered_W(Xp)$W >= W_obs)
  }
  (1 + cnt) / (B + 1)
}

#' @keywords internal
.bh_select <- function(pv, ut, q) {
  pu <- pv[ut]; m <- length(pu); ps <- sort(pu)
  crit <- seq_len(m) * q / m
  k <- which(ps <= crit)
  bh <- if (length(k)) ps[max(k)] else 0
  pv <= bh & ut
}

#' @keywords internal
.beta_select <- function(W, P, E0, tau, floor,
                         grid = c(0, .1, .25, .5, .75, 1, 1.25, 1.5, 2, 3, 5)) {
  nE0 <- sum(E0)
  if (nE0 == 0) return(list(beta = 0, final = E0))
  ret <- vapply(grid, function(b) sum(E0 & (pmax(0, W - b * P) > tau)), 0L)
  ok <- grid[ret >= ceiling(floor * nE0)]
  beta <- if (length(ok)) max(ok) else 0
  list(beta = beta, final = E0 & (pmax(0, W - beta * P) > tau))
}
