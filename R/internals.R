#' @keywords internal
.centered_W <- function(M, baseline = NULL) {
  C <- suppressWarnings(stats::cor(M, use = "pairwise.complete.obs"))
  C <- as.matrix(Matrix::nearPD(C, corr = TRUE)$mat)
  if (is.null(baseline)) {
    cb <- mean(C[upper.tri(C)]); Q <- C; Q[upper.tri(Q)] <- Q[upper.tri(Q)] - cb
    Q[lower.tri(Q)] <- t(Q)[lower.tri(Q)]
  } else { Q <- C - baseline }
  diag(Q) <- 0; A <- abs(Q); diag(A) <- 0
  k <- rowSums(A); L <- A %*% A
  W <- (L + A) / (outer(k, k, pmin) + 1 - A); diag(W) <- 0
  list(W = W, A = A, Q = Q)
}

#' Model-implied null baseline of residual correlations
#'
#' Population correlation matrix of Bartlett-score residuals under local
#' independence for a fitted lavaan model: cov2cor(Theta - Lambda (Lambda'
#' Theta^-1 Lambda)^-1 Lambda'). Used to center each pair by its own baseline.
#' @param fit a fitted lavaan object. @return p x p matrix with zero diagonal.
#' @export
model_baseline <- function(fit) {
  est <- lavaan::lavInspect(fit, "est")
  L <- as.matrix(est$lambda); Th <- diag(diag(as.matrix(est$theta)))
  Ti <- solve(Th)
  cov0 <- Th - L %*% solve(t(L) %*% Ti %*% L) %*% t(L)
  C0 <- stats::cov2cor(cov0); diag(C0) <- 0; C0
}

#' @keywords internal
.permute_pvalues <- function(Xres, groups, W_obs, B, seed) {
  set.seed(seed); p <- ncol(Xres); cnt <- matrix(0, p, p)
  gidx <- split(seq_len(nrow(Xres)), groups)
  for (b in seq_len(B)) {
    Xp <- Xres
    for (rows in gidx) for (j in seq_len(p)) Xp[rows, j] <- Xres[sample(rows), j]
    cnt <- cnt + (.centered_W(Xp)$W >= W_obs)
  }
  (1 + cnt) / (B + 1)
}

#' @keywords internal
.parametric_pvalues <- function(fit, n, W_obs, B, seed, baseline, refit = FALSE,
                                model = NULL, Q_obs = NULL, items = NULL) {
  set.seed(seed); est <- lavaan::lavInspect(fit, "est")
  L <- as.matrix(est$lambda); Th <- diag(diag(as.matrix(est$theta))); Ph <- as.matrix(est$psi)
  p <- nrow(L); m <- ncol(L); cntW <- matrix(0, p, p); cntQ <- matrix(0, p, p)
  Wb <- solve(t(L) %*% solve(Th) %*% L) %*% t(L) %*% solve(Th)
  Phalf <- t(chol(Ph)); Thalf <- sqrt(diag(Th)); nu <- as.numeric(est$nu)
  for (b in seq_len(B)) {
    eta <- matrix(stats::rnorm(n * m), n, m) %*% t(Phalf)
    X <- eta %*% t(L) + sweep(matrix(stats::rnorm(n * p), n, p), 2, Thalf, `*`)
    X <- sweep(X, 2, nu, `+`); colnames(X) <- items
    if (refit) {
      fb <- tryCatch(lavaan::cfa(model, data = as.data.frame(X), estimator = "MLR", missing = "fiml"),
                     error = function(e) NULL)
      if (is.null(fb) || !lavaan::lavInspect(fb, "converged")) { fs <- X %*% t(Wb); bl <- baseline }
      else { fs <- lavaan::lavPredict(fb, method = "Bartlett"); bl <- model_baseline(fb) }
    } else { fs <- X %*% t(Wb); bl <- baseline }
    E <- X
    for (j in seq_len(p)) E[, j] <- stats::resid(stats::lm(X[, j] ~ fs))
    cw <- .centered_W(E, bl)
    cntW <- cntW + (cw$W >= W_obs)
    if (!is.null(Q_obs)) cntQ <- cntQ + (abs(cw$Q) >= abs(Q_obs))
  }
  list(pW = (1 + cntW) / (B + 1), pQ = if (is.null(Q_obs)) NULL else (1 + cntQ) / (B + 1))
}

#' @keywords internal
.bh_select <- function(pv, ut, q) {
  pu <- pv[ut]; m <- length(pu); ps <- sort(pu); crit <- seq_len(m) * q / m
  k <- which(ps <= crit); bh <- if (length(k)) ps[max(k)] else 0
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

#' @keywords internal
.heterogeneity_pvalues <- function(Xres, groups, W, E0, baseline, B, seed) {
  set.seed(seed); ut <- upper.tri(W); idx <- which(E0 & ut, arr.ind = TRUE)
  if (!nrow(idx)) return(NULL)
  gl <- levels(groups); n <- nrow(Xres)
  obsP <- function(g) {
    Wg <- lapply(split(seq_len(n), g), function(ix) .centered_W(Xres[ix, , drop = FALSE], baseline)$W)
    Reduce(pmax, lapply(Wg, function(w) abs(w - W)))
  }
  P0 <- obsP(groups); cnt <- numeric(nrow(idx))
  for (b in seq_len(B)) {
    Pb <- obsP(groups[sample(n)])
    cnt <- cnt + (Pb[idx] >= P0[idx])
  }
  data.frame(i = idx[, 1], j = idx[, 2], P = P0[idx], p_het = (1 + cnt) / (B + 1))
}
