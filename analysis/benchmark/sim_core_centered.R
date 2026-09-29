library(tidyr)

# ============================================================================
# PART 1: UTILITY FUNCTIONS
# ============================================================================

#' Convert continuous data to dichotomous (median split)
to_dichotomous <- function(x) {
 as.integer(x > stats::median(x, na.rm = TRUE))
}

#' Convert continuous data to 5-point polytomous scale
to_polytomous5 <- function(x) {
 xj <- x + stats::rnorm(length(x), sd = 1e-6) # Break ties
 br <- stats::quantile(xj, probs = seq(0, 1, length.out = 6), na.rm = TRUE)
 br <- unique(br)
 while (length(br) < 6) br <- sort(unique(c(br, br + 1e-6 * seq_along(br))))
 as.integer(cut(xj, breaks = br, include.lowest = TRUE))
}

#' Safe correlation matrix computation
safe_cor <- function(X, method = c("pearson", "spearman")) {
 method <- match.arg(method)
 X <- as.matrix(X)
 sds <- apply(X, 2, sd, na.rm = TRUE)
 const <- which(is.na(sds) | sds == 0)
 if (length(const)) X[, const] <- X[, const] + rnorm(nrow(X) * length(const), sd = 1e-6)
 suppressWarnings(cor(X, use = "pairwise.complete.obs", method = method))
}

#' Make matrix positive semi-definite
make_psd <- function(S, corr = TRUE) {
 S <- (S + t(S)) / 2
 S[is.na(S)] <- 0
 as.matrix(Matrix::nearPD(S, corr = corr)$mat)
}

#' Extract upper triangle values
upper_vals <- function(M) M[upper.tri(M, diag = FALSE)]

# ============================================================================
# PART 2: RESIDUALIZATION (CFA-based)
# ============================================================================

#' Residualize data using CFA factor scores
#' @param X Data matrix
#' @param dtype Data type: "continuous", "polytomous", or "dichotomous"
#' @return Matrix of residualized scores
residualize_from_cfa <- function(X, dtype = c("continuous", "polytomous", "dichotomous")) {
 dtype <- match.arg(dtype)
 p <- ncol(X)
 vn <- colnames(X)
 if (is.null(vn)) vn <- paste0("x", 1:p)
 df <- as.data.frame(X)
 colnames(df) <- vn
 
 # Single-factor CFA model
 model <- paste0("F =~ ", paste(vn, collapse = " + "))
 est <- if (dtype == "continuous") "MLR" else "WLSMV"
 ord <- if (dtype == "continuous") NULL else vn
 
 # Try CFA first
 fit <- try(suppressWarnings(lavaan::cfa(model, data = df, std.lv = TRUE,
 estimator = est, ordered = ord)), silent = TRUE)
 
 if (!inherits(fit, "try-error")) {
 methods_try <- if (dtype == "continuous") c("Bartlett", "regression") else c("regression", "Bartlett")
 for (m in methods_try) {
 fs <- try(suppressWarnings(as.numeric(lavaan::lavPredict(fit, method = m))), silent = TRUE)
 if (!inherits(fs, "try-error") && all(is.finite(fs)) && sd(fs) > 0) {
 # Residualize each item
 R <- matrix(NA_real_, nrow(df), p)
 for (j in 1:p) {
 y <- df[[j]]
 ok <- is.finite(fs) & is.finite(y)
 if (sum(ok) > 2 && sd(fs[ok]) > 0) {
 R[ok, j] <- stats::resid(stats::lm(y ~ fs))
 miss <- which(!ok | is.na(R[, j]))
 if (length(miss)) R[miss, j] <- y[miss] - mean(y[ok], na.rm = TRUE)
 } else {
 R[, j] <- y - mean(y, na.rm = TRUE)
 }
 }
 colnames(R) <- vn
 return(R)
 }
 }
 }
 
 # Fallback: PCA
 Z <- scale(df)
 pc <- try(stats::prcomp(Z, center = FALSE, scale. = FALSE), silent = TRUE)
 if (!inherits(pc, "try-error") && is.matrix(pc$x) && sd(pc$x[, 1]) > 0) {
 fs <- as.numeric(pc$x[, 1])
 R <- matrix(NA_real_, nrow(df), p)
 for (j in 1:p) {
 y <- df[[j]]
 ok <- is.finite(fs) & is.finite(y)
 if (sum(ok) > 2) {
 R[ok, j] <- stats::resid(stats::lm(y ~ fs))
 } else {
 R[, j] <- y - mean(y, na.rm = TRUE)
 }
 }
 colnames(R) <- vn
 return(R)
 }
 
 # Last resort: mean centering
 as.matrix(scale(df, center = TRUE, scale = FALSE))
}

# ============================================================================
# PART 3: WEIGHTED TOPOLOGICAL OVERLAP (wTO)
# ============================================================================

#' Compute Q3 matrix (residual correlations after CFA residualization)
#' Used by sensitivity analysis demo scripts in this file
#' @param X Item response data
#' @param dtype Data type for residualization
compute_Q3_simple <- function(X, dtype = "continuous") {
 X <- as.data.frame(X)
 # Auto-detect dtype if not specified
 if (missing(dtype)) {
 n_unique <- length(unique(unlist(X)))
 dtype <- if (n_unique == 2) "dichotomous" else if (n_unique <= 7) "polytomous" else "continuous"
 }
 Xres <- residualize_from_cfa(X, dtype = dtype)
 C <- make_psd(safe_cor(Xres, method = "pearson"), corr = TRUE)
 diag(C) <- 0
 C
}

#' Compute adjacency matrix from Q3* bias-corrected residual correlations
adjacency_from_cor <- function(Xres) {
 C <- make_psd(safe_cor(Xres, method = "pearson"), corr = TRUE)
 # Q3* bias correction (Yen, 1993): subtract mean off-diagonal Q3
 q3_mean <- mean(C[row(C) != col(C)])
 C_star <- C - q3_mean
 A <- abs(C_star)
 diag(A) <- 0
 A
}

#' Compute Weighted Topological Overlap (wTO) matrix
#' @param C_or_A Correlation or adjacency matrix
#' @param is_adj TRUE if input is already adjacency matrix
#' @return wTO matrix
compute_wto <- function(C_or_A, is_adj = FALSE) {
 A <- if (is_adj) C_or_A else {
 q3_mean <- mean(C_or_A[row(C_or_A) != col(C_or_A)])
 C_star <- C_or_A - q3_mean
 A <- abs(C_star)
 diag(A) <- 0
 A
 }
 p <- ncol(A)
 k <- rowSums(A) # Node connectivity
 
 # Vectorized: shared neighbors via matrix multiplication
 shared <- A %*% A # shared[i,j] = sum_k A[i,k]*A[k,j]
 num <- shared + A # Add direct connection
 
 # Min degree matrix: min(k[i], k[j]) for all pairs
 min_k <- outer(k, k, pmin)
 den <- min_k + 1 - A
 
 W <- ifelse(den > 0, num / den, 0)
 diag(W) <- 0
 # Ensure symmetry
 W <- (W + t(W)) / 2
 W
}

# ============================================================================
# PART 4: BIAS PENALTY COMPUTATION
# ============================================================================

#' Compute cross-group bias penalty matrix
#' @param W_all Pooled wTO matrix
#' @param W_groups List of group-specific wTO matrices
#' @return Penalty matrix P where P[i,j] = max deviation across groups
compute_bias_penalty <- function(W_all, W_groups) {
 p <- ncol(W_all)
 P <- matrix(0, p, p)
 
 for (i in 1:(p - 1)) {
 for (j in (i + 1):p) {
 # Maximum absolute deviation from pooled wTO
 devs <- sapply(W_groups, function(Wg) abs(Wg[i, j] - W_all[i, j]))
 P[i, j] <- max(devs)
 P[j, i] <- P[i, j]
 }
 }
 diag(P) <- 0
 P
}

#' Apply bias penalty to wTO matrix
#' @param W wTO matrix
#' @param P Penalty matrix
#' @param beta Penalty weight parameter
#' @return Adjusted wTO* matrix
apply_bias_penalty <- function(W, P, beta = 1) {
 Wb <- W - beta * P
 Wb[Wb < 0] <- 0
 Wb[Wb > 1] <- 1
 Wb
}

# ============================================================================
# PART 5: PERMUTATION TESTING & BH-FDR
# ============================================================================

#' Compute permutation p-values for wTO
#' @param Xres Residualized data matrix
#' @param groups Group membership vector
#' @param wto_obs Observed wTO matrix
#' @param B Number of permutations
#' @param seed Random seed
#' @return Matrix of p-values
permute_wto_pvals <- function(Xres, groups, wto_obs, B = 1000, seed = NULL) {
 if (!is.null(seed)) set.seed(seed)
 p <- ncol(Xres)
 utm <- upper.tri(wto_obs, diag = FALSE)
 obs_vals <- wto_obs[utm]
 null_mat <- matrix(NA_real_, B, sum(utm))
 
 for (b in 1:B) {
 Xb <- Xres
 # Permute within each group
 for (gi in unique(groups)) {
 idx <- which(groups == gi)
 for (j in 1:p) {
 Xb[idx, j] <- sample(Xb[idx, j], length(idx), replace = FALSE)
 }
 }
 Ab <- adjacency_from_cor(Xb)
 Wb <- compute_wto(Ab, is_adj = TRUE)
 null_mat[b, ] <- Wb[utm]
 }
 
 # Two-sided p-value: proportion of |null| >= |observed|
 pvals <- (colSums(abs(null_mat) >= matrix(abs(obs_vals), nrow = B, ncol = length(obs_vals),
 byrow = TRUE), na.rm = TRUE) + 1) / (B + 1)
 P <- matrix(1, p, p)
 P[utm] <- pvals
 P <- t(P)
 P[utm] <- pvals
 diag(P) <- 1
 P
}

#' Benjamini-Hochberg FDR threshold
#' @param pvals_ut Upper triangle p-values
#' @param alpha FDR level
#' @return BH-adjusted threshold
bh_fdr_cut <- function(pvals_ut, alpha = 0.05) {
 p <- pvals_ut[is.finite(pvals_ut)]
 m <- length(p)
 if (m == 0) return(0)
 p_sorted <- sort(p, na.last = NA)
 crit <- (seq_len(m) * alpha) / m
 pass <- which(p_sorted <= crit)
 if (length(pass) == 0) return(0)
 p_sorted[max(pass)]
}

# ============================================================================
# PART 6: GRID SEARCH WITH POWER FLOOR (KEY METHODOLOGICAL CONTRIBUTION)
# ============================================================================

#' Select optimal beta using grid search with power floor
#' @param W_full Pooled wTO matrix
#' @param Pmat Penalty matrix
#' @param sel_bh BH-FDR selection matrix (edges meeting significance threshold)
#' @param beta_grid Vector of candidate beta values
#' @param power_floor Minimum proportion of BH edges to retain (default 0.90)
#' @param wto_floor Minimum wTO* effect size for an edge to be retained (default 0.10)
#' @return List with optimal beta and selection matrix
grid_search_beta <- function(W_full, Pmat, sel_bh, 
 beta_grid = c(0, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 5),
 power_floor = 0.90,
 wto_floor = 0.10) {
 
 # Count BH-FDR significant edges (reference set)
 ut <- upper.tri(sel_bh, diag = FALSE)
 n_bh_edges <- sum(sel_bh[ut])
 
 if (n_bh_edges == 0) {
 return(list(
 beta_star = 0,
 edges_retained = 0,
 overlap_pct = NA,
 sel_final = sel_bh
 ))
 }
 
 results <- data.frame(
 beta = beta_grid,
 edges = NA_integer_,
 overlap_pct = NA_real_,
 mean_bias = NA_real_,
 meets_floor = NA
 )
 
 for (i in seq_along(beta_grid)) {
 b <- beta_grid[i]
 
 # Apply penalty
 W_pen <- apply_bias_penalty(W_full, Pmat, beta = b)
 
 # Count edges with wTO* > wto_floor among BH-significant edges (dual criterion)
 sel_pen <- (W_pen > wto_floor) & sel_bh
 n_retained <- sum(sel_pen[ut])
 overlap <- n_retained / n_bh_edges
 
 # Mean bias among retained edges
 mean_bias <- if (n_retained > 0) mean(Pmat[sel_pen & ut]) else NA
 
 results$edges[i] <- n_retained
 results$overlap_pct[i] <- round(overlap * 100, 1)
 results$mean_bias[i] <- round(mean_bias, 4)
 results$meets_floor[i] <- (overlap >= power_floor) | (n_retained >= n_bh_edges - 1)
 }
 
 # Select highest beta meeting power floor
 valid <- which(results$meets_floor)
 if (length(valid) == 0) {
 beta_star <- beta_grid[1] # Default to no penalty
 } else {
 beta_star <- max(beta_grid[valid])
 }
 
 # Final selection with optimal beta (dual criterion: significant AND above wTO floor)
 W_final <- apply_bias_penalty(W_full, Pmat, beta = beta_star)
 sel_final <- (W_final > wto_floor) & sel_bh
 
 list(
 beta_star = beta_star,
 grid_results = results,
 edges_retained = sum(sel_final[ut]),
 overlap_pct = results$overlap_pct[which(beta_grid == beta_star)],
 sel_final = sel_final,
 W_final = W_final
 )
}

# ============================================================================
# PART 7: DATA GENERATION
# ============================================================================

#' Generate safe (disjoint) local dependency pairs
safe_dep_pairs <- function(p, k, disjoint = TRUE) {
 k <- max(1, k)
 if (disjoint) {
 k_max <- floor(p / 2)
 k <- min(k, k_max)
 perm <- sample.int(p, p, replace = FALSE)
 pairs <- split(perm[1:(2 * k)], rep(1:k, each = 2))
 lapply(pairs, function(v) sort(as.integer(v)))
 } else {
 replicate(k, sort(sample.int(p, 2)), simplify = FALSE)
 }
}

#' Simulate data for a single group with local dependencies
#' @param n Sample size
#' @param p Number of items
#' @param loadings_vec Factor loadings
#' @param resid_var_vec Residual variances
#' @param dep_pairs List of item pairs with local dependency
#' @param dep_rho Correlation for local dependency
#' @param extra_dep Additional dependency (for group-specific DIF)
simulate_group <- function(n, p, loadings_vec, resid_var_vec,
 dep_pairs = list(c(1, 2)), dep_rho = 0.3, extra_dep = 0.0) {
 lambda <- loadings_vec
 Psi <- diag(resid_var_vec)
 
 # Add local dependencies to residual covariance
 for (pair in dep_pairs) {
 i <- pair[1]
 j <- pair[2]
 rho <- min(0.95, dep_rho + extra_dep)
 Psi[i, j] <- rho
 Psi[j, i] <- rho
 }
 diag(Psi) <- 1
 Psi <- make_psd(Psi, corr = TRUE)
 
 # Generate factor scores and errors
 eta <- rnorm(n, 0, 1)
 E <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Psi)
 
 # Observed scores = factor contribution + residuals
 Y <- outer(eta, lambda) + E
 as.data.frame(Y)
}

# ============================================================================
# PART 8: STANDARD UVA (Christensen et al., 2023 comparison)
# ============================================================================

#' Standard UVA using EBICglasso
uva_wto <- function(X) {
 X <- as.matrix(X)
 S <- safe_cor(X, method = "pearson")
 G <- try(qgraph::EBICglasso(S, n = nrow(X)), silent = TRUE)
 if (inherits(G, "try-error")) {
 fit <- BGGM::estimate(S, n = nrow(X), analytic = TRUE)
 G <- BGGM::pcor_mat(fit)
 }
 A <- abs(G)
 diag(A) <- 0
 compute_wto(A, is_adj = TRUE)
}

# ============================================================================
# PART 9: SINGLE SIMULATION RUN
# ============================================================================

#' Run a single simulation condition
#' @param N Total sample size
#' @param p Number of items
#' @param dep_rho Local dependency effect size
#' @param n_groups Number of demographic groups
#' @param dtype Data type
#' @param beta_grid Grid of beta values for search
#' @param power_floor Power floor for beta selection
#' @param B Number of permutations
#' @param alpha FDR significance level
#' @param seed Random seed
run_single_simulation <- function(
 N = 500, p = 20, dep_rho = 0.3,
 n_groups = 2, group_props = NULL,
 dtype = c("continuous", "polytomous", "dichotomous"),
 beta_grid = c(0, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 5),
 power_floor = 0.90,
 wto_floor = 0.10,
 B = 1000, alpha = 0.05, seed = NULL
) {
 dtype <- match.arg(dtype)
 if (!is.null(seed)) set.seed(seed)
 
 # Number of planted dependency pairs
 k_pairs <- max(1, floor(p / 5))
 dep_pairs <- safe_dep_pairs(p, k_pairs, disjoint = TRUE)
 
 # Group sizes
 if (is.null(group_props)) group_props <- rep(1 / n_groups, n_groups)
 group_props <- group_props / sum(group_props)
 group_sizes <- as.integer(round(N * group_props))
 
 # Adjust to match N exactly
 while (sum(group_sizes) < N) {
 group_sizes[which.max(group_props)] <- group_sizes[which.max(group_props)] + 1
 }
 while (sum(group_sizes) > N) {
 group_sizes[which.max(group_sizes)] <- group_sizes[which.max(group_sizes)] - 1
 }
 
 # Simulate data with group-specific DIF
 loadings_vec <- rep(0.7, p)
 resid_var_vec <- rep(1, p)
  # Designed truth classes (Decision 2, locked 2026-08-12): half consistent, half
  # group-specific. Group-specific pairs carry dep_rho ONLY in one randomly chosen
  # affected group (manuscript definition: present in one group, ~zero elsewhere).
  n_dep <- length(dep_pairs)
  n_gs_planted <- floor(n_dep / 2)
  gs_sel <- if (n_gs_planted > 0) sample(seq_len(n_dep), n_gs_planted) else integer(0)
  gs_pairs_true <- dep_pairs[gs_sel]
  cs_pairs_true <- dep_pairs[setdiff(seq_len(n_dep), gs_sel)]
  affected_group <- sample.int(n_groups, 1)

  dfl <- vector("list", n_groups)
  for (gi in seq_len(n_groups)) {
    pairs_gi <- c(cs_pairs_true, if (gi == affected_group) gs_pairs_true else list())
    dfl[[gi]] <- simulate_group(group_sizes[gi], p, loadings_vec, resid_var_vec,
                                dep_pairs = pairs_gi, dep_rho = dep_rho,
                                extra_dep = 0)
    colnames(dfl[[gi]]) <- paste0("X", 1:p)
    dfl[[gi]]$group <- factor(gi)
  }
 df <- do.call(rbind, dfl)
 items <- paste0("X", 1:p)
 
 # Apply discretization if needed
 if (dtype == "polytomous") {
 for (it in items) df[[it]] <- to_polytomous5(df[[it]])
 }
 if (dtype == "dichotomous") {
 for (it in items) df[[it]] <- to_dichotomous(df[[it]])
 }
 
 # ==================== BIAS-AWARE METHOD ====================
 
 # Step 1: Residualize (pooled)
 Xres_full <- residualize_from_cfa(df[, items, drop = FALSE], dtype = dtype)
 
 # Step 2: Compute pooled wTO
 A_full <- adjacency_from_cor(Xres_full)
 W_full <- compute_wto(A_full, is_adj = TRUE)
 
 # Step 3: Compute group-specific wTO
 W_groups <- lapply(seq_len(n_groups), function(gi) {
 Xg <- df[df$group == gi, items, drop = FALSE]
 Xrg <- residualize_from_cfa(Xg, dtype = dtype)
 Ag <- adjacency_from_cor(Xrg)
 compute_wto(Ag, is_adj = TRUE)
 })
 
 # Step 4: Compute bias penalty
 Pmat <- compute_bias_penalty(W_full, W_groups)
 
 # Step 5: Permutation p-values
 G <- as.integer(df$group)
 Pvals <- permute_wto_pvals(Xres_full, G, W_full, B = B, seed = sample.int(1e7, 1))
 
 # Step 6: BH-FDR threshold
 pv_ut <- Pvals[upper.tri(Pvals, diag = FALSE)]
 pv_cut <- bh_fdr_cut(pv_ut, alpha = alpha)
 sel_bh <- (Pvals <= pv_cut)
 
 # Step 7: Grid search with power floor
 grid_result <- grid_search_beta(W_full, Pmat, sel_bh, 
 beta_grid = beta_grid, 
 power_floor = power_floor,
 wto_floor = wto_floor)
 
 # ==================== STANDARD UVA (Comparison) ====================
 W_uva <- uva_wto(df[, items, drop = FALSE])
 sel_uva_025 <- (W_uva > 0.25) # Fixed threshold (strict)
 sel_uva_020 <- (W_uva > 0.20) # Fixed threshold (liberal, for sensitivity)
 
 # ==================== Q3 BENCHMARK (Yen 1984; Christensen et al., 2017) ====================
 # A_full = |Q3*| (already computed above, no extra cost)
 # Cutoffs span Christensen et al. (2017) critical-value range
 sel_q3_020 <- A_full > 0.20 # Yen (1984) traditional cutoff
 sel_q3_015 <- A_full > 0.15 # Christensen et al. (2017) liberal
 sel_q3_025 <- A_full > 0.25 # Christensen et al. (2017) strict
 diag(sel_q3_020) <- FALSE
 diag(sel_q3_015) <- FALSE
 diag(sel_q3_025) <- FALSE
 
 # ==================== EVALUATE DETECTION ====================
 
 # True positive = correctly identified planted dependency
 evaluate_detection <- function(sel_mat, true_pairs) {
 detected <- 0
 for (pair in true_pairs) {
 if (sel_mat[pair[1], pair[2]] || sel_mat[pair[2], pair[1]]) {
 detected <- detected + 1
 }
 }
 detected
 }
 
 ut <- upper.tri(sel_bh, diag = FALSE)
 
  # Truth labels are assigned at generation (designed classes), not post hoc.
  gs_pairs <- gs_pairs_true
  cs_pairs <- cs_pairs_true
 
 # Overall detection counts
 TP_biasaware <- evaluate_detection(grid_result$sel_final, dep_pairs)
 TP_standard <- evaluate_detection(sel_uva_025, dep_pairs)
 TP_standard_020 <- evaluate_detection(sel_uva_020, dep_pairs)
 TP_q3_020 <- evaluate_detection(sel_q3_020, dep_pairs)
 TP_q3_015 <- evaluate_detection(sel_q3_015, dep_pairs)
 TP_q3_025 <- evaluate_detection(sel_q3_025, dep_pairs)
 
 # Group-specific detection
 TP_gs_biasaware <- evaluate_detection(grid_result$sel_final, gs_pairs)
 TP_gs_standard <- evaluate_detection(sel_uva_025, gs_pairs)
 TP_gs_standard_020 <- evaluate_detection(sel_uva_020, gs_pairs)
 TP_gs_q3_020 <- evaluate_detection(sel_q3_020, gs_pairs)
 TP_gs_q3_015 <- evaluate_detection(sel_q3_015, gs_pairs)
 TP_gs_q3_025 <- evaluate_detection(sel_q3_025, gs_pairs)
 
 # Consistent detection
 TP_cs_biasaware <- evaluate_detection(grid_result$sel_final, cs_pairs)
 TP_cs_standard <- evaluate_detection(sel_uva_025, cs_pairs)
 TP_cs_standard_020 <- evaluate_detection(sel_uva_020, cs_pairs)
 TP_cs_q3_020 <- evaluate_detection(sel_q3_020, cs_pairs)
 TP_cs_q3_015 <- evaluate_detection(sel_q3_015, cs_pairs)
 TP_cs_q3_025 <- evaluate_detection(sel_q3_025, cs_pairs)
 
 # False positive counts
 dep_set <- matrix(FALSE, p, p)
 for (pair in dep_pairs) {
 dep_set[pair[1], pair[2]] <- dep_set[pair[2], pair[1]] <- TRUE
 }
 
 FP_biasaware <- sum(grid_result$sel_final[ut] & !dep_set[ut])
 FP_standard <- sum(sel_uva_025[ut] & !dep_set[ut])
 FP_standard_020 <- sum(sel_uva_020[ut] & !dep_set[ut])
 FP_q3_020 <- sum(sel_q3_020[ut] & !dep_set[ut])
 FP_q3_015 <- sum(sel_q3_015[ut] & !dep_set[ut])
 FP_q3_025 <- sum(sel_q3_025[ut] & !dep_set[ut])
 
 n_gs <- length(gs_pairs)
 n_cs <- length(cs_pairs)
 
  ed <- getOption("sim_export_dir", NULL)
  if (!is.null(ed)) try(.sim_export_all(ed, environment()), silent = TRUE)

 # Return results
 list(
 # Simulation parameters
 N = N, p = p, dep_rho = dep_rho, n_groups = n_groups, dtype = dtype,
 n_true_pairs = length(dep_pairs),
 n_gs_pairs = n_gs,
 n_cs_pairs = n_cs,
    affected_group = affected_group,
    gs_pairs_true = gs_pairs,
    cs_pairs_true = cs_pairs,
 
 # Beta selection
 beta_star = grid_result$beta_star,
 grid_results = grid_result$grid_results,
 
 # Overall detection results
 TP_biasaware = TP_biasaware,
 TP_standard = TP_standard,
 TP_standard_020 = TP_standard_020,
 TP_q3_020 = TP_q3_020,
 TP_q3_015 = TP_q3_015,
 TP_q3_025 = TP_q3_025,
 FP_biasaware = FP_biasaware,
 FP_standard = FP_standard,
 FP_standard_020 = FP_standard_020,
 FP_q3_020 = FP_q3_020,
 FP_q3_015 = FP_q3_015,
 FP_q3_025 = FP_q3_025,
 
 # Group-specific detection (the key metric for Bias-Aware)
 TP_gs_biasaware = TP_gs_biasaware,
 TP_gs_standard = TP_gs_standard,
 TP_gs_standard_020 = TP_gs_standard_020,
 TP_gs_q3_020 = TP_gs_q3_020,
 TP_gs_q3_015 = TP_gs_q3_015,
 TP_gs_q3_025 = TP_gs_q3_025,
 
 # Consistent detection
 TP_cs_biasaware = TP_cs_biasaware,
 TP_cs_standard = TP_cs_standard,
 TP_cs_standard_020 = TP_cs_standard_020,
 TP_cs_q3_020 = TP_cs_q3_020,
 TP_cs_q3_015 = TP_cs_q3_015,
 TP_cs_q3_025 = TP_cs_q3_025,
 
 # Edge counts
 edges_biasaware = sum(grid_result$sel_final[ut]),
 edges_standard = sum(sel_uva_025[ut]),
 edges_standard_020 = sum(sel_uva_020[ut]),
 edges_q3_020 = sum(sel_q3_020[ut]),
 edges_q3_015 = sum(sel_q3_015[ut]),
 edges_q3_025 = sum(sel_q3_025[ut]),
 edges_bh_fdr = sum(sel_bh[ut]),
 
 # Overall detection rates
 det_rate_biasaware = TP_biasaware / length(dep_pairs),
 det_rate_standard = TP_standard / length(dep_pairs),
 det_rate_standard_020 = TP_standard_020 / length(dep_pairs),
 det_rate_q3_020 = TP_q3_020 / length(dep_pairs),
 det_rate_q3_015 = TP_q3_015 / length(dep_pairs),
 det_rate_q3_025 = TP_q3_025 / length(dep_pairs),
 
 # Group-specific detection rates
 det_rate_gs_biasaware = if (n_gs > 0) TP_gs_biasaware / n_gs else NA,
 det_rate_gs_standard = if (n_gs > 0) TP_gs_standard / n_gs else NA,
 det_rate_gs_standard_020 = if (n_gs > 0) TP_gs_standard_020 / n_gs else NA,
 det_rate_gs_q3_020 = if (n_gs > 0) TP_gs_q3_020 / n_gs else NA,
 det_rate_gs_q3_015 = if (n_gs > 0) TP_gs_q3_015 / n_gs else NA,
 det_rate_gs_q3_025 = if (n_gs > 0) TP_gs_q3_025 / n_gs else NA,
 
 # Consistent detection rates
 det_rate_cs_biasaware = if (n_cs > 0) TP_cs_biasaware / n_cs else NA,
 det_rate_cs_standard = if (n_cs > 0) TP_cs_standard / n_cs else NA,
 det_rate_cs_standard_020 = if (n_cs > 0) TP_cs_standard_020 / n_cs else NA,
 det_rate_cs_q3_020 = if (n_cs > 0) TP_cs_q3_020 / n_cs else NA,
 det_rate_cs_q3_015 = if (n_cs > 0) TP_cs_q3_015 / n_cs else NA,
 det_rate_cs_q3_025 = if (n_cs > 0) TP_cs_q3_025 / n_cs else NA,
 
 # Matrices (for plotting if needed)
 W_full = W_full,
 W_final = grid_result$W_final,
 Pmat = Pmat,
 sel_biasaware = grid_result$sel_final,
 sel_standard = sel_uva_025,
 sel_q3 = sel_q3_020,
 dep_pairs = dep_pairs
 )
}

# ============================================================================
# PART 10: FULL SIMULATION STUDY
# ============================================================================

#' Run full Monte Carlo simulation study
#' @param replications Number of replications per condition
#' @param B Permutations per replication
#' @param alpha FDR level
#' @param power_floor Power floor for beta selection
#' @param verbose Print progress
#' @param parallel Use parallel processing
run_simulation_study <- function(
 replications = 100,
 B = 1000,
 alpha = 0.05,
 power_floor = 0.90,
 wto_floor = 0.10,
 verbose = TRUE,
 parallel = FALSE
) {
 
 # ---- FULL SIMULATION DESIGN (162 conditions) ----
 Ns <- c(250, 500, 1000) # Sample sizes
 Ps <- c(10, 20, 30) # Number of items
 rhos <- c(0.1, 0.3, 0.5) # Effect sizes
 groups <- getOption("sim_groups", c(2, 4))  # FLAG: manuscript says 2,3 -- Fatih to confirm
 types <- c("continuous", "polytomous", "dichotomous") # Data types
 
 # Create condition grid
 conditions <- expand.grid(
 N = Ns, p = Ps, rho = rhos, n_groups = groups, dtype = types,
 stringsAsFactors = FALSE
 )
 n_conditions <- nrow(conditions)
 
 if (verbose) {
 cat("=== SIMULATION STUDY ===\n")
 cat(sprintf("Conditions: %d\n", n_conditions))
 cat(sprintf("Replications: %d\n", replications))
 cat(sprintf("Total runs: %d\n", n_conditions * replications))
 cat(sprintf("Power floor: %.0f%%\n", power_floor * 100))
 cat("========================\n\n")
 }
 
 # Results storage
 all_results <- list()
 
 # Pre-generate random seeds for reproducibility
 rep_seeds <- sample.int(1e7, replications)
 
 # Helper function to run one replication
 run_one_rep <- function(r, cond, rep_seeds, power_floor, wto_floor, B, alpha) {
 run_single_simulation(
 N = cond$N,
 p = cond$p,
 dep_rho = cond$rho,
 n_groups = cond$n_groups,
 dtype = cond$dtype,
 power_floor = power_floor,
 wto_floor = wto_floor,
 B = B,
 alpha = alpha,
 seed = rep_seeds[r]
 )
 }
 
 for (i in 1:n_conditions) {
 cond <- conditions[i, ]
 
 if (verbose) {
 cat(sprintf("\nCondition %d/%d: N=%d, p=%d, rho=%.1f, groups=%d, type=%s\n",
 i, n_conditions, cond$N, cond$p, cond$rho, cond$n_groups, cond$dtype))
 }
 
 # Replications - parallel or sequential based on parallel flag
 if (parallel && .Platform$OS.type != "windows") {
 # Unix/Linux: use mclapply
 rep_results <- parallel::mclapply(1:replications, function(r) {
 run_one_rep(r, cond, rep_seeds, power_floor, wto_floor, B, alpha)
 }, mc.cores = parallel::detectCores() - 1)
 } else if (parallel) {
 # Windows: use parLapply
 cl <- parallel::makeCluster(parallel::detectCores() - 1)
 parallel::clusterExport(cl, c("run_single_simulation", "run_one_rep", "cond", "rep_seeds",
 "power_floor", "wto_floor", "B", "alpha"),
 envir = environment())
 rep_results <- parallel::parLapply(cl, 1:replications, function(r) {
 run_one_rep(r, cond, rep_seeds, power_floor, wto_floor, B, alpha)
 })
 parallel::stopCluster(cl)
 } else {
 # Sequential execution
 rep_results <- lapply(1:replications, function(r) {
 if (verbose && r %% 10 == 0) cat(sprintf(" Rep %d/%d\n", r, replications))
 run_one_rep(r, cond, rep_seeds, power_floor, wto_floor, B, alpha)
 })
 }
 
 # Aggregate results
 all_results[[i]] <- data.frame(
 N = cond$N,
 p = cond$p,
 rho = cond$rho,
 n_groups = cond$n_groups,
 dtype = cond$dtype,
 
 # Detection rates (proportion of replications with at least one TP)
 det_rate_biasaware = mean(sapply(rep_results, function(x) x$det_rate_biasaware > 0)),
 det_rate_standard = mean(sapply(rep_results, function(x) x$det_rate_standard > 0)),
 det_rate_standard_020 = mean(sapply(rep_results, function(x) x$det_rate_standard_020 > 0)),
 
 # Mean TP counts
 mean_TP_biasaware = mean(sapply(rep_results, function(x) x$TP_biasaware)),
 mean_TP_standard = mean(sapply(rep_results, function(x) x$TP_standard)),
 mean_TP_standard_020 = mean(sapply(rep_results, function(x) x$TP_standard_020)),
 
 # Mean FP counts
 mean_FP_biasaware = mean(sapply(rep_results, function(x) x$FP_biasaware)),
 mean_FP_standard = mean(sapply(rep_results, function(x) x$FP_standard)),
 mean_FP_standard_020 = mean(sapply(rep_results, function(x) x$FP_standard_020)),
 
 # Group-specific TP (the key advantage of Bias-Aware)
 mean_TP_gs_biasaware = mean(sapply(rep_results, function(x) x$TP_gs_biasaware)),
 mean_TP_gs_standard = mean(sapply(rep_results, function(x) x$TP_gs_standard)),
 mean_TP_gs_standard_020 = mean(sapply(rep_results, function(x) x$TP_gs_standard_020)),
 
 # Consistent TP
 mean_TP_cs_biasaware = mean(sapply(rep_results, function(x) x$TP_cs_biasaware)),
 mean_TP_cs_standard = mean(sapply(rep_results, function(x) x$TP_cs_standard)),
 mean_TP_cs_standard_020 = mean(sapply(rep_results, function(x) x$TP_cs_standard_020)),
 
 # Group-specific detection rates
 det_rate_gs_biasaware = mean(sapply(rep_results, function(x) {
 r <- x$det_rate_gs_biasaware; if (is.na(r)) 0 else as.numeric(r > 0)
 })),
 det_rate_gs_standard = mean(sapply(rep_results, function(x) {
 r <- x$det_rate_gs_standard; if (is.na(r)) 0 else as.numeric(r > 0)
 })),
 
 # Consistent detection rates
 det_rate_cs_biasaware = mean(sapply(rep_results, function(x) {
 r <- x$det_rate_cs_biasaware; if (is.na(r)) 0 else as.numeric(r > 0)
 })),
 det_rate_cs_standard = mean(sapply(rep_results, function(x) {
 r <- x$det_rate_cs_standard; if (is.na(r)) 0 else as.numeric(r > 0)
 })),
 
 # Mean beta selected
 mean_beta = mean(sapply(rep_results, function(x) x$beta_star)),
 
 # Mean edges
 mean_edges_biasaware = mean(sapply(rep_results, function(x) x$edges_biasaware)),
 mean_edges_standard = mean(sapply(rep_results, function(x) x$edges_standard)),
 mean_edges_standard_020 = mean(sapply(rep_results, function(x) x$edges_standard_020)),
 
 # Standard errors
 se_det_biasaware = sd(sapply(rep_results, function(x) x$det_rate_biasaware > 0)) / sqrt(replications),
 se_det_standard = sd(sapply(rep_results, function(x) x$det_rate_standard > 0)) / sqrt(replications),
 se_det_standard_020 = sd(sapply(rep_results, function(x) x$det_rate_standard_020 > 0)) / sqrt(replications)
 )
 }
 
 # Combine all results
 results_df <- do.call(rbind, all_results)
 
 if (verbose) {
 cat("\n=== SIMULATION COMPLETE ===\n")
 }
 
 results_df
}

# ============================================================================
# PART 11: RESULTS VISUALIZATION
# ============================================================================

#' Plot detection rates by factor
plot_detection_rates <- function(results_df) {
 
 # By effect size
 p1 <- results_df %>%
 group_by(rho) %>%
 summarise(
 `Bias-Aware` = mean(det_rate_biasaware),
 `Standard UVA` = mean(det_rate_standard),
 .groups = "drop"
 ) %>%
 pivot_longer(cols = c(`Bias-Aware`, `Standard UVA`), 
 names_to = "Method", values_to = "Detection") %>%
 ggplot(aes(x = factor(rho), y = Detection, fill = Method)) +
 geom_bar(stat = "identity", position = "dodge") +
 labs(x = "Effect Size (ρ)", y = "Detection Rate", 
 title = "Detection by Effect Size") +
 theme_minimal() +
 scale_fill_manual(values = c("Bias-Aware" = "#2E86AB", "Standard UVA" = "#A23B72"))
 
 # By sample size
 p2 <- results_df %>%
 group_by(N) %>%
 summarise(
 `Bias-Aware` = mean(det_rate_biasaware),
 `Standard UVA` = mean(det_rate_standard),
 .groups = "drop"
 ) %>%
 pivot_longer(cols = c(`Bias-Aware`, `Standard UVA`), 
 names_to = "Method", values_to = "Detection") %>%
 ggplot(aes(x = factor(N), y = Detection, fill = Method)) +
 geom_bar(stat = "identity", position = "dodge") +
 labs(x = "Sample Size (N)", y = "Detection Rate", 
 title = "Detection by Sample Size") +
 theme_minimal() +
 scale_fill_manual(values = c("Bias-Aware" = "#2E86AB", "Standard UVA" = "#A23B72"))
 
 # By number of items
 p3 <- results_df %>%
 group_by(p) %>%
 summarise(
 `Bias-Aware` = mean(det_rate_biasaware),
 `Standard UVA` = mean(det_rate_standard),
 .groups = "drop"
 ) %>%
 pivot_longer(cols = c(`Bias-Aware`, `Standard UVA`), 
 names_to = "Method", values_to = "Detection") %>%
 ggplot(aes(x = factor(p), y = Detection, fill = Method)) +
 geom_bar(stat = "identity", position = "dodge") +
 labs(x = "Number of Items (p)", y = "Detection Rate", 
 title = "Detection by Number of Items") +
 theme_minimal() +
 scale_fill_manual(values = c("Bias-Aware" = "#2E86AB", "Standard UVA" = "#A23B72"))
 
 list(by_effect = p1, by_sample = p2, by_items = p3)
}

#' Create heatmap of detection rates
plot_detection_heatmap <- function(results_df) {
 results_df %>%
 group_by(N, rho) %>%
 summarise(Detection = mean(det_rate_biasaware), .groups = "drop") %>%
 ggplot(aes(x = factor(rho), y = factor(N), fill = Detection)) +
 geom_tile() +
 geom_text(aes(label = sprintf("%.2f", Detection)), color = "white", size = 4) +
 scale_fill_gradient(low = "#F0F0F0", high = "#2E86AB") +
 labs(x = "Effect Size (ρ)", y = "Sample Size (N)", 
 title = "Bias-Aware Detection Rate Heatmap") +
 theme_minimal()
}

# ============================================================================
# PART 12: EXAMPLE USAGE
# ============================================================================

if (FALSE) { # Quick test - set to TRUE if you just want to verify code works
 
 # Run small test (quick, ~5-10 minutes)
 cat("\n>>> Running quick test (5 reps, B=50)...\n\n")
 test_results <- run_simulation_study(
 replications = 5,
 B = 50,
 power_floor = 0.90,
 verbose = TRUE
 )
 print(test_results)
 write.csv(test_results, "test_results.csv", row.names = FALSE)
 cat("\n>>> Test complete! Results saved to test_results.csv\n")
}

if (FALSE) { # FULL SIMULATION - 162 conditions × 100 reps × B=1000
 
 cat("\n>>> Running FULL simulation on HPC...\n")
 cat(">>> 162 conditions × 100 replications × B=1,000 permutations\n")
 cat(">>> Total iterations: 16,200,000\n\n")
 
 full_results <- run_simulation_study(
 replications = 100,
 B = 1000,
 power_floor = 0.90,
 wto_floor = 0.10,
 parallel = TRUE,
 verbose = TRUE
 )
 
 # Save results
 write.csv(full_results, "simulation_results.csv", row.names = FALSE)
 
 # Create plots
 plots <- plot_detection_rates(full_results)
 ggsave("Figure_4_1_detection_by_effect.png", plots$by_effect, width = 8, height = 6)
 ggsave("Figure_4_2_detection_by_sample.png", plots$by_sample, width = 8, height = 6)
 
 heatmap <- plot_detection_heatmap(full_results)
 ggsave("Figure_4_3_heatmap.png", heatmap, width = 8, height = 6)
}

# ============================================================================
# END OF SIMULATION CODE
