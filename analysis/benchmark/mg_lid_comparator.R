# ============================================================================
# Multigroup Limited-Information CFA Comparator for Group-Specific Local
# Dependence (LID)
# ----------------------------------------------------------------------------
# Reviewer-requested benchmark for the bias-aware wTO method (journal version).
#
# WHAT IT DOES
#   Fits a CONFIGURAL multigroup CFA with the WLSMV estimator (limited-
#   information; uses POLYCHORIC correlations for ordinal/dichotomous items),
#   then estimates residual (local) dependence between item pairs WITHIN each
#   group from the Bollen correlation residuals:
#         residual_ij(g) = observed polychoric cor_ij(g) - model-implied cor_ij(g)
#   A large positive residual = the factor model under-predicts the pair's
#   association in that group = local dependence (the Yen's Q3 analog).
#
#   A pair flagged as locally dependent in SOME groups but not ALL groups is
#   classified "group-specific LID" -- the same estimand the bias-aware
#   penalty targets, reached by a non-network, factor-based route.
#
#   Two detection routes are reported:
#     (1) Bollen correlation residuals      (magnitude + z-test, primary)
#     (2) Modification indices (op == "~~")  (classical specification search)
#
# FAIR COMPARISON
#   Per-group Benjamini-Hochberg FDR is applied to the residual z-tests, so the
#   comparator gets the SAME inferential footing (multiplicity control) as the
#   bias-aware method. After that, the remaining contribution of the bias-aware
#   method is the EQUITY PENALTY and the data-tuned beta -- not merely FDR.
#
# CAVEAT (state in the manuscript)
#   This comparator -- like the bias-aware penalty's per-group wTO term -- relies
#   on per-group estimates, so it shares the small-cell instability problem
#   (e.g., PISA-gifted: 6 gender x SES cells, N = 570). The bias-aware method's
#   DETECTION borrows strength from the pooled network; its PENALTY does not.
#
# Requires: lavaan (>= 0.6-12)
# ============================================================================

suppressPackageStartupMessages(library(lavaan))

# ----------------------------------------------------------------------------
# Core comparator
# ----------------------------------------------------------------------------
mg_lid_comparator <- function(data,
                              items,
                              group,
                              ordered     = TRUE,        # TRUE = all items ordered (WLSMV/polychoric); or a vector of ordered item names
                              estimator   = "WLSMV",     # use "MLR" for continuous items
                              model       = NULL,        # optional custom lavaan measurement model; default = single factor
                              factor_name = "F",
                              res_cutoffs = c(.10, .20), # |residual cor| flags: .10 parallels your tau floor; .20 = Yen's Q3 convention
                              z_alpha     = .05,         # alpha for residual z-tests
                              fdr         = TRUE,        # BH-FDR on residual z p-values, within each group
                              mi_cutoff   = 3.84,        # MI flag (chi-square_1, p < .05)
                              classify_by = c("mag20", "mag10", "z", "z_fdr"),
                              fit_pooled  = TRUE,        # also fit pooled model -> residual-cor disparity (parallels P_ij)
                              std.lv      = TRUE,
                              verbose     = TRUE) {

  classify_by <- match.arg(classify_by)
  stopifnot(length(items) >= 3, group %in% names(data))

  # --- single-factor model (edit lhs if your construct is multidimensional) ---
  model   <- if (!is.null(model)) model else paste0(factor_name, " =~ ", paste(items, collapse = " + "))
  ord_arg <- if (isTRUE(ordered)) items else ordered

  # --- configural multigroup fit (loadings/thresholds free per group) --------
  fit <- tryCatch(
    lavaan::cfa(model, data = data, group = group,
                ordered = ord_arg, estimator = estimator, std.lv = std.lv),
    error = function(e)
      stop("Multigroup CFA failed: ", conditionMessage(e),
           "\n  -> With small/sparse cells, collapse groups or categories,",
           "\n     or use fit_groups_separately() below.", call. = FALSE))
  if (!lavInspect(fit, "converged"))
    warning("Multigroup model did not converge cleanly; interpret with caution.")

  G       <- lavInspect(fit, "ngroups")
  glabels <- lavInspect(fit, "group.label")
  gn      <- lavInspect(fit, "nobs"); names(gn) <- glabels

  fitmeas <- tryCatch(
    fitMeasures(fit, c("cfi.scaled", "rmsea.scaled", "srmr")),
    error = function(e) fitMeasures(fit, c("cfi", "rmsea", "srmr")))

  # --- Bollen correlation residuals per group -------------------------------
  res <- lavResiduals(fit, type = "cor.bollen")
  get_res <- function(g) {
    if (G == 1) list(cov = res$cov,      z = res$cov.z)
    else        list(cov = res[[g]]$cov, z = res[[g]]$cov.z)
  }

  # --- modification indices for residual covariances ------------------------
  mi_all <- tryCatch(modindices(fit, op = "~~"), error = function(e) NULL)
  if (!is.null(mi_all)) {
    mi_all <- mi_all[mi_all$lhs %in% items & mi_all$rhs %in% items &
                       mi_all$lhs != mi_all$rhs, , drop = FALSE]
    if ("group" %in% names(mi_all)) mi_all$group_label <- glabels[mi_all$group]
  }

  # --- pooled residual correlations (for disparity, parallels P_ij) ---------
  rescor_pooled <- NULL
  if (isTRUE(fit_pooled)) {
    fitp <- tryCatch(
      lavaan::cfa(model, data = data, ordered = ord_arg,
                  estimator = estimator, std.lv = std.lv),
      error = function(e) NULL)
    if (!is.null(fitp))
      rescor_pooled <- lavResiduals(fitp, type = "cor.bollen")$cov
  }

  # --- assemble long (pair x group) table -----------------------------------
  prs  <- utils::combn(items, 2)
  rows <- vector("list", ncol(prs) * G); r <- 0L
  for (g in seq_len(G)) {
    rg <- get_res(g); M <- rg$cov; Z <- rg$z
    for (k in seq_len(ncol(prs))) {
      i <- prs[1, k]; j <- prs[2, k]
      rc <- M[i, j]
      zz <- if (!is.null(Z)) Z[i, j] else NA_real_
      pp <- if (is.na(zz)) NA_real_ else 2 * stats::pnorm(-abs(zz))
      mi_ij <- NA_real_; epc_ij <- NA_real_
      if (!is.null(mi_all)) {
        grp_match <- if ("group" %in% names(mi_all)) mi_all$group == g else TRUE
        sel <- which(((mi_all$lhs == i & mi_all$rhs == j) |
                      (mi_all$lhs == j & mi_all$rhs == i)) & grp_match)
        if (length(sel)) { mi_ij <- mi_all$mi[sel[1]]; epc_ij <- mi_all$epc[sel[1]] }
      }
      r <- r + 1L
      rows[[r]] <- data.frame(
        item_i = i, item_j = j, pair = paste(i, j, sep = "--"),
        group = glabels[g], n_group = unname(gn[g]),
        rescor = rc, z = zz, p = pp, mi = mi_ij, epc = epc_ij,
        stringsAsFactors = FALSE)
    }
  }
  by_group <- do.call(rbind, rows)

  # --- per-group BH-FDR on residual z p-values ------------------------------
  by_group$p_fdr <- NA_real_
  if (isTRUE(fdr)) {
    for (g in glabels) {
      idx <- which(by_group$group == g & !is.na(by_group$p))
      if (length(idx))
        by_group$p_fdr[idx] <- stats::p.adjust(by_group$p[idx], method = "BH")
    }
  }

  # --- flags ----------------------------------------------------------------
  by_group$flag_mag10 <- abs(by_group$rescor) > res_cutoffs[1]
  by_group$flag_mag20 <- abs(by_group$rescor) > res_cutoffs[length(res_cutoffs)]
  by_group$flag_z     <- !is.na(by_group$p)     & by_group$p     < z_alpha
  by_group$flag_z_fdr <- !is.na(by_group$p_fdr) & by_group$p_fdr < z_alpha
  by_group$flag_mi    <- !is.na(by_group$mi)    & by_group$mi    > mi_cutoff

  # --- collapse to pair level + classify ------------------------------------
  flagcol <- switch(classify_by, mag10 = "flag_mag10", mag20 = "flag_mag20",
                    z = "flag_z", z_fdr = "flag_z_fdr")
  classify <- function(nflag, G)
    if (nflag == 0) "none" else if (nflag >= G) "consistent" else "group-specific"

  pair_ids <- unique(by_group$pair)
  pr <- vector("list", length(pair_ids)); m <- 0L
  for (pp_ in pair_ids) {
    sub <- by_group[by_group$pair == pp_, , drop = FALSE]
    i <- sub$item_i[1]; j <- sub$item_j[1]
    disp <- if (is.null(rescor_pooled)) NA_real_
            else max(abs(sub$rescor - rescor_pooled[i, j]), na.rm = TRUE)
    m <- m + 1L
    pr[[m]] <- data.frame(
      pair = pp_, item_i = i, item_j = j,
      rescor_pooled    = if (is.null(rescor_pooled)) NA_real_ else rescor_pooled[i, j],
      max_group_rescor = sub$rescor[which.max(abs(sub$rescor))],
      disparity        = disp,
      n_flag_mag10 = sum(sub$flag_mag10), n_flag_mag20 = sum(sub$flag_mag20),
      n_flag_z = sum(sub$flag_z),         n_flag_zfdr = sum(sub$flag_z_fdr),
      class_mag10 = classify(sum(sub$flag_mag10), G),
      class_mag20 = classify(sum(sub$flag_mag20), G),
      class_z     = classify(sum(sub$flag_z),     G),
      class_zfdr  = classify(sum(sub$flag_z_fdr), G),
      primary_class = classify(sum(sub[[flagcol]]), G),
      stringsAsFactors = FALSE)
  }
  by_pair <- do.call(rbind, pr)
  by_pair <- by_pair[order(-abs(by_pair$max_group_rescor)), ]

  summ <- as.data.frame(table(factor(by_pair$primary_class,
                    levels = c("none", "group-specific", "consistent"))))
  names(summ) <- c("classification", "n_pairs")

  if (isTRUE(verbose)) {
    cat("Multigroup limited-information CFA comparator\n")
    cat("  Estimator:", estimator, "| Groups:", G,
        "(", paste(glabels, collapse = ", "), ")\n")
    cat("  Per-group n:", paste0(glabels, "=", gn, collapse = ", "), "\n")
    cat("  Configural fit:",
        paste0(names(fitmeas), "=", round(as.numeric(fitmeas), 3), collapse = ", "), "\n")
    cat("  Classification (primary =", classify_by, "):\n")
    print(summ, row.names = FALSE)
  }

  list(fit = fit, n_groups = G, group_labels = glabels, n_per_group = gn,
       fit_measures = fitmeas,
       residuals_by_group = setNames(lapply(seq_len(G), get_res), glabels),
       by_group = by_group, by_pair = by_pair, summary = summ,
       rescor_pooled = rescor_pooled, mi_table = mi_all,
       settings = list(res_cutoffs = res_cutoffs, z_alpha = z_alpha, fdr = fdr,
                       mi_cutoff = mi_cutoff, classify_by = classify_by))
}

# ----------------------------------------------------------------------------
# Agreement with the bias-aware method
#   ba_group_specific_pairs: character vector of "itemA--itemB" pairs the
#   bias-aware method flagged as group-specific (i.e., penalized / removed).
#   Pair strings are order-insensitive (canonicalized internally).
# ----------------------------------------------------------------------------
compare_to_bias_aware <- function(comparator, ba_group_specific_pairs,
                                 class_rule = "primary_class") {
  canon <- function(x) vapply(strsplit(x, "--|__|:| "), function(p)
    paste(sort(p[nzchar(p)]), collapse = "--"), character(1))
  comp_gs <- canon(comparator$by_pair$pair[comparator$by_pair[[class_rule]] == "group-specific"])
  ba_gs   <- canon(ba_group_specific_pairs)

  both    <- intersect(comp_gs, ba_gs)
  only_ba <- setdiff(ba_gs, comp_gs)
  only_cf <- setdiff(comp_gs, ba_gs)
  uni     <- union(comp_gs, ba_gs)
  jacc    <- if (length(uni)) length(both) / length(uni) else NA_real_

  cat("Group-specific LID: bias-aware vs multigroup-CFA comparator\n")
  cat("  Bias-aware flagged :", length(ba_gs), "\n")
  cat("  Comparator flagged :", length(comp_gs), "\n")
  cat("  Agreement (both)   :", length(both), "\n")
  cat("  Only bias-aware    :", length(only_ba),
      if (length(only_ba)) paste0(" [", paste(only_ba, collapse = ", "), "]") else "", "\n")
  cat("  Only comparator    :", length(only_cf),
      if (length(only_cf)) paste0(" [", paste(only_cf, collapse = ", "), "]") else "", "\n")
  cat("  Jaccard overlap    :", round(jacc, 3), "\n")
  invisible(list(both = both, only_bias_aware = only_ba,
                 only_comparator = only_cf, jaccard = jacc))
}

# ----------------------------------------------------------------------------
# Simulation scoring: comparator detection vs known truth
#   true_pairs            : character vector of TRUE LID pairs ("i--j")
#   true_group_specific   : subset that is group-specific in the truth
#   Returns detection rate / precision / F1 for any-LID and group-specific-LID.
# ----------------------------------------------------------------------------
score_comparator_against_truth <- function(comparator, true_pairs,
                                          true_group_specific = character(0),
                                          class_rule = "primary_class") {
  canon <- function(x) vapply(strsplit(x, "--|__|:| "), function(p)
    paste(sort(p[nzchar(p)]), collapse = "--"), character(1))
  bp <- comparator$by_pair
  detected_any <- canon(bp$pair[bp[[class_rule]] != "none"])
  detected_gs  <- canon(bp$pair[bp[[class_rule]] == "group-specific"])
  T_any <- canon(true_pairs); T_gs <- canon(true_group_specific)

  prf <- function(detected, truth) {
    tp <- length(intersect(detected, truth))
    fp <- length(setdiff(detected, truth))
    fn <- length(setdiff(truth, detected))
    prec <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
    rec  <- if (tp + fn > 0) tp / (tp + fn) else NA_real_      # detection rate
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0)
              2 * prec * rec / (prec + rec) else NA_real_
    c(detection_rate = rec, precision = prec, F1 = f1, TP = tp, FP = fp, FN = fn)
  }
  data.frame(target = c("any_LID", "group_specific_LID"),
             rbind(prf(detected_any, T_any), prf(detected_gs, T_gs)),
             row.names = NULL)
}

# ----------------------------------------------------------------------------
# Fallback for small/sparse cells: fit the factor model in each group
# separately (mathematically the configural model, but robust to one cell
# failing). Returns a named list of per-group residual correlation matrices.
# ----------------------------------------------------------------------------
fit_groups_separately <- function(data, items, group, ordered = TRUE,
                                  estimator = "WLSMV", factor_name = "F",
                                  std.lv = TRUE) {
  model   <- paste0(factor_name, " =~ ", paste(items, collapse = " + "))
  ord_arg <- if (isTRUE(ordered)) items else ordered
  gl <- unique(as.character(data[[group]]))
  setNames(lapply(gl, function(g) {
    d <- data[as.character(data[[group]]) == g, , drop = FALSE]
    f <- tryCatch(lavaan::cfa(model, data = d, ordered = ord_arg,
                              estimator = estimator, std.lv = std.lv),
                  error = function(e) {
                    warning("Group ", g, " failed: ", conditionMessage(e)); NULL })
    if (is.null(f)) NULL else lavResiduals(f, type = "cor.bollen")$cov
  }), gl)
}

# ============================================================================
# SELF-CONTAINED RUNNABLE DEMO  (synthetic ordinal data with a planted
# group-specific local dependence). Run with:  RUN_MG_DEMO=1 Rscript this.R
# or just source() interactively. Delete before integrating.
# ============================================================================
if (identical(Sys.getenv("RUN_MG_DEMO"), "1")) {
  set.seed(2026)
  sim_one <- function(n, ld_pair = NULL, ld = 0) {
    p <- 8; theta <- rnorm(n); lam <- rep(.6, p)
    Z <- sapply(seq_len(p), function(i) lam[i] * theta + sqrt(1 - lam[i]^2) * rnorm(n))
    if (!is.null(ld_pair) && ld != 0) {                 # shared residual factor
      e <- rnorm(n)
      Z[, ld_pair[1]] <- Z[, ld_pair[1]] + ld * e
      Z[, ld_pair[2]] <- Z[, ld_pair[2]] + ld * e
    }
    apply(Z, 2, function(z)                              # 4-category ordinal
      as.integer(cut(z, breaks = c(-Inf, -1, 0, 1, Inf), labels = 1:4)))
  }
  A <- sim_one(400, ld_pair = c(3, 4), ld = 1.2)        # group A: LID on i3,i4
  B <- sim_one(400, ld_pair = NULL,    ld = 0)          # group B: none
  dat <- as.data.frame(rbind(A, B)); names(dat) <- paste0("i", 1:8)
  dat$grp <- rep(c("A", "B"), each = 400)

  out <- mg_lid_comparator(dat, items = paste0("i", 1:8), group = "grp",
                           ordered = TRUE, estimator = "WLSMV",
                           classify_by = "mag20")
  cat("\nTop pairs by |max group residual cor| (expect i3--i4 group-specific):\n")
  print(head(out$by_pair[, c("pair", "max_group_rescor", "disparity",
                             "primary_class")], 5), row.names = FALSE)

  sc <- score_comparator_against_truth(out, true_pairs = "i3--i4",
                                       true_group_specific = "i3--i4")
  cat("\nComparator vs truth:\n"); print(sc, row.names = FALSE)
}

# ============================================================================
# EMPIRICAL APPLICATION TEMPLATE  (TIMSS / PISA-gifted / IMPACT)
# ============================================================================
# timss_items <- c("MSC1", "MSC2", "MSC3", "MSC4", "MSC5", "MSC6", "MSC7", "MSC8", "MSC9")
# timss <- mg_lid_comparator(
#   data        = timss_df,
#   items       = timss_items,
#   group       = "gender_ses",      # your 6-level gender x SES factor
#   ordered     = TRUE,              # Likert -> polychoric / WLSMV (limited-information)
#   estimator   = "WLSMV",
#   res_cutoffs = c(.10, .20),
#   classify_by = "mag20")
# print(timss$summary)
# write.csv(timss$by_pair, "timss_mgcfa_comparator_bypair.csv", row.names = FALSE)
#
# # Head-to-head with the edges your bias-aware penalty REMOVED in TIMSS:
# ba_removed <- c("MSC1--MSC4", "MSC1--MSC7", "MSC3--MSC8")   # <- your 3 removed edges
# compare_to_bias_aware(timss, ba_removed, class_rule = "primary_class")
#
# PISA-gifted (N = 570, 6 cells): if the configural fit will not converge,
#   res_by_group <- fit_groups_separately(pisa_df, pisa_items, "gender_ses")
#   and build the pair x group table from res_by_group the same way.
#
# ============================================================================
# SIMULATION INTEGRATION (drop into your single_simulation() in Appendix A)
# ----------------------------------------------------------------------------
# Inside each replication, after you generate `sim_data` (with column `group`)
# and know the planted truth, add the comparator alongside BA and Standard UVA:
#
#   cmp <- mg_lid_comparator(sim_data, items = item_cols, group = "group",
#                            ordered = (item_type != "continuous"),
#                            estimator = if (item_type == "continuous") "MLR" else "WLSMV",
#                            classify_by = "mag20", verbose = FALSE)
#   cmp_score <- score_comparator_against_truth(
#                  cmp,
#                  true_pairs          = true_lid_pairs,          # from your generator
#                  true_group_specific = true_group_specific_pairs)
#
# Then bind cmp_score$detection_rate / precision / F1 into your results table
# next to the bias-aware and UVA columns. (Pass me your generator's truth
# encoding and I'll wire `true_lid_pairs` / `true_group_specific_pairs` exactly.)
# ============================================================================
