# rerun_benchmark_v3.R -- targeted benchmark under engine v0.3.1 (model-implied centering, parametric reference,
# direct-dependence stage, heterogeneity test) vs the MG-CFA comparator rules. 3 formats x 100 reps.
# Run from ~/sim_centered (needs benchmark_generator.R and mg_lid_comparator.R from the August benchmark folder):
#   Rscript -e 'install.packages("biasawareLID_0.3.1.tar.gz", repos=NULL, type="source")'
#   nohup caffeinate -i Rscript rerun_benchmark_v3.R > log_bench_v3.txt 2>&1 &
suppressPackageStartupMessages({library(biasawareLID); library(lavaan); library(Matrix)})
source("sim_core_centered.R"); source("benchmark_helpers.R"); source("benchmark_generator.R"); source("mg_lid_comparator.R")
dir.create("bench_v3", showWarnings = FALSE); unlink(list.files("bench_v3", full.names = TRUE))
ap <- function(path, df){ h<-!file.exists(path); suppressWarnings(write.table(df,path,sep=",",row.names=FALSE,col.names=h,append=!h)) }
canon <- function(i, j) paste(pmin(i, j), pmax(i, j), sep = "--")
for (dtype in c("continuous", "polytomous", "dichotomous")) for (rep in 1:100) {
  seed <- 7000L + rep; gen <- generate_targeted_benchmark_data(dtype = dtype, seed = seed)
  X <- as.matrix(gen$data[, gen$items]); grp <- factor(gen$data$group); aff <- as.character(gen$affected_group)
  it <- gen$items
  cs <- vapply(gen$consistent_pairs, function(p) canon(it[p[1]], it[p[2]]), ""); gs <- vapply(gen$group_specific_pairs, function(p) canon(it[p[1]], it[p[2]]), "")
  planted <- c(cs, gs)
  res <- tryCatch(bias_aware_lid(X, grp, B = 2000, seed = seed, centering = "model", null = "parametric",
                  direct_test = TRUE, heterogeneity_test = TRUE), error = function(e) NULL)
  if (is.null(res)) { cat("FAIL", dtype, rep, "\n", file = "bench_v3/failures.log", append = TRUE); next }
  e <- res$edges; key <- canon(e$item_i, e$item_j)
  fin <- key[e$final]; e0 <- key[e$E0]
  tp <- sum(fin %in% planted); fp <- sum(!fin %in% planted)
  # affected-group localization for planted group-specific pairs in E0: argmax_g |W^(g) - W|
  loc_hits <- 0L; loc_n <- 0L
  for (pr in gen$group_specific_pairs) { i <- pr[1]; j <- pr[2]
    if (canon(it[i], it[j]) %in% e0) { loc_n <- loc_n + 1L
      dev <- vapply(res$W_groups, function(w) abs(w[i, j] - res$W[i, j]), 0)
      if (names(which.max(dev)) == aff) loc_hits <- loc_hits + 1L } }
  het <- res$heterogeneity; hk <- if (!is.null(het)) canon(het$item_i, het$item_j) else character(0)
  het_gs <- if (length(hk)) sum(het$het_BH & hk %in% gs) else 0; het_cs_false <- if (length(hk)) sum(het$het_BH & hk %in% cs) else 0
  row <- data.frame(dtype = dtype, rep = rep, affected_group = aff, n_bh = res$n_bh, n_E0 = res$n_E0, n_final = res$n_final, beta = res$beta_star,
    recall_any = tp/4, precision_any = ifelse(tp+fp > 0, tp/(tp+fp), NA), F1_any = ifelse(tp+fp > 0, 2*tp/(2*tp+fp+(4-tp)), NA),
    recall_gs_final = sum(fin %in% gs)/2, gs_in_E0 = loc_n, localization_hits = loc_hits, het_flags_gs = het_gs, het_false_cs = het_cs_false,
    direct_confirmed = sum(e$direct_BH, na.rm = TRUE), removed = res$n_E0 - res$n_final)
  # MG-CFA comparator rules (one call; all four rules recovered from the flag counts)
  comp <- tryCatch(mg_lid_comparator(data = gen$data, items = gen$items, group = "group",
            ordered = (dtype != "continuous"), estimator = if (dtype == "continuous") "MLR" else "WLSMV",
            res_cutoffs = c(.10, .20), z_alpha = .05, fdr = TRUE, classify_by = "mag20", fit_pooled = TRUE, std.lv = TRUE, verbose = FALSE),
            error = function(e) NULL)
  if (!is.null(comp)) { bp <- comp$by_pair; pk <- canon(bp$item_i, bp$item_j)
    for (rule in c("mag10", "mag20", "z", "zfdr")) { nf <- bp[[paste0("n_flag_", rule)]]; det <- pk[nf > 0]
      tpc <- sum(det %in% planted); fpc <- sum(!det %in% planted)
      row[[paste0("cmp_", rule, "_recall")]] <- tpc/4; row[[paste0("cmp_", rule, "_precision")]] <- ifelse(tpc+fpc > 0, tpc/(tpc+fpc), NA)
      # affected-group flag for planted gs pairs
      bg <- comp$by_group; fcol <- c(mag10 = "flag_mag10", mag20 = "flag_mag20", z = "flag_z", zfdr = "flag_z_fdr")[rule]
      hits <- 0L
      for (pr in gen$group_specific_pairs) { sub <- bg[canon(bg$item_i, bg$item_j) == canon(it[pr[1]], it[pr[2]]) & as.character(bg$group) == aff, ]
        if (nrow(sub) && any(sub[[fcol]])) hits <- hits + 1L }
      row[[paste0("cmp_", rule, "_affected_recall")]] <- hits/2 } }
  ap("bench_v3/replication_level.csv", row)
  if (rep %% 10 == 0) cat(sprintf("[%s] %s rep %d done\n", format(Sys.time(), "%H:%M:%S"), dtype, rep))
}
d <- read.csv("bench_v3/replication_level.csv")
s <- aggregate(. ~ dtype, d[, c("dtype", grep("recall|precision|F1|localization|gs_in_E0|het_|direct|removed", names(d), value = TRUE))], mean, na.rm = TRUE, na.action = na.pass)
write.csv(s, "bench_v3/summary.csv", row.names = FALSE); print(t(s))
utils::zip("bench_v3_results.zip", list.files("bench_v3", full.names = TRUE)); cat("\nALL DONE -> upload bench_v3_results.zip\n")
