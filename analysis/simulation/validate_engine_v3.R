# validate_engine_v3.R -- the validity studies the review requires (engine v0.3.0). Run from ~/sim_centered:
#   install.packages("biasawareLID_0.3.0.tar.gz", repos = NULL, type = "source")
#   nohup caffeinate -i Rscript validate_engine_v3.R > log_validate_v3.txt 2>&1 &
# Part A: equivalence of mean vs model-implied centering on IDENTICAL datasets (equal population loadings).
# Part B: null calibration under model-implied centering with equal and heterogeneous loadings, both references.
suppressPackageStartupMessages({library(biasawareLID); library(lavaan); library(Matrix)})
dir.create("validate_v3", showWarnings = FALSE)
ap <- function(path, df) { h <- !file.exists(path); suppressWarnings(write.table(df, path, sep = ",", row.names = FALSE, col.names = h, append = !h)) }
key <- function(e, col) paste(e$item_i[e[[col]]], e$item_j[e[[col]]], sep = "-")
jac <- function(a, b) { u <- union(a, b); if (!length(u)) 1 else length(intersect(a, b)) / length(u) }
# ---- Part A ----
gridA <- expand.grid(N = c(250, 500, 1000), p = c(10, 20, 30), rho = c(0, .3), stringsAsFactors = FALSE)
for (r in seq_len(nrow(gridA))) for (rep in 1:50) {
  cd <- gridA[r, ]; s <- simulate_lid_data(N = cd$N, p = cd$p, rho = cd$rho, n_groups = 2, loading = .7, seed = 10000 * r + rep)
  m <- bias_aware_lid(s$data, s$group, B = 2000, seed = 500 + rep, centering = "mean")
  g <- bias_aware_lid(s$data, s$group, B = 2000, seed = 500 + rep, centering = "model")
  ap("validate_v3/partA_equivalence.csv", data.frame(N = cd$N, p = cd$p, rho = cd$rho, rep = rep,
     bh_mean = m$n_bh, bh_model = g$n_bh, e0_mean = m$n_E0, e0_model = g$n_E0, final_mean = m$n_final, final_model = g$n_final,
     jaccard_E0 = jac(key(m$edges, "E0"), key(g$edges, "E0")), jaccard_final = jac(key(m$edges, "final"), key(g$edges, "final")),
     max_abs_Qdiff = max(abs(m$edges$Q3_star_signed - g$edges$Q3_star_signed))))
  cat(sprintf("[A] N=%d p=%d rho=%.1f rep %d: final mean/model %d/%d\n", cd$N, cd$p, cd$rho, rep, m$n_final, g$n_final))
}
# ---- Part B ----
gridB <- expand.grid(N = c(250, 500, 1000), p = c(10, 20), loadings = c("equal", "spread_30_90", "spread_50_80"), stringsAsFactors = FALSE)
lamf <- function(kind, p) switch(kind, equal = rep(.7, p), spread_30_90 = seq(.3, .9, length.out = p), spread_50_80 = seq(.5, .8, length.out = p))
for (r in seq_len(nrow(gridB))) for (rep in 1:200) {
  cd <- gridB[r, ]; s <- simulate_lid_data(N = cd$N, p = cd$p, rho = 0, n_groups = 2, loading = lamf(cd$loadings, cd$p), seed = 900000 + 1000 * r + rep)
  pm <- bias_aware_lid(s$data, s$group, B = 2000, seed = 700 + rep, centering = "model", null = "permutation")
  pp <- bias_aware_lid(s$data, s$group, B = 2000, seed = 700 + rep, centering = "model", null = "parametric")
  ap("validate_v3/partB_null.csv", data.frame(N = cd$N, p = cd$p, loadings = cd$loadings, rep = rep,
     bh_perm = pm$n_bh, e0_perm = pm$n_E0, final_perm = pm$n_final, bh_param = pp$n_bh, e0_param = pp$n_E0, final_param = pp$n_final,
     max_wTO = max(pm$edges$wTO)))
  if (rep %% 25 == 0) cat(sprintf("[B] N=%d p=%d %s rep %d\n", cd$N, cd$p, cd$loadings, rep))
}
a <- read.csv("validate_v3/partA_equivalence.csv"); b <- read.csv("validate_v3/partB_null.csv")
sa <- aggregate(cbind(jaccard_E0, jaccard_final, max_abs_Qdiff) ~ N + p + rho, a, mean)
sb <- aggregate(cbind(any_bh_perm = bh_perm > 0, any_final_perm = final_perm > 0, any_bh_param = bh_param > 0, any_final_param = final_param > 0) ~ N + p + loadings, b, mean)
write.csv(sa, "validate_v3/summary_partA.csv", row.names = FALSE); write.csv(sb, "validate_v3/summary_partB.csv", row.names = FALSE)
print(sa); print(sb)
utils::zip("validate_v3_results.zip", list.files("validate_v3", full.names = TRUE)); cat("\nALL DONE -> upload validate_v3_results.zip\n")
