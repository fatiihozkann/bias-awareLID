# rerun_empirical_v2.R -- ONE script: installs engine v0.2.0 and reruns the three empirical baselines
# under model-implied centering, with the heterogeneity test and a parametric-null robustness check.
# Run from ~/sim_centered (needs biasawareLID_0.2.0.tar.gz, timss_bootstrap_prep.R, pisa_bootstrap_prep.R,
# and your data files at the paths those prep scripts already use; Impact.csv in ~/Downloads).
if (!requireNamespace("biasawareLID", quietly = TRUE) ||
    packageVersion("biasawareLID") < "0.2.0")
  install.packages("biasawareLID_0.2.0.tar.gz", repos = NULL, type = "source")
suppressPackageStartupMessages({library(biasawareLID); library(haven); library(dplyr); library(lavaan); library(Matrix)})
source("timss_bootstrap_prep.R"); source("pisa_bootstrap_prep.R")
dir.create("rerun_v2", showWarnings = FALSE)
save_out <- function(x, name) {
  write.csv(x$edges, file.path("rerun_v2", paste0(name, "_edges.csv")), row.names = FALSE)
  if (!is.null(x$heterogeneity))
    write.csv(x$heterogeneity, file.path("rerun_v2", paste0(name, "_heterogeneity.csv")), row.names = FALSE)
  sink(file.path("rerun_v2", paste0(name, "_summary.txt"))); print(x); sink(); cat("\n==", name, "==\n"); print(x)
}
# ---- TIMSS (36 pairs; B = 2,000 attainable) ----
t <- prepare_timss_bootstrap_data(); Xt <- as.matrix(t$data[, t$items]); gt <- t$groups
save_out(bias_aware_lid(Xt, gt, B = 2000, seed = 42, heterogeneity_test = TRUE), "timss_model_permutation")
save_out(bias_aware_lid(Xt, gt, B = 2000, seed = 42, null = "parametric"),      "timss_model_parametric")
# ---- PISA (28 pairs) ----
p <- prepare_pisa_bootstrap_data(); Xp <- as.matrix(p$data[, p$items]); gp <- p$groups
save_out(bias_aware_lid(Xp, gp, B = 2000, seed = 42, heterogeneity_test = TRUE), "pisa_model_permutation")
save_out(bias_aware_lid(Xp, gp, B = 2000, seed = 42, null = "parametric"),      "pisa_model_parametric")
# ---- IMPACT (666 pairs; B = 20,000 needed for attainability; ~20-40 min per run) ----
d <- read.csv("~/Downloads/Impact.csv", check.names = FALSE); d <- d[, nzchar(trimws(names(d)))]
names(d) <- paste0("X", names(d)); Xi <- as.matrix(d); Xi[Xi > 5] <- NA   # seven out-of-range responses (6, 33, 43, 44) on the 1-5 scale
Z <- rowSums(Xi, na.rm = TRUE); c12 <- quantile(Z, c(1/3, 2/3))
gi <- cut(Z, c(-Inf, c12[1], c12[2], Inf), labels = c("Low", "Medium", "High"))
suf <- sub("^X\\d+", "", names(d)); smap <- split(names(d), factor(suf, levels = c("TS", "SK", "C", "PG", "CT")))
impact_model <- paste(vapply(names(smap), function(s) paste0(s, " =~ ", paste(smap[[s]], collapse = " + ")), ""), collapse = "\n")
save_out(bias_aware_lid(Xi, gi, model = impact_model, B = 20000, seed = 42, heterogeneity_test = TRUE), "impact_model_permutation")
save_out(bias_aware_lid(Xi, gi, model = impact_model, B = 20000, seed = 42, null = "parametric"),      "impact_model_parametric")
utils::zip("rerun_v2_results.zip", list.files("rerun_v2", full.names = TRUE))
cat("\nALL DONE -> upload rerun_v2_results.zip\n")
