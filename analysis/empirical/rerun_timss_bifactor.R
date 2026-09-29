# rerun_timss_bifactor.R -- TIMSS under a bifactor model (general confidence + orthogonal negative-wording method factor).
suppressPackageStartupMessages({library(biasawareLID); library(haven); library(dplyr); library(lavaan); library(Matrix)})
source("timss_bootstrap_prep.R"); dir.create("rerun_v2", showWarnings = FALSE)
fitline <- function(fit, label) { fm <- fitMeasures(fit, c("cfi.robust","tli.robust","rmsea.robust","srmr"))
  cat(sprintf("%-28s CFI %.3f | TLI %.3f | RMSEA %.3f | SRMR %.3f\n", label, fm[1], fm[2], fm[3], fm[4])) }
save_out <- function(x, name) {
  write.csv(x$edges, file.path("rerun_v2", paste0(name, "_edges.csv")), row.names = FALSE)
  if (!is.null(x$heterogeneity)) write.csv(x$heterogeneity, file.path("rerun_v2", paste0(name, "_heterogeneity.csv")), row.names = FALSE)
  sink(file.path("rerun_v2", paste0(name, "_summary.txt"))); print(x); sink(); cat("\n==", name, "==\n"); print(x) }
t <- prepare_timss_bootstrap_data(); Xt <- as.matrix(t$data[, t$items]); gt <- t$groups
mb <- "Gen =~ BSBM19A + BSBM19B + BSBM19C + BSBM19D + BSBM19E + BSBM19F + BSBM19G + BSBM19H + BSBM19I
Neg =~ BSBM19B + BSBM19C + BSBM19E + BSBM19H + BSBM19I
Gen ~~ 0*Neg"
fb <- cfa(mb, data = as.data.frame(Xt), estimator = "MLR", missing = "fiml", std.lv = TRUE)
fitline(fb, "TIMSS bifactor")
save_out(bias_aware_lid(Xt, gt, model = mb, B = 2000, seed = 42, null = "parametric"), "timss_bifactor_parametric")
save_out(bias_aware_lid(Xt, gt, model = mb, B = 2000, seed = 42, heterogeneity_test = TRUE), "timss_bifactor_permutation")
utils::zip("rerun_v2_timss_bifactor.zip", list.files("rerun_v2", pattern = "bifactor", full.names = TRUE))
cat("\nDONE -> paste the fit line and upload rerun_v2_timss_bifactor.zip\n")
