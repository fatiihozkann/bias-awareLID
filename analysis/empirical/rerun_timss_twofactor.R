# rerun_timss_twofactor.R -- TIMSS under a two-factor wording model (positive vs negative items),
# plus fit indices for all three measurement models. Run from ~/sim_centered after rerun_empirical_v2.R.
suppressPackageStartupMessages({library(biasawareLID); library(haven); library(dplyr); library(lavaan); library(Matrix)})
source("timss_bootstrap_prep.R"); source("pisa_bootstrap_prep.R")
dir.create("rerun_v2", showWarnings = FALSE)
fitline <- function(fit, label) { fm <- fitMeasures(fit, c("cfi.robust","tli.robust","rmsea.robust","srmr"))
  cat(sprintf("%-28s CFI %.3f | TLI %.3f | RMSEA %.3f | SRMR %.3f\n", label, fm[1], fm[2], fm[3], fm[4])) }
save_out <- function(x, name) {
  write.csv(x$edges, file.path("rerun_v2", paste0(name, "_edges.csv")), row.names = FALSE)
  if (!is.null(x$heterogeneity)) write.csv(x$heterogeneity, file.path("rerun_v2", paste0(name, "_heterogeneity.csv")), row.names = FALSE)
  sink(file.path("rerun_v2", paste0(name, "_summary.txt"))); print(x); sink(); cat("\n==", name, "==\n"); print(x) }
# ---- TIMSS: one-factor vs two-factor wording model ----
t <- prepare_timss_bootstrap_data(); Xt <- as.matrix(t$data[, t$items]); gt <- t$groups
m1 <- paste0("F1 =~ ", paste(t$items, collapse = " + "))
m2 <- "Pos =~ BSBM19A + BSBM19D + BSBM19F + BSBM19G
Neg =~ BSBM19B + BSBM19C + BSBM19E + BSBM19H + BSBM19I"
fitline(cfa(m1, data = as.data.frame(Xt), estimator = "MLR", missing = "fiml"), "TIMSS one-factor")
fitline(cfa(m2, data = as.data.frame(Xt), estimator = "MLR", missing = "fiml"), "TIMSS two-factor wording")
save_out(bias_aware_lid(Xt, gt, model = m2, B = 2000, seed = 42, heterogeneity_test = TRUE), "timss_twofactor_permutation")
save_out(bias_aware_lid(Xt, gt, model = m2, B = 2000, seed = 42, null = "parametric"),      "timss_twofactor_parametric")
# ---- fit indices for the PISA and IMPACT models used in rerun_empirical_v2.R ----
p <- prepare_pisa_bootstrap_data(); Xp <- as.matrix(p$data[, p$items])
fitline(cfa(paste0("F1 =~ ", paste(p$items, collapse = " + ")), data = as.data.frame(Xp), estimator = "MLR", missing = "fiml"), "PISA one-factor")
d <- read.csv("~/Downloads/Impact.csv", check.names = FALSE); d <- d[, nzchar(trimws(names(d)))]; names(d) <- paste0("X", names(d)); d[d > 5] <- NA   # seven out-of-range responses on the 1-5 scale
suf <- sub("^X\\d+", "", names(d)); smap <- split(names(d), factor(suf, levels = c("TS","SK","C","PG","CT")))
impact_model <- paste(vapply(names(smap), function(s) paste0(s, " =~ ", paste(smap[[s]], collapse = " + ")), ""), collapse = "\n")
fitline(cfa(impact_model, data = d, estimator = "MLR", missing = "fiml"), "IMPACT five-factor")
utils::zip("rerun_v2_timss2f.zip", list.files("rerun_v2", pattern = "twofactor", full.names = TRUE))
cat("\nDONE -> paste the four fit-index lines here and upload rerun_v2_timss2f.zip\n")
