# rerun_bootstraps_v2.R -- PISA and IMPACT bootstrap stability under engine v0.2.0 (model centering).
# PISA: 200 reps x B = 2,000 (~1-2 h). IMPACT: 200 reps x B = 10,000 (attainable with >= 2 tied pairs; ~overnight).
suppressPackageStartupMessages({library(biasawareLID); library(haven); library(dplyr); library(lavaan); library(Matrix)})
source("pisa_bootstrap_prep.R"); dir.create("rerun_v2", showWarnings = FALSE)
p <- prepare_pisa_bootstrap_data(); Xp <- as.matrix(p$data[, p$items]); gp <- p$groups
bp <- lid_bootstrap(Xp, gp, reps = 200, B = 2000, seed = 880000)
write.csv(bp$edge_stability, "rerun_v2/pisa_boot_edge_stability.csv", row.names = FALSE)
write.csv(bp$replications,  "rerun_v2/pisa_boot_replications.csv",  row.names = FALSE)
cat("PISA bootstrap done: stable >=80%:", sum(bp$edge_stability$stable_80), "| >=90%:", sum(bp$edge_stability$stable_90), "\n")
d <- read.csv("~/Downloads/Impact.csv", check.names = FALSE); d <- d[, nzchar(trimws(names(d)))]; names(d) <- paste0("X", names(d)); Xi <- as.matrix(d); Xi[Xi > 5] <- NA   # seven out-of-range responses (6, 33, 43, 44) on the 1-5 scale
Z <- rowSums(Xi, na.rm = TRUE); c12 <- quantile(Z, c(1/3, 2/3)); gi <- cut(Z, c(-Inf, c12[1], c12[2], Inf), labels = c("Low","Medium","High"))
suf <- sub("^X\\d+", "", names(d)); smap <- split(names(d), factor(suf, levels = c("TS","SK","C","PG","CT")))
impact_model <- paste(vapply(names(smap), function(s) paste0(s, " =~ ", paste(smap[[s]], collapse = " + ")), ""), collapse = "\n")
bi <- lid_bootstrap(Xi, gi, model = impact_model, reps = 200, B = 10000, seed = 880000)
write.csv(bi$edge_stability, "rerun_v2/impact_boot_edge_stability.csv", row.names = FALSE)
write.csv(bi$replications,  "rerun_v2/impact_boot_replications.csv",  row.names = FALSE)
cat("IMPACT bootstrap done: stable >=80%:", sum(bi$edge_stability$stable_80), "| >=90%:", sum(bi$edge_stability$stable_90), "\n")
utils::zip("rerun_v2_bootstraps.zip", list.files("rerun_v2", pattern = "boot", full.names = TRUE))
cat("\nALL DONE -> upload rerun_v2_bootstraps.zip\n")
