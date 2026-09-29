# rerun_impact_bootstrap_clean.R -- IMPACT bootstrap stability on the cleaned data (200 resamples, B = 10,000).
# Overnight on a laptop. Run from this folder:  nohup caffeinate -i Rscript rerun_impact_bootstrap_clean.R > log_impact_boot_clean.txt 2>&1 &
suppressPackageStartupMessages({library(biasawareLID); library(lavaan); library(Matrix)})
d <- read.csv("~/Downloads/Impact.csv", check.names = FALSE); d <- d[, nzchar(trimws(names(d)))]
names(d) <- paste0("X", names(d)); Xi <- as.matrix(d); Xi[Xi > 5] <- NA   # seven out-of-range responses on the 1-5 scale
Z <- rowSums(Xi, na.rm = TRUE); c12 <- quantile(Z, c(1/3, 2/3)); gi <- cut(Z, c(-Inf, c12[1], c12[2], Inf), labels = c("Low","Medium","High"))
suf <- sub("^X\\d+", "", colnames(Xi)); smap <- split(colnames(Xi), factor(suf, levels = c("TS","SK","C","PG","CT")))
impact_model <- paste(vapply(names(smap), function(s) paste0(s, " =~ ", paste(smap[[s]], collapse = " + ")), ""), collapse = "\n")
dir.create("impact_boot_clean", showWarnings = FALSE)
bi <- lid_bootstrap(Xi, gi, model = impact_model, reps = 200, B = 10000, seed = 880000)
write.csv(bi$edge_stability, "impact_boot_clean/impact_boot_edge_stability.csv", row.names = FALSE)
write.csv(bi$replications,  "impact_boot_clean/impact_boot_replications.csv",  row.names = FALSE)
cat("IMPACT (cleaned) bootstrap done: stable >=80%:", sum(bi$edge_stability$stable_80), "| >=90%:", sum(bi$edge_stability$stable_90), "\n")
utils::zip("impact_boot_clean.zip", list.files("impact_boot_clean", full.names = TRUE)); cat("ALL DONE -> upload impact_boot_clean.zip\n")
