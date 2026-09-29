# rerun_timss_bifactor_bootstrap.R -- run AFTER rerun_bootstraps_v2.R finishes (~1 h).
suppressPackageStartupMessages({library(biasawareLID); library(haven); library(dplyr); library(lavaan); library(Matrix)})
source("timss_bootstrap_prep.R"); dir.create("rerun_v2", showWarnings = FALSE)
t <- prepare_timss_bootstrap_data(); Xt <- as.matrix(t$data[, t$items]); gt <- t$groups
mb <- "Gen =~ BSBM19A + BSBM19B + BSBM19C + BSBM19D + BSBM19E + BSBM19F + BSBM19G + BSBM19H + BSBM19I
Neg =~ BSBM19B + BSBM19C + BSBM19E + BSBM19H + BSBM19I
Gen ~~ 0*Neg"
bt <- lid_bootstrap(Xt, gt, model = mb, reps = 200, B = 2000, seed = 880000)
write.csv(bt$edge_stability, "rerun_v2/timss_boot_edge_stability.csv", row.names = FALSE)
write.csv(bt$replications,  "rerun_v2/timss_boot_replications.csv",  row.names = FALSE)
cat("TIMSS bifactor bootstrap done: stable >=80%:", sum(bt$edge_stability$stable_80), "| >=90%:", sum(bt$edge_stability$stable_90), "\n")
utils::zip("rerun_v2_timss_boot.zip", list.files("rerun_v2", pattern = "timss_boot", full.names = TRUE))
cat("DONE -> upload rerun_v2_timss_boot.zip\n")
