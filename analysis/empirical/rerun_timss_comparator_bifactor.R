# rerun_timss_comparator_bifactor.R -- MG-CFA comparator for TIMSS under the bifactor model (matches the final analysis).
suppressPackageStartupMessages({library(lavaan); library(haven); library(dplyr)})
source("timss_bootstrap_prep.R"); source("mg_lid_comparator.R"); dir.create("rerun_v2", showWarnings = FALSE)
t <- prepare_timss_bootstrap_data(); dat <- as.data.frame(t$data); dat$Group <- t$groups
mb <- "Gen =~ BSBM19A + BSBM19B + BSBM19C + BSBM19D + BSBM19E + BSBM19F + BSBM19G + BSBM19H + BSBM19I
Neg =~ BSBM19B + BSBM19C + BSBM19E + BSBM19H + BSBM19I
Gen ~~ 0*Neg"
comp <- mg_lid_comparator(data = dat, items = t$items, group = "Group", ordered = TRUE, estimator = "WLSMV",
                          model = mb, res_cutoffs = c(.10, .20), z_alpha = .05, fdr = TRUE, classify_by = "mag20",
                          fit_pooled = TRUE, std.lv = TRUE, verbose = FALSE)
bp <- comp$by_pair; write.csv(bp, "rerun_v2/timss_comparator_bifactor_pairs.csv", row.names = FALSE)
gs20 <- sum(bp$n_flag_mag20 > 0 & bp$n_flag_mag20 < comp$n_groups); gs10 <- sum(bp$n_flag_mag10 > 0 & bp$n_flag_mag10 < comp$n_groups)
cat(sprintf("TIMSS comparator under bifactor: group-specific pairs at .20 = %d of 36; at .10 = %d of 36\n", gs20, gs10))
cat("flagged pairs at .20:", paste(bp$pair[bp$n_flag_mag20 > 0 & bp$n_flag_mag20 < comp$n_groups], collapse = ", "), "\n")
cat("DONE -> paste these lines and upload rerun_v2/timss_comparator_bifactor_pairs.csv\n")
