# impact_cleaned_analysis.R -- IMPACT baseline on the cleaned data (out-of-range responses set to missing),
# permutation and parametric references, heterogeneity test, UVA and centered Q3* counts, sensitivity inputs,
# and the multigroup CFA comparator. Needs Impact.csv in ~/Downloads and the comparator script alongside.
suppressPackageStartupMessages({library(biasawareLID); library(lavaan); library(Matrix); library(qgraph)})
source("../benchmark/mg_lid_comparator.R"); source("../benchmark/sim_core_centered.R")
d <- read.csv("~/Downloads/Impact.csv", check.names = FALSE); d <- d[, nzchar(trimws(names(d)))]
names(d) <- paste0("X", names(d)); X <- as.matrix(d)
cat("out-of-range responses set to missing:", sum(X > 5, na.rm = TRUE), "\n"); X[X > 5] <- NA
Z <- rowSums(X, na.rm = TRUE); c12 <- quantile(Z, c(1/3, 2/3)); g <- cut(Z, c(-Inf, c12[1], c12[2], Inf), labels = c("Low","Medium","High"))
suf <- sub("^X\\d+", "", colnames(X)); smap <- split(colnames(X), factor(suf, levels = c("TS","SK","C","PG","CT")))
model <- paste(vapply(names(smap), function(s) paste0(s, " =~ ", paste(smap[[s]], collapse = " + ")), ""), collapse = "\n")
dir.create("impact_clean", showWarnings = FALSE)
perm <- bias_aware_lid(X, g, model = model, B = 20000, seed = 42, heterogeneity_test = TRUE); print(perm)
write.csv(perm$edges, "impact_clean/impact_permutation_edges.csv", row.names = FALSE)
write.csv(perm$heterogeneity, "impact_clean/impact_heterogeneity.csv", row.names = FALSE)
par <- bias_aware_lid(X, g, model = model, B = 20000, seed = 42, null = "parametric"); print(par)
write.csv(par$edges, "impact_clean/impact_parametric_edges.csv", row.names = FALSE)
W <- uva_wto(X); ut <- upper.tri(W); cat("UVA .25:", sum(W[ut] > .25), "| UVA .20:", sum(W[ut] > .20), "\n")
df <- as.data.frame(X); df$grp <- g
cmp <- mg_lid_comparator(df, colnames(X), "grp", ordered = FALSE, estimator = "MLR", model = model, res_cutoffs = c(.10, .20), classify_by = "mag20", verbose = FALSE)
bp <- cmp$by_pair; G <- cmp$n_groups
cat("comparator group-specific pairs: .20 =", sum(bp$n_flag_mag20 > 0 & bp$n_flag_mag20 < G), "| .10 =", sum(bp$n_flag_mag10 > 0 & bp$n_flag_mag10 < G), "\n")
write.csv(bp, "impact_clean/impact_comparator_bypair.csv", row.names = FALSE)
