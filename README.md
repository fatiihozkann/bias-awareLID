# biasawareLID

Bias-aware detection of group-specific local item dependence (Ozkan, 2026).
Centered residual correlations, weighted topological overlap, within-group
permutation inference with BH-FDR control, a qualified reference set E0, a
retention-constrained group-disparity penalty, and signed reversal diagnostics.

## Install
```r
remotes::install_github("fatiihozkann/bias-awareLID")
```

## Quick start
```r
library(biasawareLID)
sim <- simulate_lid_data(N = 500, p = 20, rho = .3, n_groups = 3)
fit <- bias_aware_lid(sim$data, sim$group, B = 1000)
print(fit); summary(fit)
boot <- lid_bootstrap(sim$data, sim$group, reps = 200)   # stability
```
For real data, pass your lavaan measurement model via `model=`.
All defaults (B = 1000, q = .05, tau = .10, floor = .90, the beta grid) match
the published pipeline; the seed makes every result exactly reproducible.
