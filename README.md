# biasawareLID

Bias-aware detection of group-specific local item dependence (Ozkan, 2026).
Centered residual correlations, weighted topological overlap, within-group
permutation inference with BH-FDR control, a qualified reference set E0, a
retention-constrained group-disparity penalty, and signed reversal diagnostics.

## Version 0.2.0 (post-review engine)
- `centering = "model"` (default): each residual correlation is centered by its model-implied
  null baseline, cov2cor(Theta - Lambda(Lambda'Theta^-1 Lambda)^-1 Lambda'), so the null is exact
  for any loading pattern (mean-centering is exact only for equal loadings).
- `null = "parametric"`: reference distribution simulated from the fitted measurement model.
- `heterogeneity_test = TRUE`: group-label permutation test of each qualified pair's disparity
  P_ij with BH control; `$heterogeneity` lists P, p_het, and het_BH.
- `simulate_lid_data(loading = <vector>)` plants heterogeneous loadings.
- Attainability: with m pairs the BH cutoff .05/m must exceed 1/(B+1); use B >= 20,000 for 37-item scales.

## Install
```r
remotes::install_github("fatiihozkann/bias-awareLID")
```

## Quick start
```r
library(biasawareLID)
sim <- simulate_lid_data(N = 500, p = 20, rho = .3, n_groups = 3)
fit <- bias_aware_lid(sim$data, sim$group, B = 5000)   # 190 pairs: B must exceed 190/.05
print(fit); summary(fit)
boot <- lid_bootstrap(sim$data, sim$group, reps = 200)   # stability
```
For real data, pass your lavaan measurement model via `model=`.
All defaults (B = 1000, q = .05, tau = .10, floor = .90, the beta grid) match
the published pipeline; the seed makes every result exactly reproducible.
