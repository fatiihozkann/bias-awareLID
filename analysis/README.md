# Analysis scripts

These scripts reproduce every result in the article. Install the package first:

    remotes::install_github("fatiihozkann/bias-awareLID")

Run each script from its own folder. Data files are read from `~/Downloads`:
TIMSS 2019 grade 8 U.S. student file (`bsgusam7.sav`, IEA TIMSS 2019 International Database),
the PISA 2022 student questionnaire file (`CY08MSP_STU_QQQ.SAV`, OECD PISA 2022 Database), and
`Impact.csv` (not public; see the article's data availability statement). In IMPACT, seven
out-of-range responses (values above 5 on the 1-5 scale) are set to missing in every script.

| Folder | Script | Produces |
|---|---|---|
| simulation | `rerun_main_sim_v3.R` | Main 108-condition simulation (Tables 2-3, Figure 2, Figure S2) |
| simulation | `validate_engine_v3.R` | Global-null study and centering-equivalence study (Figure 4) |
| benchmark | `rerun_benchmark_v3.R` | Targeted benchmark (Table 4, Figure 3) |
| empirical | `rerun_timss_twofactor.R`, `rerun_timss_bifactor.R` | TIMSS one-, two-factor, and bifactor analyses; fit indices |
| empirical | `rerun_empirical_v2.R` | PISA baseline (TIMSS one-factor, IMPACT) |
| empirical | `impact_cleaned_analysis.R` | IMPACT baseline, parametric check, UVA, comparator |
| empirical | `rerun_bootstraps_v2.R`, `rerun_timss_bifactor_bootstrap.R`, `rerun_impact_bootstrap_clean.R` | Bootstrap stability (Table 5, Figure 5) |
| empirical | `rerun_timss_comparator_bifactor.R` | TIMSS multigroup CFA comparator (Table 6) |

Seeds are fixed in every script. The simulation scripts are chunked and restartable
(`--chunk=i/k --skip-done`).

## Saved output: IMPACT baseline on the cleaned data

`Impact.csv` cannot be shared, so the pair-level output of `empirical/impact_cleaned_analysis.R`
(B = 20,000, seed = 42, package version 0.3.1) is saved in `empirical/impact_clean/`:

| File | Content |
|---|---|
| `impact_permutation_edges.csv` | All 666 item pairs under the within-group permutation reference: unpenalized wTO, signed centered Q3*, permutation p value, group disparity, and the BH, E0 (qualified), and final (retained) indicators |
| `impact_heterogeneity.csv` | Heterogeneity test for the qualified pairs: disparity P, unadjusted p value, BH decision |
| `sessionInfo.txt` | R and package versions used for this saved run |

These files give the IMPACT row of the article: groups of 340, 335, and 301 after seven
out-of-range responses are set to missing; 8 BH-significant pairs, 4 qualified, 4 retained at
beta* = 0 (X4TS-X22C, X26PG-X27PG, X23PG-X25PG, X24PG-X28PG); mean group disparity .036 over
all pairs and .039 over qualified pairs; smallest unadjusted heterogeneity p value .11.
