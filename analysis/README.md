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
