# DEP pipeline map

| Order | Script | Scientific role | Main output class |
|---:|---|---|---|
| 1 | `01_DEP_primary_analysis.R` | Input harmonization, primary limma DEP, primary sensitivities and enrichment | Private analysis results/workspace |
| 2 | `02_DEP_rebuild_sensitivity_fixed_map.R` | Rebuilds gene-level sensitivity outputs using the fixed primary gene–SOMAmer map | Private sensitivity results |
| 3 | `03_DEP_APOE_ATN_equal_sample_models.R` | Separates complete-case restriction from APOE/AT(N) covariate adjustment | Private sensitivity results |
| 4 | `04_DEP_country_site_robustness.R` | LOCO, LOSO, country meta-analysis, interaction and balanced resampling | Private robustness results |
| 5 | `05_DEP_matching_biomarker_sensitivity.R` | Propensity-score-matched and biomarker-compatible internal sensitivities | Private sensitivity results |
| 6 | `06_DEP_ptau217_threshold_sweep.R` | Exploratory cohort-defined p-tau217 enrichment sweep | Private sensitivity results |
| 7 | `07_DEP_demographic_balance.R` | Country/site demographic balance audit | Private robustness results |
| 8 | `08_DEP_generate_supplementary_tables.R` | Consolidated DEP Supplementary Tables | Publication candidate |
| 9 | `09_DEP_generate_source_data.R` | Source Data workbooks for quantitative DEP panels | Publication candidate |
| 10 | `10_DEP_generate_supplementary_data.R` | Complete high-dimensional Supplementary Data workbooks | Publication candidate |
| 11 | `11_DEP_generate_main_figure2.R` | Main Figure 2 from frozen Source Data | Publication candidate |
| 12 | `12_DEP_generate_extended_data_figures.R` | DEP Extended Data figures from frozen Source Data | Publication candidate |
| 13 | `13_DEP_simulation_based_detectable_effect_analysis.R` | Optional retrospective simulation-based sensitivity and detectable-effect analysis reproducing the primary limma and BH workflow | Private sensitivity results |

The legacy fragmented Supplementary Tables generator is retained under `archive/DEP_legacy/` and is not executed.