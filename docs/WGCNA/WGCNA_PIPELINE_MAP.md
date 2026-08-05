# WGCNA pipeline map

| Order | Script | Main purpose | Output class |
|---:|---|---|---|
| 1 | `01_WGCNA_prepare_input.R` | Build the outcome-independent, gene-collapsed input | Private intermediate |
| 2 | `02_WGCNA_construct_network.R` | Construct the signed WGCNA network | Private intermediate |
| 3 | `03_WGCNA_characterize_modules.R` | Compute kME, hubs, DEP burden and enrichment | Analysis output |
| 4 | `04_WGCNA_module_trait_models.R` | Test module–trait and context associations | Analysis output |
| 5 | `05_WGCNA_site_biomarker_sensitivity.R` | Run HC3, nested-site and biomarker sensitivity models | Analysis output |
| 6 | `06_WGCNA_association_stability.R` | Run LOCO, LOSO and balanced downsampling | Analysis output |
| 7 | `07_WGCNA_biomarker_fdr_correction.R` | Apply the 32-model biomarker FDR family | Analysis output |
| 8 | `08_WGCNA_structural_preservation.R` | Test fixed-gene country and site preservation | Analysis output |
| 9 | `09_WGCNA_network_quality.R` | Quantify modularity and topology | Analysis output |
| 10 | `10_WGCNA_generate_main_figure3.R` | Generate Figure 3 and source tables | Publication candidate |
| 11 | `11_WGCNA_generate_extended_data_figures.R` | Generate Extended Data figures and source tables | Publication candidate |
| 12 | `12_WGCNA_generate_supplementary_tables.R` | Generate Supplementary Tables 20–27 | Publication candidate |
| 13 | `13_WGCNA_generate_manuscript_text.R` | Generate methods, results and legends | Editorial support |
| 14 | `14_WGCNA_build_submission_package.R` | Assemble the WGCNA submission package | Publication candidate |
| 15 | `15_WGCNA_audit_submission_package.R` | Audit package completeness, science and privacy | Audit output |

The network is built before clinical association testing. Country and site analyses are internal stability analyses, not external validation.
