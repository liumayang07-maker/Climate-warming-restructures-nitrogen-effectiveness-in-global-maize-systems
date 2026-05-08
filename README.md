# Nitrogen management meta-analysis public release

This repository is a private staging area for reviewer- and publication-ready public materials from the nitrogen management meta-analysis project. It is intentionally kept private while files are screened and can be switched to public after scripts, source data, figures, and documentation have been checked for reproducibility, licensing, and sensitive information.

## Current contents

This staging version contains the code scripts from the `_v2` figure-analysis directories. These scripts cover the analytical workflow described in the manuscript `Data and code availability` and `Statistical analysis` sections, including climate-data matching, thermal-time reconstruction, EDH/GDD calculation, effect-size analysis, mixed-effects modelling, subgroup analyses, statistical testing, and most scripted figure generation.

This staging version now includes the Figure 1 public input and derived output tables, the Figure 2 main-panel model input and numerical outputs, the Figure 3 regionalized source-data extracts and supplementary-table CSV files, Extended Data Figure 1 with Supplementary Table 4, and the Extended Data Figure 2 nitrogen-strategy package. The complete project-level `meta_data_v2.csv` is not included; figure workflows use minimum public input tables instead.

## Statistical workflow covered by the scripts

All analyses were conducted in R 4.4.2. ERA5 climate data were downloaded with a Python script using `cdsapi` to access the Copernicus Climate Data Store API. Growing degree days (GDD) and extreme degree hours (EDH) were calculated from local hourly temperature after UTC-to-local conversion.

Linear mixed-effects models were fitted with `lme` in the `nlme` package. Alternative fixed-effect structures were compared using AIC and ANOVA under maximum likelihood, and selected models were refitted using restricted maximum likelihood for estimation. Estimated marginal means and 95% confidence intervals were obtained with `emmeans`. Pairwise model contrasts were Tukey adjusted, and omnibus tests within environment groups were based on joint tests of `emmeans` contrasts. Raw-yield and other non-model summaries were assessed with Wilcoxon tests with Benjamini-Hochberg adjustment. Figure 4 strategy comparisons used `rstatix` for Games-Howell tests of mean yield differences and Dunn tests with Benjamini-Hochberg correction for rank-based distributional differences.

The site-distribution map in Figure 1A was prepared manually in ArcGIS, and the final rendering of Figure 4C,E was completed in Origin 2025. The underlying plotting data for these figures are expected to be provided separately in the manuscript Source Data file.

## Script inventory

### Figure 1

- `Figure1_v2/scripts/download_era5_by_gp_year.py`
- `Figure1_v2/scripts/convert_utc_to_local_and_standardize_units.R`
- `Figure1_v2/scripts/calculate_gdd_and_pre_r1_edh_metrics.R`
- `Figure1_v2/scripts/calculate_relative_thermal_progress_edh.R`
- `Figure1_v2/scripts/plot_edh_heatmap_and_trends.R`
- `Figure1_v2/scripts/plot_lnrr_histograms.R`
- `Figure1_v2/scripts/plot_pre_r1_gdd_ratio_boxplot.R`

Public Figure 1 inputs and derived tables currently included:

- `Figure1_v2/input_data/coords.csv`
- `Figure1_v2/input_data/figure1_lnrr_histogram_input.csv`
- `Figure1_v2/output_data/env_EDH_by_relGDD_100bins.csv`
- `Figure1_v2/output_data/env_GDD_EDH_results_from_local.csv`

The full `meta_data_v2.csv` is not included in this public-release staging repository. `figure1_lnrr_histogram_input.csv` keeps only `lnRR` and `Env_Type`, which are the columns required by `plot_lnrr_histograms.R`.

ERA5 hourly files are not included in GitHub. The expected local paths are documented in `Figure1_v2/external_data/README.md`; they can be regenerated from `coords.csv` with `download_era5_by_gp_year.py` and `convert_utc_to_local_and_standardize_units.R` after configuring Copernicus CDS credentials outside this repository.

### Figure 1 reproducibility and Source Data notes

The manuscript Source Data workbook contains a `Figure 1` sheet with the final plotted or panel-ready data for Figure 1. The sheet is organized as separate side-by-side blocks:

- `Figure 1a_Data set`: site coordinates and point counts for the manually prepared site-distribution map.
- `Figure 1b`: summary statistics for the overall nitrogen response distribution.
- `Figure 1c_Data set`: EDH-by-relative-thermal-progress data corresponding to `Figure1_v2/output_data/env_EDH_by_relGDD_100bins.csv`.
- `Figure 1d_Data set`: cluster-level binned EDH curve summaries derived from `plot_edh_heatmap_and_trends.R`.
- `Figure 1f`: environment-group summary statistics.

The public repository supports Figure 1 review in two complementary ways:

- Numerical verification: reviewers can compare the `Figure 1` sheet in the manuscript Source Data file against the public CSV files. In particular, `Figure 1c_Data set` matches the columns and rows of the public `env_EDH_by_relGDD_100bins.csv`, except that Excel may display date columns as date values while the CSV stores them as ISO-style text.
- Scripted reproduction: the public scripts can regenerate the scripted Figure 1 outputs from the public inputs that are included here. The histogram script uses `figure1_lnrr_histogram_input.csv`, a minimal row-level file containing only `lnRR` and `Env_Type`, rather than the full project-level `meta_data_v2.csv`.

Important boundaries:

- The complete `meta_data_v2.csv` is not included because it is the curated project data set. Only the minimum row-level information needed for the Figure 1 histogram is included.
- Raw hourly ERA5 files and converted local-hour ERA5 files are not bundled in GitHub. They can be regenerated from `coords.csv` using the CDS download and UTC-to-local conversion scripts, or provided later through an external archive if required.
- Figure 1A map layout was prepared manually in ArcGIS. The public materials provide the plotted coordinate/count data, but not an ArcGIS project file or a fully scripted map workflow.
- The public repository is intended to support reviewer inspection and reproducibility of the scripted Figure 1 panels, while preserving the non-public full curated data set and avoiding large regenerated ERA5 files in Git history.

### Figure 2

- `Figure2_v2/scripts/lnrr_unweighted_mixed_effects.R`
- `Figure2_v2/scripts/aen_unweighted_mixed_effects.R`
- `Figure2_v2/scripts/year_ge5_supplementary_effects.R`

Public Figure 2 main-panel inputs and numerical outputs currently included:

- `Figure2_v2/input_data/figure2_model_input.csv`
- `Figure2_v2/output_data/lnrr/*.csv`
- `Figure2_v2/output_data/lnrr/*.txt`
- `Figure2_v2/output_data/aen/*.csv`
- `Figure2_v2/output_data/aen/*.txt`

The full `meta_data_v2.csv` is not included in this public-release staging repository. `figure2_model_input.csv` keeps only the row-level columns needed by the main Figure 2 lnRR and AEN mixed-effects scripts: `lnRR`, `AEN`, `Env_Type`, `reference_ID`, `site_year_ID`, `env_ID`, `EDH35_total`, `N_rate`, and `N_rate_cluster`.

The manuscript Source Data workbook contains a `Figure 2` sheet with the final panel-ready data for Figure 2. The public CSV/TXT outputs in `output_data/lnrr/` and `output_data/aen/` are the script-generated numerical outputs used to trace the main Figure 2 panels, including the model data set, environment-level estimates, nitrogen-rate-cluster estimates, AIC model comparisons, pairwise contrasts, and QC summaries.

Important Figure 2 boundaries:

- `year_ge5_supplementary_effects.R` is retained as a script, but its input and `output_data/year_ge5/` results are not included in this Figure 2 main-panel package. They will be handled with the relevant supplementary figure materials.
- Figure files (`*.pdf`, `*.ai`, `*.jpg`) are not included at this stage. If journal or reviewer requirements make them necessary, they can be screened and added separately or archived externally.

### Figure 3

- `Figure3_v2/scripts/calculate_area_weighted_regional_n_input.R`
- `Figure3_v2/scripts/plot_regional_edh_intensity_stack.R`
- `Figure3_v2/scripts/analyze_regional_lnrr_and_yield.R`
- `Figure3_v2/scripts/analyze_regional_aen_and_yield.R`

Public Figure 3 inputs currently included:

- `Figure3_v2/input_data/Figure3_edh_stack_input.csv`
- `Figure3_v2/input_data/Figure3_regional_model_input.csv`
- `Figure3_v2/input_data/Figure3_regional_n_input_1990_2020.csv`
- `Figure3_v2/input_data/Figure3_regional_n_input_2011_2020.csv`

The full `meta_data_v2.csv` and `region_meta_data_v2.csv` are not included. The public input files keep only the minimum columns needed for the public Figure 3 scripts and for reviewer inspection of the regionalized workflow.

Public Figure 3 and regional supplementary outputs currently included:

- `Figure3_v2/output_data/Figure3*.csv`
- `Figure3_v2/output_data/Supplementary_Figure3*.csv`
- `Figure3_v2/output_data/Supplementary_Figure4*.csv`
- `Figure3_v2/output_data/Supplementary_Figure5*.csv`
- `Figure3_v2/supplementary_tables/Supplementary_Table1_model_specifications.csv`
- `Figure3_v2/supplementary_tables/Supplementary_Table2_regional_observation_balance.csv`
- `Figure3_v2/supplementary_tables/Supplementary_Table3_regional_conventional_n_input_windows.csv`

The file names follow the current manuscript figure and table titles in `Figure.docx`, `Supplemental Figure.docx`, and `Supplemental Table.docx`. Historical internal output names from earlier scripts are not used as public file names. Some original Source Data block labels contain spelling errors such as `Fgiure`; those labels are retained only in the mapping below for auditability.

#### Figure 3 Source Data mapping

| Public file | Source Data sheet | Source Data block label | Rows x cols | Source |
| --- | --- | --- | --- | --- |
| `Figure3a_data_set.csv` | `Figure3` | `Figure3a_data set` | 110 x 14 | Source Data block |
| `Figure3a_region_environment_counts.csv` | `Figure3` | `Figure3a` | 11 x 10 | Source Data block |
| `Figure3b_edh_intensity_stack.csv` | `Figure3` | `Figure3b` | 77 x 4 | Source Data block / `plot_regional_edh_intensity_stack.R` |
| `Figure3c_data_set.csv` | `Figure3` | `Figure3c` occurrence 1 | 4170 x 25 | Source Data block |
| `Figure3c_model_results.csv` | `Figure3` | `Figure3c` occurrence 2 | 39 x 14 | Source Data block |
| `Figure3d_yield_wilcoxon_tests.csv` | `Figure3` | `Figure3d` | 31 x 10 | Source Data block |
| `Figure3e_data_set.csv` | `Figure3` | `Figure3e` occurrence 1 | 1601 x 25 | Source Data block |
| `Figure3e_model_results.csv` | `Figure3` | `Figure3e` occurrence 2 | 39 x 14 | Source Data block |
| `Figure3f_yield_wilcoxon_tests.csv` | `Figure3` | `Figure3f` | 31 x 10 | Source Data block |
| `Supplementary_Figure3a_n_rate_all_observations.csv` | `Supplement Figure 3` | `Supplement Figure 3a` | 15 x 10 | Source Data block |
| `Supplementary_Figure3b_n_rate_regionalized_dataset.csv` | `Supplement Figure 3` | `Supplement Fgiure 3b` | 31 x 10 | Source Data block |
| `Supplementary_Figure3c_zero_n_yield_all_observations.csv` | `Supplement Figure 3` | `Supplement Fgiure 3c` | 15 x 10 | Source Data block |
| `Supplementary_Figure3d_zero_n_yield_regionalized_dataset.csv` | `Supplement Figure 3` | `Supplement Fgiure 3d` | 31 x 10 | Source Data block |
| `Supplementary_Figure4a_aen_all_observations_model_results.csv` | `Supplement Figure 4` | `Supplement Fgiure 4a` | 39 x 11 | Source Data block |
| `Supplementary_Figure4b_aen_regionalized_dataset_model_results.csv` | `Supplement Figure 4` | `Supplement Fgiure 4b` | 39 x 11 | Source Data block |
| `Supplementary_Figure5_data_set.csv` | `Supplement Figure 5` | `Supplement Fgiure 5_data_set` | 1164 x 21 | Source Data block |
| `Supplementary_Figure5a_lnrr_china_usa_model_results.csv` | `Supplement Figure 5` | `Supplement Fgiure 5a` | 21 x 14 | Source Data block |
| `Supplementary_Figure5b_aen_china_usa_model_results.csv` | `Supplement Figure 5` | `Supplement Fgiure 5b` | 21 x 14 | Source Data block |
| `Supplementary_Figure5c_fertilized_yield_wilcoxon_tests.csv` | `Supplement Figure 5` | `Supplement Fgiure 5c` | 13 x 10 | Source Data block |
| `Supplementary_Figure5d_zero_n_yield_wilcoxon_tests.csv` | `Supplement Figure 5` | `Supplement Fgiure 5d` | 13 x 10 | Source Data block |

Important Figure 3 boundaries:

- `Supplementary Table 4` is not included in this Figure 3 package because it belongs to the plant-density / Extended Data Figure 1 workflow.
- `Supplementary Figure 2` is not included here because it is part of the thermal-progress heat-exposure workflow rather than the regionalized Figure 3 workflow.
- Figure image files (`*.pdf`, `*.ai`, `*.jpg`) are not included at this stage. The public package provides the current Source Data-aligned CSV files and the scripts needed to inspect or regenerate the scripted components.

### Figure 4

- `Figure4_v2/scripts/compare_adaptation_strategies_yield.R`

### Extended Data Figure 1

- `Extended Data Figure1_v2/scripts/analyze_density_lnrr_and_yield.R`
- `Extended Data Figure1_v2/scripts/analyze_density_lnrr_and_yield_china.R`
- `Extended Data Figure1_v2/scripts/analyze_density_lnrr_and_yield_usa.R`
- `Extended Data Figure1_v2/scripts/analyze_density_lnrr_and_yield_china_usa.R`
- `Extended Data Figure1_v2/scripts/analyze_density_aen.R`
- `Extended Data Figure1_v2/scripts/analyze_density_aen_china.R`
- `Extended Data Figure1_v2/scripts/analyze_density_aen_usa.R`
- `Extended Data Figure1_v2/scripts/analyze_density_aen_china_usa.R`
- `Extended Data Figure1_v2/scripts/analyze_density_env_contrasts_china.R`
- `Extended Data Figure1_v2/scripts/analyze_density_env_contrasts_usa.R`
- `Extended Data Figure1_v2/scripts/analyze_density_env_contrasts_china_usa.R`
- `Extended Data Figure1_v2/scripts/leave_one_reference_sensitivity_china_usa_lnRR.R`

Public Extended Data Figure 1 inputs currently included:

- `Extended Data Figure1_v2/input_data/Extended_Data_Figure1_density_model_input.csv`

The full `meta_data_v2.csv` is not included. `Extended_Data_Figure1_density_model_input.csv` keeps only the China/USA observations within the representative regional N-input window (`region_N_label == RN`) and the columns required by the plant-density mixed-effects, yield-summary and sensitivity scripts. The public scripts read this minimum input table and use `plant_density_cluster` as the primary planting-density class, while retaining compatibility with the older `plant_density_3cluster` column name if needed.

Public Extended Data Figure 1 Source Data extracts currently included:

- `Extended Data Figure1_v2/output_data/Extended_Data_Figure1*.csv`
- `Extended Data Figure1_v2/supplementary_tables/Supplementary_Table4_plant_density_windows.csv`

The `Extended Data Figure 1` sheet in the manuscript Source Data workbook contains final panel-ready data for Extended Data Fig. 1. Some original block labels in that sheet still use legacy `Figure4*` names; those labels are preserved only in the mapping table below for exact auditability and are not used as public file names.

#### Extended Data Figure 1 Source Data mapping

| Public file | Source Data sheet | Source Data block label | Rows x cols | Source |
| --- | --- | --- | --- | --- |
| `Extended_Data_Figure1_data_set.csv` | `Extended Data Figure 1` | `Figure 4 Selected Data set` | 1066 x 19 | Source Data block |
| `Extended_Data_Figure1a_lnrr_china_usa_model_results.csv` | `Extended Data Figure 1` | `Figure 4a_China+USA` | 27 x 16 | Source Data block |
| `Extended_Data_Figure1b_aen_china_usa_model_results.csv` | `Extended Data Figure 1` | `Figure4b_China+USA` | 27 x 16 | Source Data block |
| `Extended_Data_Figure1c_lnrr_china_model_results.csv` | `Extended Data Figure 1` | `Figure4c_China` | 27 x 16 | Source Data block |
| `Extended_Data_Figure1d_aen_china_model_results.csv` | `Extended Data Figure 1` | `Figure4d_China` | 27 x 15 | Source Data block |
| `Extended_Data_Figure1e_lnrr_usa_model_results.csv` | `Extended Data Figure 1` | `Figure4e_USA` | 27 x 16 | Source Data block |
| `Extended_Data_Figure1f_aen_usa_model_results.csv` | `Extended Data Figure 1` | `Figure4f_USA` | 27 x 16 | Source Data block |
| `Extended_Data_Figure1g_yield_china_usa_wilcoxon_tests.csv` | `Extended Data Figure 1` | `Figure4g_China+USA` | 20 x 10 | Source Data block |
| `Extended_Data_Figure1h_yield_china_wilcoxon_tests.csv` | `Extended Data Figure 1` | `Figure4h_China` | 20 x 11 | Source Data block |
| `Extended_Data_Figure1i_yield_usa_wilcoxon_tests.csv` | `Extended Data Figure 1` | `Figure4i_USA` | 20 x 9 | Source Data block |
| `Supplementary_Table4_plant_density_windows.csv` | `Supplemental Table.docx` | `Supplementary Table 4` | 15 x 5 | Latest supplementary table |

Important Extended Data Figure 1 boundaries:

- Figure image files (`*.pdf`, `*.ai`, `*.jpg`) are not included at this stage.
- Older local files named `Supplementary_Figure6_China.*` and `Supplementary_Figure7_USA.*` are not included because the latest `Supplemental Figure.docx` currently contains Supplementary Figure 1-5 only.
- Large sensitivity output tables are not staged here unless separately required for reviewer QA.

### Extended Data Figure 2

- `Extended Data Figure2_v2/scripts/analyze_n_strategy_lnrr_and_yield.R`
- `Extended Data Figure2_v2/scripts/analyze_n_strategy_aen.R`
- `Extended Data Figure2_v2/scripts/analyze_n_strategy_env_contrasts_within_strategy.R`

Public Extended Data Figure 2 inputs currently included:

- `Extended Data Figure2_v2/input_data/Extended_Data_Figure2_n_strategy_model_input.csv`

The full `meta_data_v2.csv` is not included. `Extended_Data_Figure2_n_strategy_model_input.csv` keeps only China/USA observations within the representative regional N-input window (`region_N_label == RN`) and the nitrogen-strategy classes used by the Extended Data Figure 2 scripts (`S`, `SV`, `SR`, `SVR`). The public scripts read this minimum input table directly.

Public Extended Data Figure 2 Source Data extracts currently included:

- `Extended Data Figure2_v2/output_data/Extended_Data_Figure2*.csv`

The `Extended Data Figure 2` sheet in the manuscript Source Data workbook contains final panel-ready data for Extended Data Fig. 2. Some original block labels in that sheet still use legacy `Figure 5*` names; those labels are preserved only in the mapping table below for exact auditability and are not used as public file names.

#### Extended Data Figure 2 Source Data mapping

| Public file | Source Data sheet | Source Data block label | Rows x cols | Source |
| --- | --- | --- | --- | --- |
| `Extended_Data_Figure2_data_set.csv` | `Extended Data Figure 2` | `Figure 5 Selected Data set` | 1152 x 21 | Source Data block |
| `Extended_Data_Figure2a_lnrr_china_usa_model_results.csv` | `Extended Data Figure 2` | `Figure 5a_China+USA` | 24 x 18 | Source Data block |
| `Extended_Data_Figure2b_aen_china_usa_model_results.csv` | `Extended Data Figure 2` | `Figure 5b_China+USA` | 24 x 15 | Source Data block |
| `Extended_Data_Figure2c_lnrr_china_model_results.csv` | `Extended Data Figure 2` | `Figure 5c_China` | 24 x 18 | Source Data block |
| `Extended_Data_Figure2d_aen_china_model_results.csv` | `Extended Data Figure 2` | `Figure 5d_China` | 24 x 15 | Source Data block |
| `Extended_Data_Figure2e_lnrr_usa_model_results.csv` | `Extended Data Figure 2` | `Figure 5e_USA` | 18 x 18 | Source Data block |
| `Extended_Data_Figure2f_aen_usa_model_results.csv` | `Extended Data Figure 2` | `Figure 5f_USA` | 18 x 15 | Source Data block |
| `Extended_Data_Figure2g_yield_china_usa_wilcoxon_tests.csv` | `Extended Data Figure 2` | `Figure 5g_China+USA` | 32 x 11 | Source Data block |
| `Extended_Data_Figure2h_yield_china_wilcoxon_tests.csv` | `Extended Data Figure 2` | `Figure5h_China` | 32 x 12 | Source Data block |
| `Extended_Data_Figure2i_yield_usa_wilcoxon_tests.csv` | `Extended Data Figure 2` | `Figure5i_USA` | 11 x 11 | Source Data block |

Important Extended Data Figure 2 boundaries:

- Figure image files (`*.pdf`, `*.ai`, `*.jpg`) are not included at this stage.
- The latest `Supplemental Figure.docx` contains Supplementary Figure 1-5 only, so there is no separate Supplementary Figure 6/7 package for Extended Data Figure 2.
- Extended Data Figure 2 uses model M6 from `Supplementary Table 1`; `Figure3_v2/supplementary_tables/Supplementary_Table1_model_specifications.csv` is already included in this staging repository.

## Not yet included

- Full Source Data workbook for Figure 2-4, Extended Data Figure 1-2 and Supplementary Figure 2-5
- Later figure-specific supplementary-table packages not yet staged
- Final figure files and editable figure layouts
- Draft manuscript, cover letter and internal review files
- External raw/local-hour ERA5 files and local credentials for Copernicus CDS access

## Excluded from public release

- Draft manuscripts and cover letters
- Internal review notes
- Raw credentials or local environment files
- API keys, tokens, and `.cdsapirc`
- Local virtual environments
- Copyright-restricted literature PDFs or third-party images
