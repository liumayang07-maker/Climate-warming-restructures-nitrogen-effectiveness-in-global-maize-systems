# Extended Data Figure1

This directory is the public-release package for `Extended Data Figure 1` and the linked `Supplementary Table 4`.

## Scope

Extended Data Figure 1 evaluates plant-density responses within China and USA observations under representative regional nitrogen-input windows. The package contains minimum public input data, Source Data-aligned panel CSV files, and the plant-density supplementary table.

The public file names follow the current manuscript numbering. Some original Source Data block labels still use legacy `Figure4*` names from an earlier draft; those labels are retained only in the mapping CSV for auditability.

## Script Inventory

- `scripts/analyze_density_lnrr_and_yield.R`
- `scripts/analyze_density_lnrr_and_yield_china.R`
- `scripts/analyze_density_lnrr_and_yield_usa.R`
- `scripts/analyze_density_lnrr_and_yield_china_usa.R`
- `scripts/analyze_density_aen.R`
- `scripts/analyze_density_aen_china.R`
- `scripts/analyze_density_aen_usa.R`
- `scripts/analyze_density_aen_china_usa.R`
- `scripts/analyze_density_env_contrasts_china.R`
- `scripts/analyze_density_env_contrasts_usa.R`
- `scripts/analyze_density_env_contrasts_china_usa.R`
- `scripts/leave_one_reference_sensitivity_china_usa_lnRR.R`

## Packaged Input Data

- `input_data/Extended_Data_Figure1_density_model_input.csv`

The full `meta_data_v2.csv` is not included. The public input keeps only the China/USA observations within the representative regional N-input window (`region_N_label == RN`) and the columns required by the plant-density mixed-effects, yield-summary, environment-contrast, and sensitivity scripts.

The public scripts use `plant_density_cluster` as the primary density class and retain compatibility with the older `plant_density_3cluster` name where needed.

## Packaged Output Data

`output_data/` contains current canonical CSV exports for:

- `Extended_Data_Figure1_data_set.csv`
- `Extended_Data_Figure1a_lnrr_china_usa_model_results.csv`
- `Extended_Data_Figure1b_aen_china_usa_model_results.csv`
- `Extended_Data_Figure1c_lnrr_china_model_results.csv`
- `Extended_Data_Figure1d_aen_china_model_results.csv`
- `Extended_Data_Figure1e_lnrr_usa_model_results.csv`
- `Extended_Data_Figure1f_aen_usa_model_results.csv`
- `Extended_Data_Figure1g_yield_china_usa_wilcoxon_tests.csv`
- `Extended_Data_Figure1h_yield_china_wilcoxon_tests.csv`
- `Extended_Data_Figure1i_yield_usa_wilcoxon_tests.csv`
- `Extended_Data_Figure1_source_data_mapping.csv`

`supplementary_tables/` contains:

- `Supplementary_Table4_plant_density_windows.csv`

## Source Data Mapping

The manuscript Source Data workbook contains final panel-ready data in the `Extended Data Figure 1` sheet. The full mapping is provided in:

- `output_data/Extended_Data_Figure1_source_data_mapping.csv`

The linked supplementary table is exported from the latest `Supplemental Table.docx`:

| Public file | Source |
| --- | --- |
| `supplementary_tables/Supplementary_Table4_plant_density_windows.csv` | Supplementary Table 4 |

## Exclusions

- complete `meta_data_v2.csv`
- rendered figure files and editable figure layouts
- older local `Supplementary_Figure6_China.*` and `Supplementary_Figure7_USA.*` files
- large sensitivity-output tables not required for the current manuscript Source Data package

The latest `Supplemental Figure.docx` contains Supplementary Figures 1-5 only, so no Supplementary Figure 6/7 public package is included here.
