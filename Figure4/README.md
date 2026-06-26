# Figure4

This directory is the public-release package for `Figure 4` and the linked `Supplementary Table 4`.

## Scope

Figure 4 evaluates plant-density responses within China and USA observations under representative regional nitrogen-input windows. The package contains minimum public input data, Source Data-aligned panel CSV files, and the plant-density supplementary table.

The public file names and Source Data block labels follow the current manuscript numbering for Figure 4.

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

- `input_data/Figure4_density_model_input.csv`

The full `meta_data_v2.csv` is not included. The public input keeps only the China/USA observations within the representative regional N-input window (`region_N_label == RN`) and the columns required by the plant-density mixed-effects, yield-summary, environment-contrast, and sensitivity scripts.

The public scripts use `plant_density_cluster` as the primary density class and retain compatibility with the older `plant_density_3cluster` name where needed.

## Packaged Output Data

`output_data/` contains current canonical CSV exports for:

- `Figure4_data_set.csv`
- `Figure4a_lnrr_china_usa_model_results.csv`
- `Figure4b_aen_china_usa_model_results.csv`
- `Figure4c_lnrr_china_model_results.csv`
- `Figure4d_aen_china_model_results.csv`
- `Figure4e_lnrr_usa_model_results.csv`
- `Figure4f_aen_usa_model_results.csv`
- `Figure4g_yield_china_usa_wilcoxon_tests.csv`
- `Figure4h_yield_china_wilcoxon_tests.csv`
- `Figure4i_yield_usa_wilcoxon_tests.csv`
- `Figure4_source_data_mapping.csv`

`supplementary_tables/` contains:

- `Supplementary_Table4_plant_density_windows.csv`

## Source Data Mapping

The manuscript Source Data workbook contains final panel-ready data in the `Figure 4` sheet. The full mapping is provided in:

- `output_data/Figure4_source_data_mapping.csv`

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
