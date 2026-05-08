# Figure3

This directory is the public-release package for the current `Figure 3` regionalized analysis and its directly linked supplementary materials.

## Scope

Figure 3 links regional heat exposure, regional nitrogen-input windows, `lnRR`, `AEN`, and observed-yield contrasts. This package also contains the current Source Data-aligned outputs for Supplementary Figures 3-5 and Supplementary Tables 1-3.

File names follow the latest manuscript figure and table numbering. Historical internal names or misspellings in Source Data block labels are preserved only in mapping files for auditability.

## Script Inventory

- `scripts/calculate_area_weighted_regional_n_input.R`
- `scripts/plot_regional_edh_intensity_stack.R`
- `scripts/analyze_regional_lnrr_and_yield.R`
- `scripts/analyze_regional_aen_and_yield.R`

## Packaged Input Data

- `input_data/Figure3_edh_stack_input.csv`
- `input_data/Figure3_regional_model_input.csv`
- `input_data/Figure3_regional_n_input_1990_2020.csv`
- `input_data/Figure3_regional_n_input_2011_2020.csv`

The full `meta_data_v2.csv` and `region_meta_data_v2.csv` are not included. The public inputs keep only the minimum columns needed for the regionalized Figure 3 scripts and reviewer inspection.

## Packaged Output Data

`output_data/` contains current canonical CSV exports for:

- `Figure3*.csv`
- `Supplementary_Figure3*.csv`
- `Supplementary_Figure4*.csv`
- `Supplementary_Figure5*.csv`
- `Figure3_source_data_mapping.csv`

`supplementary_tables/` contains:

- `Supplementary_Table1_model_specifications.csv`
- `Supplementary_Table2_regional_observation_balance.csv`
- `Supplementary_Table3_regional_conventional_n_input_windows.csv`

## Source Data Mapping

The manuscript Source Data workbook provides the panel-ready data in sheets `Figure3`, `Supplement Figure 3`, `Supplement Figure 4`, and `Supplement Figure 5`. The full mapping is provided in:

- `output_data/Figure3_source_data_mapping.csv`

Key public outputs include:

| Public file group | Source Data sheet |
| --- | --- |
| `Figure3a_*` to `Figure3f_*` | `Figure3` |
| `Supplementary_Figure3a_*` to `Supplementary_Figure3d_*` | `Supplement Figure 3` |
| `Supplementary_Figure4a_*` and `Supplementary_Figure4b_*` | `Supplement Figure 4` |
| `Supplementary_Figure5_*` | `Supplement Figure 5` |
| `Supplementary_Table1_*` to `Supplementary_Table3_*` | latest `Supplemental Table.docx` |

Some original Source Data block labels contain legacy spelling such as `Fgiure`. Those labels are retained in the mapping CSV only so that reviewers can match the public files to the exact workbook blocks.

## Exclusions

- complete `meta_data_v2.csv`
- complete `region_meta_data_v2.csv`
- raw spatial inputs used before regional aggregation
- rendered figure files and editable figure layouts
- `Supplementary Table 4`, which belongs to the Extended Data Figure 1 plant-density package

## Reproducibility Notes

The public package provides already aggregated regional nitrogen-input tables. The script `calculate_area_weighted_regional_n_input.R` is retained for transparency, but its original spatial input files are not included in this GitHub package.
