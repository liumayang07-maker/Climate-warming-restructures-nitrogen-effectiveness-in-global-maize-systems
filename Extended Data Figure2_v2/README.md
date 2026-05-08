# Extended Data Figure2_v2

This directory is the public-release package for `Extended Data Figure 2` and its directly linked country-specific follow-up analyses.

## Scope

The current manuscript places these analyses under `Nitrogen strategy and integrated management strategy analyses`, but the extracted Methods file still contains only the section heading rather than finalized text. Based on the current manuscript structure, figure assets, and executable scripts, the `Extended Data Figure 2` workflow retained here is the nitrogen-management-strategy branch, not the broader integrated-management branch.

The public workflow in this package is therefore aligned to the following analysis logic:

- restriction to the regional conventional-N subset (`region_N_label == "RN"`)
- combined analysis using observations from China and the USA
- comparison among four nitrogen-management strategies in the combined and China-only subsets: `S`, `SV`, `SR`, and `SVR`
- USA-specific follow-up analysis with sparse enhanced strategies collapsed into `SS`
- unweighted mixed-effects models for `lnRR` and `AEN`
- BH-adjusted Wilcoxon tests for observed yield under N fertilization
- environment-specific estimates within `Normal`, `PFR`, and `GFR`

## Audit Summary

After review against the current manuscript and the original `Extended Data Figure 2/` directory:

- the two primary scripts retained here match the current `Extended Data Figure 2` figure chain
- the public scripts now use repository-relative paths instead of machine-specific absolute paths
- `EDH` terminology has been synchronized in the packaged scripts and input table
- panel-linked CSV outputs were renamed in the scripts to the current `Extended Data Figure 2a-i` naming scheme
- several files from the original directory were intentionally excluded because they correspond to historical exploratory branches rather than the current figure

## Excluded Legacy Items

The following asset families were not included in the public figure package because they do not match the current `Extended Data Figure 2` manuscript workflow:

- `SEM_*` tables and SEM path-diagram PDFs
- `lavaan_*` outputs
- `HY_strategy` / `HRN` / `LRN` / `CP` / `HP` yield-comparison outputs
- the legacy integrated-management trend script from the original directory

These items belong to older integrated-management or exploratory branches and should not be treated as the current source of truth for `Extended Data Figure 2`.

## Directory Layout

- `scripts/`: public-facing analysis scripts for the nitrogen-management-strategy workflow
- `input_data/`: packaged input table used by the public scripts
- `output_data/`: Source Data-aligned numerical CSV exports corresponding to the current figure
- `figures/`: excluded from this public staging package unless requested later
- `docs/`: reserved for future release notes

## Script Inventory

### `scripts/analyze_n_strategy_lnrr_and_yield.R`

This script covers:

- the combined China+USA `lnRR` analysis
- the China-only `lnRR` follow-up analysis
- the USA-only `lnRR` follow-up analysis with `SV`, `SR`, and `SVR` collapsed into `SS`
- observed-yield comparisons for the combined, China-only, and USA-only subsets

### `scripts/analyze_n_strategy_aen.R`

This script covers:

- the combined China+USA `AEN` analysis
- the China-only `AEN` follow-up analysis
- the USA-only `AEN` follow-up analysis with `SV`, `SR`, and `SVR` collapsed into `SS`

## Packaged Input Data

- `input_data/Extended_Data_Figure2_n_strategy_model_input.csv`

The full `meta_data_v2.csv` is not included. `Extended_Data_Figure2_n_strategy_model_input.csv` keeps only the China/USA observations within the representative regional N-input window (`region_N_label == RN`) and the columns required by the public N-strategy scripts. It uses the `EDH35_total` and `EDH35_label` terminology adopted in the current manuscript-facing public directories.

## Packaged Output Data

`output_data/` contains the canonical panel-linked CSV files exported from `稿件/Source data.xlsx` (sheet `Extended Data Figure 2`) for the current `Extended Data Figure 2` package:

- `Extended_Data_Figure2_data_set.csv` — model input dataset
- `Extended_Data_Figure2a_lnrr_china_usa_model_results.csv`
- `Extended_Data_Figure2b_aen_china_usa_model_results.csv`
- `Extended_Data_Figure2c_lnrr_china_model_results.csv`
- `Extended_Data_Figure2d_aen_china_model_results.csv`
- `Extended_Data_Figure2e_lnrr_usa_model_results.csv`
- `Extended_Data_Figure2f_aen_usa_model_results.csv`
- `Extended_Data_Figure2g_yield_china_usa_wilcoxon_tests.csv`
- `Extended_Data_Figure2h_yield_china_wilcoxon_tests.csv`
- `Extended_Data_Figure2i_yield_usa_wilcoxon_tests.csv`

These are the final Source Data panel exports. Script-run CSV/PDF outputs are retained locally but not duplicated in this public output_data directory.

### Extended Data Figure 2 Source Data mapping

| Public file | Source Data sheet | Source Data block label | Rows x cols | Source |
| --- | --- | --- | --- | --- |
| `Extended_Data_Figure2_data_set.csv` | `Extended Data Figure 2` | `Figure 5 Selected Data set` | 1152 x 20 | Source Data block |
| `Extended_Data_Figure2a_lnrr_china_usa_model_results.csv` | `Extended Data Figure 2` | `Figure 5a_China+USA` | 24 x 17 | Source Data block |
| `Extended_Data_Figure2b_aen_china_usa_model_results.csv` | `Extended Data Figure 2` | `Figure 5b_China+USA` | 24 x 14 | Source Data block |
| `Extended_Data_Figure2c_lnrr_china_model_results.csv` | `Extended Data Figure 2` | `Figure 5c_China` | 24 x 17 | Source Data block |
| `Extended_Data_Figure2d_aen_china_model_results.csv` | `Extended Data Figure 2` | `Figure 5d_China` | 24 x 14 | Source Data block |
| `Extended_Data_Figure2e_lnrr_usa_model_results.csv` | `Extended Data Figure 2` | `Figure 5e_USA` | 18 x 17 | Source Data block |
| `Extended_Data_Figure2f_aen_usa_model_results.csv` | `Extended Data Figure 2` | `Figure 5f_USA` | 18 x 14 | Source Data block |
| `Extended_Data_Figure2g_yield_china_usa_wilcoxon_tests.csv` | `Extended Data Figure 2` | `Figure 5g_China+USA` | 32 x 9 | Source Data block |
| `Extended_Data_Figure2h_yield_china_wilcoxon_tests.csv` | `Extended Data Figure 2` | `Figure5h_China` | 32 x 11 | Source Data block |
| `Extended_Data_Figure2i_yield_usa_wilcoxon_tests.csv` | `Extended Data Figure 2` | `Figure5i_USA` | 11 x 11 | Source Data block |

Some original Source Data block labels use legacy `Figure 5` naming from an earlier manuscript draft. Those labels are preserved only in the mapping table for auditability and are not used as public file names.

Important Extended Data Figure 2 boundaries:

- Figure image files (`*.ai`, `*.jpg`, `*.pdf`) are not included at this stage.
- Large intermediate script-run outputs and legacy exploratory outputs are not staged here.
- The latest `Supplemental Figure.docx` contains Supplementary Figure 1-5 only, so there is no separate Supplementary Figure 6/7 package for Extended Data Figure 2.
- Extended Data Figure 2 uses model M6 from `Supplementary Table 1`; `../Figure3_v2/supplementary_tables/Supplementary_Table1_model_specifications.csv` is already included in this staging repository.
- This package provides the current Source Data-aligned CSV files and the scripts needed to inspect or regenerate the scripted components.

## Reproducibility Notes

- The scripts in this directory now use repository-relative paths and read from `input_data/Extended_Data_Figure2_n_strategy_model_input.csv` instead of the full `meta_data_v2.csv`.
- The public scripts are intended to reproduce the current nitrogen-management-strategy figure branch.
- The original project directory still contains legacy exploratory files that are deliberately excluded here.
- This packaging step did not rerun the analyses; it reorganized and cleaned the current manuscript-matched workflow for public release.

