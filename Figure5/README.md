# Figure5

This directory is the public-release package for `Figure 5` and its directly linked country-specific follow-up analyses.

## Scope

The current public package is aligned to the nitrogen-management-strategy branch of the manuscript. It covers:

- restriction to the regional conventional-N subset (`region_N_label == "RN"`)
- combined analysis using observations from China and the USA
- comparison among four nitrogen-management strategies in the combined and China-only subsets: `S`, `SV`, `SR`, and `SVR`
- USA-specific follow-up analysis with sparse enhanced strategies collapsed into `SS`
- unweighted mixed-effects models for `lnRR` and `AEN`
- BH-adjusted Wilcoxon tests for observed yield under N fertilization
- environment-specific estimates within `Normal`, `PFR`, and `GFR`

## Script Inventory

- `scripts/analyze_n_strategy_lnrr_and_yield.R`
- `scripts/analyze_n_strategy_aen.R`
- `scripts/analyze_n_strategy_env_contrasts_within_strategy.R`

## Packaged Input Data

- `input_data/Figure5_n_strategy_model_input.csv`

The full `meta_data_v2.csv` is not included. The public input keeps only China/USA observations within the representative regional N-input window (`region_N_label == RN`) and the columns required by the public N-strategy scripts. It uses the `EDH35_total` and `EDH35_label` terminology adopted in the current manuscript-facing public directories.

## Packaged Output Data

`output_data/` contains the canonical panel-linked CSV files exported from the manuscript Source Data workbook, sheet `Figure 5`:

- `Figure5_data_set.csv`
- `Figure5a_lnrr_china_usa_model_results.csv`
- `Figure5b_aen_china_usa_model_results.csv`
- `Figure5c_lnrr_china_model_results.csv`
- `Figure5d_aen_china_model_results.csv`
- `Figure5e_lnrr_usa_model_results.csv`
- `Figure5f_aen_usa_model_results.csv`
- `Figure5g_yield_china_usa_wilcoxon_tests.csv`
- `Figure5h_yield_china_wilcoxon_tests.csv`
- `Figure5i_yield_usa_wilcoxon_tests.csv`
- `Figure5_source_data_mapping.csv`

## Source Data Mapping

The full mapping is provided in:

- `output_data/Figure5_source_data_mapping.csv`

The public file names and Source Data block labels follow the current manuscript numbering for Figure 5.

## Exclusions

- complete `meta_data_v2.csv`
- rendered figure files and editable figure layouts
- legacy SEM, `lavaan`, and integrated-management exploratory outputs
- large intermediate script-run outputs not required for the current Source Data package

The latest `Supplemental Figure.docx` contains Supplementary Figures 1-5 only, so there is no separate Supplementary Figure 6/7 package for Figure 5. Figure 5 uses model M6 from `Supplementary Table 1`; `../Figure3/supplementary_tables/Supplementary_Table1_model_specifications.csv` is already included in this staging repository.

## Reproducibility Notes

The scripts in this directory use repository-relative paths and read from `input_data/Figure5_n_strategy_model_input.csv` instead of the full `meta_data_v2.csv`. This package provides the current Source Data-aligned CSV files and scripts needed to inspect or regenerate the scripted components.
