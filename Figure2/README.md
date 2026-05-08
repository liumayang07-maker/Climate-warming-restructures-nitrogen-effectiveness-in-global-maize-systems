# Figure2

This directory is the public-release package for the main `Figure 2` mixed-effects model workflow.

## Scope

Figure 2 evaluates nitrogen-response effect sizes and agronomic efficiency across environmental contexts and nitrogen-rate clusters. This package covers the main Figure 2 `lnRR` and `AEN` analyses.

The long-term supplementary script is retained for traceability, but its complete input table and output files are not included in this main Figure 2 package. It should be treated as a supplementary-workflow script unless the long-term experiment input package is released separately.

## Script Inventory

- `scripts/lnrr_unweighted_mixed_effects.R`
- `scripts/aen_unweighted_mixed_effects.R`
- `scripts/year_ge5_supplementary_effects.R`

## Packaged Input Data

- `input_data/figure2_model_input.csv`

The full `meta_data_v2.csv` is not included. `figure2_model_input.csv` keeps only the row-level columns needed by the main Figure 2 mixed-effects scripts:

- `lnRR`
- `AEN`
- `Env_Type`
- `reference_ID`
- `site_year_ID`
- `env_ID`
- `EDH35_total`
- `N_rate`
- `N_rate_cluster`

## Packaged Output Data

`output_data/lnrr/` contains the `lnRR` model data, AIC comparisons, fitted-model summaries, estimated marginal means, pairwise contrasts, and QC summaries.

`output_data/aen/` contains the corresponding `AEN` outputs.

These CSV/TXT files are script-generated numerical outputs used to trace the main Figure 2 panels against the `Figure 2` sheet in the manuscript Source Data workbook.

## Source Data Notes

The manuscript Source Data workbook contains a `Figure 2` sheet with final panel-ready data. The public output directories provide the underlying numerical model outputs for auditing:

- model data set and filtering checks
- overall and environment-specific estimates
- nitrogen-rate-cluster estimates
- AIC model comparisons
- pairwise environment contrasts
- row and trial count QC summaries

## Exclusions

- complete `meta_data_v2.csv`
- `meta_data_v2_filled_Data_year_duration.csv`, `output_data/year_ge5/`, and long-term supplementary outputs
- rendered figure files and editable figure layouts
- manuscript workbooks or draft manuscript files

The `year_ge5_supplementary_effects.R` script should be handled with the relevant supplementary figure package if those long-term materials are released later.
