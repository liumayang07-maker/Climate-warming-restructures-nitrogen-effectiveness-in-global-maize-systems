# Figure1

This directory is the public-release package for the current `Figure 1` workflow.

## Scope

Figure 1 summarizes the environmental heat-exposure classification, the overall nitrogen-response distribution, and the thermal-progress EDH pattern used to define pre-flowering and grain-filling risk contexts.

The public package includes the scripted ERA5/GDD/EDH workflow, the minimum `lnRR` histogram input, and derived CSV outputs needed to audit the Figure 1 numerical panels.

## Script Inventory

- `scripts/download_era5_by_gp_year.py`
- `scripts/convert_utc_to_local_and_standardize_units.R`
- `scripts/calculate_gdd_and_pre_r1_edh_metrics.R`
- `scripts/calculate_relative_thermal_progress_edh.R`
- `scripts/plot_edh_heatmap_and_trends.R`
- `scripts/plot_lnrr_histograms.R`
- `scripts/plot_pre_r1_gdd_ratio_boxplot.R`

## Packaged Input Data

- `input_data/coords.csv`
- `input_data/figure1_lnrr_histogram_input.csv`

The full `meta_data_v2.csv` is not included. `figure1_lnrr_histogram_input.csv` keeps only `lnRR` and `Env_Type`, which are the columns required by `plot_lnrr_histograms.R`.

## Packaged Output Data

- `output_data/env_EDH_by_relGDD_100bins.csv`
- `output_data/env_GDD_EDH_results_from_local.csv`

These files are derived outputs from the thermal-time and EDH workflow and can be compared with the relevant blocks in the manuscript Source Data workbook.

## Source Data Mapping

The manuscript Source Data workbook contains a `Figure 1` sheet with final plotted or panel-ready data. The relevant blocks are:

| Source Data block label | Public file or source |
| --- | --- |
| `Figure 1a_Data set` | site coordinates and counts derived from `coords.csv`; the final map was prepared in ArcGIS |
| `Figure 1b` | nitrogen-response distribution summary derived from the Figure 1 histogram workflow |
| `Figure 1c_Data set` | `output_data/env_EDH_by_relGDD_100bins.csv` |
| `Figure 1d_Data set` | binned EDH curve summaries from `plot_edh_heatmap_and_trends.R` |
| `Figure 1f` | environment-group summary statistics |

## ERA5 Handling

Raw and converted hourly ERA5 files are not stored in this GitHub package. The expected external directory layout is documented in `external_data/README.md`.

Reviewers can regenerate the climate inputs from `coords.csv` using:

1. `scripts/download_era5_by_gp_year.py`
2. `scripts/convert_utc_to_local_and_standardize_units.R`
3. the GDD/EDH calculation scripts listed above

Copernicus CDS credentials must be configured outside this repository.

## Exclusions

- complete `meta_data_v2.csv`
- raw and local-hour ERA5 CSV files
- ArcGIS project files
- rendered figure files and editable figure layouts
- local credentials such as `.cdsapirc`
