# Climate warming restructures nitrogen effectiveness in global maize systems

[中文说明](README.zh-CN.md)

This repository provides code and figure-level data exports for the manuscript **Climate warming restructures nitrogen effectiveness in global maize systems**.

The study combines paired maize nitrogen-response observations from the literature with hourly ERA5 climate data to evaluate how phenology-defined heat exposure changes nitrogen effectiveness. Site-years are classified into normal, pre-flowering risk (PFR), and grain-filling risk (GFR) regimes, and the analyses quantify how these heat-risk regimes alter yield response, agronomic nitrogen efficiency, regional response patterns, and candidate management options.

## Repository Contents

The repository is organized by manuscript figure. Each figure directory contains its own `README.md` describing the scripts, public input files, generated outputs, and Source Data mapping for that figure.

| Directory | Contents |
| --- | --- |
| `Figure1/` | ERA5/GDD/EDH processing workflow, Figure 1 histogram input, and derived thermal-exposure outputs. |
| `Figure2/` | Main Figure 2 `lnRR` and `AEN` mixed-effects model inputs and numerical outputs. |
| `Figure3/` | Regionalized Figure 3 analyses, Supplementary Figures 3-5 outputs, and Supplementary Tables 1-3. |
| `Figure4/` | Adaptive-management strategy yield-comparison package for Figure 4. |
| `Extended Data Figure1/` | Plant-density analyses for Extended Data Figure 1 and Supplementary Table 4. |
| `Extended Data Figure2/` | Nitrogen-strategy analyses for Extended Data Figure 2. |

## Data Availability

The complete manuscript Source Data workbook is provided separately with the manuscript. This repository provides figure-specific CSV exports and analysis scripts so that the numerical results can be inspected and reproduced at the figure level.

The full project-level curated data table (`meta_data_v2.csv`) is not included. Instead, figure-specific minimum input tables are provided with only the columns needed by the public scripts.

Raw hourly ERA5 files are not stored in this GitHub repository because they are large regenerated climate files. They can be rebuilt from the public coordinate table and the ERA5 download/conversion scripts after configuring Copernicus Climate Data Store credentials locally.

## Workflow Summary

The analytical workflow includes:

- literature-derived paired maize nitrogen-response observations
- ERA5 climate extraction and UTC-to-local hourly conversion
- growing degree day (GDD) and extreme degree hour (EDH) calculation
- phenology-defined heat-risk classification
- effect-size calculation using `lnRR` and agronomic nitrogen efficiency (`AEN`)
- mixed-effects modelling with `nlme`
- estimated marginal means and contrasts with `emmeans`
- Wilcoxon, Dunn, and Games-Howell tests for selected non-model summaries
- figure-level CSV export and scripted numerical checks

Figure 1A was prepared in ArcGIS, and the final rendering of Figure 4C,E was completed in Origin 2025. The underlying numerical data for these panels are provided as CSV files.

## Notes

This repository does not include editable figure artwork, manuscript documents, local credentials, virtual environments, copyright-restricted literature files, or large generated climate files.

File names follow the current manuscript figure, extended-data figure, supplementary-figure, and supplementary-table numbering. Where older Source Data block labels were retained for traceability, the relevant figure-level README documents the mapping.
