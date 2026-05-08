# Climate warming restructures nitrogen effectiveness in global maize systems

This repository is the screened public-release package for the manuscript **Climate warming restructures nitrogen effectiveness in global maize systems**. It is currently kept as a private staging repository while the files are checked for reproducibility, naming consistency, licensing, and sensitive information before public release.

The project combines paired maize nitrogen-response observations from the literature with hourly ERA5 climate data to evaluate how phenology-defined heat exposure changes nitrogen effectiveness. Site-years are classified into normal, pre-flowering risk (PFR), and grain-filling risk (GFR) regimes; the figure packages document the corresponding thermal-time reconstruction, effect-size analyses, mixed-effects models, regional subset analyses, and adaptive-management comparisons.

This repository contains figure-linked scripts, minimal public input tables, and Source Data-aligned CSV outputs. It does not contain the complete curated project data set, draft manuscripts, editable figure files, local credentials, or large regenerated climate files.

## Repository structure

| Directory | Contents |
| --- | --- |
| `Figure1/` | ERA5/GDD/EDH workflow, Figure 1 histogram input, and derived thermal-exposure outputs. |
| `Figure2/` | Main Figure 2 `lnRR` and `AEN` mixed-effects model inputs and numerical outputs. |
| `Figure3/` | Regionalized Figure 3 package, Supplementary Figures 3-5 outputs, and Supplementary Tables 1-3. |
| `Figure4/` | Adaptive-strategy yield comparison package for Figure 4. |
| `Extended Data Figure1/` | Plant-density package for Extended Data Figure 1 and Supplementary Table 4. |
| `Extended Data Figure2/` | Nitrogen-strategy package for Extended Data Figure 2. |

Each figure directory contains its own `README.md` with the relevant script inventory, public input files, Source Data mapping, and figure-specific exclusions.

## Workflow summary

The analytical workflow combined literature-derived agronomic observations, ERA5 climate extraction, thermal-time reconstruction, heat-exposure classification, effect-size calculation, mixed-effects modelling, subgroup analyses, statistical testing, and scripted figure-data generation.

Analyses were conducted in R. ERA5 climate data were downloaded with a Python script using `cdsapi` to access the Copernicus Climate Data Store API. Growing degree days (GDD) and extreme degree hours (EDH) were calculated from local hourly temperature after UTC-to-local conversion.

Mixed-effects models were fitted with `lme` from `nlme`. Estimated marginal means and contrasts were obtained with `emmeans`; Wilcoxon, Dunn, and Games-Howell tests were used for non-model summaries or Figure 4 strategy comparisons where appropriate. Figure 1A was prepared manually in ArcGIS, and the final rendering of Figure 4C,E was completed in Origin 2025; the underlying numerical data are provided as CSV files.

## Public-data policy

The full project-level `meta_data_v2.csv` is not included. Instead, this repository provides figure-specific minimum input tables containing only the columns needed to rerun or audit the public scripts.

The complete manuscript Source Data workbook is provided separately with the manuscript. This repository provides figure-specific CSV exports and scripts for reviewer inspection and reproducibility.

Large ERA5 hourly files are not stored in GitHub. They can be regenerated from the public coordinate table and the ERA5 download/conversion scripts after configuring Copernicus CDS credentials outside this repository. If no-download reproduction is required later, these climate files should be deposited in an external archive such as Zenodo, figshare, a journal data package, or Git LFS after screening.

## Exclusions

The public-release package intentionally excludes:

- complete curated project data tables such as `meta_data_v2.csv`
- draft manuscripts, cover letters, internal review notes, and manuscript Word/PDF files
- full Source Data workbooks and local helper workbooks
- editable or rendered figure assets such as `*.ai`, `*.jpg`, `*.pdf`, `*.svg`, and `*.opju`
- local credentials, `.cdsapirc`, API keys, tokens, and virtual environments
- copyright-restricted literature PDFs or third-party images

## Reproducibility notes

The CSV files in each `output_data/` directory are intended to be compared with the manuscript Source Data workbook. Files named with current figure, extended-data, supplementary-figure, or supplementary-table labels follow the latest manuscript numbering. Where older internal Source Data block labels were retained for auditability, the local figure README documents that mapping.

This repository is designed for reviewer inspection and publication transparency, while protecting the non-public curated data set and avoiding large generated files in Git history.
