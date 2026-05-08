# Nitrogen management meta-analysis public release

This repository is a private staging area for reviewer- and publication-ready public materials from the nitrogen management meta-analysis project. It is intentionally kept private while files are screened and can be switched to public after scripts, source data, figures, and documentation have been checked for reproducibility, licensing, and sensitive information.

## Current contents

This staging version contains the code scripts from the `_v2` figure-analysis directories. These scripts cover the analytical workflow described in the manuscript `Data and code availability` and `Statistical analysis` sections, including climate-data matching, thermal-time reconstruction, EDH/GDD calculation, effect-size analysis, mixed-effects modelling, subgroup analyses, statistical testing, and most scripted figure generation.

This staging version now includes the Figure 1 public input and derived output tables needed for the scripted Figure 1 workflows. The complete project-level `meta_data_v2.csv` is not included; the Figure 1 histogram workflow uses a minimal two-column table instead.

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

### Figure 2

- `Figure2_v2/scripts/lnrr_unweighted_mixed_effects.R`
- `Figure2_v2/scripts/aen_unweighted_mixed_effects.R`
- `Figure2_v2/scripts/year_ge5_supplementary_effects.R`

### Figure 3

- `Figure3_v2/scripts/calculate_area_weighted_regional_n_input.R`
- `Figure3_v2/scripts/plot_regional_edh_intensity_stack.R`
- `Figure3_v2/scripts/analyze_regional_lnrr_and_yield.R`
- `Figure3_v2/scripts/analyze_regional_aen_and_yield.R`

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

### Extended Data Figure 2

- `Extended Data Figure2_v2/scripts/analyze_n_strategy_lnrr_and_yield.R`
- `Extended Data Figure2_v2/scripts/analyze_n_strategy_aen.R`
- `Extended Data Figure2_v2/scripts/analyze_n_strategy_env_contrasts_within_strategy.R`

## Not yet included

- Source Data file for Figure 2-4, Extended Data Figure 1-2 and Supplementary Figure 2-5
- Supplementary Tables 1-4
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
