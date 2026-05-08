# Figure1 external ERA5 data

This directory documents the expected location of ERA5 files used by the Figure 1 climate-processing scripts. The hourly ERA5 CSV files are not included in the GitHub repository because they can be regenerated from `Figure1/input_data/coords.csv` using the Copernicus Climate Data Store API and would make the repository large.

Expected subdirectories:

- `ERA5_timeseries_byGPYear_raw_utc/`
- Created by `Figure1/scripts/download_era5_by_gp_year.py`
  - Expected files: `ERA5_<lat>_<lon>_<year>.csv`
- `ERA5_local_hourly_byLocalYear/`
- Created by `Figure1/scripts/convert_utc_to_local_and_standardize_units.R`
  - Expected files: `ERA5_local_<lat>_<lon>_<year>.csv`

To regenerate these files, configure local Copernicus CDS credentials outside this repository and run the Figure 1 scripts in the order described in the root README.
