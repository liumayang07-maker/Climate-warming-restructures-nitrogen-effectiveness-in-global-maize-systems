import os
import random
import shutil
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed

import cdsapi
import pandas as pd


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, ".."))
INPUT_DIR = os.path.join(ROOT_DIR, "input_data")
EXTERNAL_DIR = os.path.join(ROOT_DIR, "external_data")
IN_CSV = os.path.join(INPUT_DIR, "coords.csv")
OUT_DIR = os.path.join(EXTERNAL_DIR, "ERA5_timeseries_byGPYear_raw_utc")
os.makedirs(OUT_DIR, exist_ok=True)

DATASET = "reanalysis-era5-single-levels-timeseries"
VARIABLES = [
    "2m_temperature",
    "total_precipitation",
    "surface_solar_radiation_downwards",
]

RES = 0.25
MAX_WORKERS = 2
RETRIES = 6
TEST_N_TASKS = None


def norm_lon(lon: float) -> float:
    return ((lon + 180.0) % 360.0) - 180.0


def snap(x: float, res: float = 0.25) -> float:
    return round(x / res) * res


def gp_id(lat_gp: float, lon_gp: float) -> str:
    return f"ERA5_{lat_gp:.2f}_{lon_gp:.2f}"


def csv_header_ok(csv_path: str) -> bool:
    if not (os.path.exists(csv_path) and os.path.getsize(csv_path) > 100):
        return False
    try:
        with open(csv_path, "r", encoding="utf-8", errors="ignore") as handle:
            header = handle.readline().strip()
        return ("valid_time" in header) or (header.lower().startswith("time")) or ("t2m" in header)
    except Exception:
        return False


def download_one_gp_year(lat_gp, lon_gp, year):
    gid = gp_id(lat_gp, lon_gp)
    csv_path = os.path.join(OUT_DIR, f"{gid}_{year}.csv")
    part_zip = os.path.join(OUT_DIR, f"{gid}_{year}.zip.part")
    tmp_dir = os.path.join(OUT_DIR, "_tmp_extract")
    os.makedirs(tmp_dir, exist_ok=True)

    if csv_header_ok(csv_path):
        return csv_path

    request = {
        "variable": VARIABLES,
        "location": {"latitude": float(lat_gp), "longitude": float(lon_gp)},
        "date": [f"{year}-01-01/{year}-12-31"],
        "data_format": "csv",
    }

    client = cdsapi.Client()

    for attempt in range(1, RETRIES + 1):
        try:
            client.retrieve(DATASET, request).download(part_zip)

            if not zipfile.is_zipfile(part_zip):
                raise RuntimeError("Downloaded file is not a valid zip archive.")

            with zipfile.ZipFile(part_zip) as zf:
                bad_entry = zf.testzip()
                if bad_entry is not None:
                    raise RuntimeError(f"Corrupted zip entry detected: {bad_entry}")

                csv_names = [name for name in zf.namelist() if name.lower().endswith(".csv")]
                if not csv_names:
                    raise RuntimeError("The downloaded zip archive does not contain a CSV file.")

                inner_csv = csv_names[0]
                zf.extract(inner_csv, tmp_dir)
                extracted = os.path.join(tmp_dir, inner_csv)

                if os.path.exists(csv_path):
                    os.remove(csv_path)
                shutil.move(extracted, csv_path)

            if os.path.exists(part_zip):
                os.remove(part_zip)

            if not csv_header_ok(csv_path):
                raise RuntimeError("The extracted CSV header is invalid.")

            return csv_path

        except Exception as exc:
            for file_path in (part_zip, csv_path):
                if os.path.exists(file_path):
                    try:
                        os.remove(file_path)
                    except Exception:
                        pass

            wait = min(180, 5 * (2 ** (attempt - 1))) + random.uniform(0, 3)
            print(f"[{gid} {year}] retry {attempt}/{RETRIES}: {exc} | sleep {wait:.1f}s")
            time.sleep(wait)

    return None


df = pd.read_csv(IN_CSV).drop_duplicates("env_ID").copy()
df["sowing_date"] = pd.to_datetime(df["sowing_date"], format="%Y/%m/%d", errors="coerce")
df["harvest_date"] = pd.to_datetime(df["harvest_date"], format="%Y/%m/%d", errors="coerce")
df = df.dropna(subset=["env_ID", "lat", "lon", "sowing_date", "harvest_date"])
df = df[df["harvest_date"] >= df["sowing_date"]].copy()

df["lon"] = df["lon"].astype(float).map(norm_lon)
df["lat"] = df["lat"].astype(float)
df["gp_lat"] = df["lat"].map(lambda x: snap(x, RES))
df["gp_lon"] = df["lon"].map(lambda x: snap(x, RES))
df["gp_id"] = df.apply(lambda row: gp_id(row["gp_lat"], row["gp_lon"]), axis=1)

rows = []
for _, row in df.iterrows():
    year_start = int(row["sowing_date"].year)
    year_end = int(row["harvest_date"].year)
    for year in range(year_start, year_end + 1):
        rows.append((row["gp_lat"], row["gp_lon"], row["gp_id"], year))

tasks = pd.DataFrame(rows, columns=["gp_lat", "gp_lon", "gp_id", "year"])
tasks = tasks.drop_duplicates(subset=["gp_id", "year"]).sort_values(["gp_id", "year"]).reset_index(drop=True)

print("Unique environments:", df["env_ID"].nunique())
print("Unique grid points:", df["gp_id"].nunique())
print("Tasks (grid point x year):", tasks.shape[0])

df_out = df[["env_ID", "lat", "lon", "gp_id", "gp_lat", "gp_lon", "sowing_date", "harvest_date"]].copy()
df_out["sowing_date"] = df_out["sowing_date"].dt.strftime("%Y-%m-%d")
df_out["harvest_date"] = df_out["harvest_date"].dt.strftime("%Y-%m-%d")
df_out.to_csv(os.path.join(OUT_DIR, "env_to_gp_mapping.csv"), index=False)
tasks.to_csv(os.path.join(OUT_DIR, "gp_year_tasks.csv"), index=False)

if TEST_N_TASKS is not None:
    tasks = tasks.head(TEST_N_TASKS).copy()
    print("TEST mode: tasks =", tasks.shape[0])

jobs = []
results = []
with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    for _, task in tasks.iterrows():
        jobs.append(executor.submit(download_one_gp_year, task["gp_lat"], task["gp_lon"], int(task["year"])))
    for future in as_completed(jobs):
        results.append(future.result())

ok = sum(1 for item in results if item)
print(f"Done: {ok}/{len(results)} | saved to: {OUT_DIR}")
