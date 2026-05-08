library(data.table)
library(lubridate)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

root_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
input_dir <- file.path(root_dir, "input_data")
output_dir <- file.path(root_dir, "output_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

local_dir <- file.path(root_dir, "external_data", "ERA5_local_hourly_byLocalYear")
coords_csv <- file.path(input_dir, "coords.csv")

res_deg <- 0.25
t_low <- 10
t_high <- 30
edh_thr <- c(30, 35)

out_bin <- file.path(output_dir, "env_EDH_by_relGDD_100bins.csv")
out_miss <- file.path(output_dir, "missing_local_files.csv")

snap <- function(x, res = 0.25) res * floor(x / res + 0.5)

make_local_path <- function(gp_lat, gp_lon, year, base_dir = local_dir) {
  file.path(base_dir, sprintf("ERA5_local_%.2f_%.2f_%d.csv", gp_lat, gp_lon, year))
}

read_local_year_hourly <- function(path) {
  dt <- fread(path, select = c("time_local", "t2m_c"))
  dt[, time_local := suppressWarnings(parse_date_time(
    time_local,
    orders = c("Y/m/d H:M", "Y/m/d H:M:S", "Y-m-d H:M", "Y-m-d H:M:S"),
    tz = "UTC"
  ))]

  if (all(is.na(dt$time_local))) {
    stop("Failed to parse time_local in file: ", path)
  }

  dt <- dt[!is.na(time_local)]
  dt[, t2m_c := as.numeric(t2m_c)]
  dt[, date_local := as.Date(time_local)]
  setorder(dt, time_local)
  dt[, .(time_local, date_local, t2m_c)]
}

gdd_inc_hour <- function(temp_c, t_low = t_low, t_high = t_high) {
  temp_star <- pmin(t_high, pmax(t_low, temp_c))
  pmax(0, (temp_star - t_low) / 24)
}

edh_hour <- function(temp_c, threshold) pmax(0, temp_c - threshold)

summarise_env_to_bins <- function(hourly_dt, sow, harv, env_id, thresholds = c(30, 35)) {
  sow_t <- as.POSIXct(sow, tz = "UTC")
  harv_t <- as.POSIXct(harv + days(1), tz = "UTC")

  d <- hourly_dt[time_local >= sow_t & time_local < harv_t]
  if (nrow(d) == 0) {
    return(data.table(
      env_ID = env_id,
      bin = 1:100,
      bin_hours = 0L,
      EDH30_sum = 0,
      EDH35_sum = 0,
      EDH30_per_hour = NA_real_,
      EDH35_per_hour = NA_real_,
      totalGDD = NA_real_,
      total_hours = 0L
    ))
  }

  d[, GDD_h := gdd_inc_hour(t2m_c)]
  total_gdd <- sum(d$GDD_h, na.rm = TRUE)

  if (!is.finite(total_gdd) || total_gdd <= 0) {
    return(data.table(
      env_ID = env_id,
      bin = 1:100,
      bin_hours = 0L,
      EDH30_sum = 0,
      EDH35_sum = 0,
      EDH30_per_hour = NA_real_,
      EDH35_per_hour = NA_real_,
      totalGDD = total_gdd,
      total_hours = nrow(d)
    ))
  }

  d[, cumGDD := cumsum(GDD_h)]
  d[, rel_progress := cumGDD / total_gdd]
  d[, bin := pmin(100L, pmax(1L, as.integer(floor(rel_progress * 100) + 1L)))]

  d[, EDH30_h := edh_hour(t2m_c, thresholds[1])]
  d[, EDH35_h := edh_hour(t2m_c, thresholds[2])]

  agg <- d[, .(
    bin_hours = .N,
    EDH30_sum = sum(EDH30_h, na.rm = TRUE),
    EDH35_sum = sum(EDH35_h, na.rm = TRUE)
  ), by = .(bin)]

  full <- data.table(bin = 1:100)[agg, on = "bin"]
  full[is.na(bin_hours), `:=`(bin_hours = 0L, EDH30_sum = 0, EDH35_sum = 0)]
  full[, EDH30_per_hour := fifelse(bin_hours > 0, EDH30_sum / bin_hours, NA_real_)]
  full[, EDH35_per_hour := fifelse(bin_hours > 0, EDH35_sum / bin_hours, NA_real_)]

  full[, `:=`(env_ID = env_id, totalGDD = total_gdd, total_hours = nrow(d))]
  setcolorder(full, c(
    "env_ID", "bin", "bin_hours",
    "EDH30_sum", "EDH35_sum", "EDH30_per_hour", "EDH35_per_hour",
    "totalGDD", "total_hours"
  ))
  full[]
}

coords <- fread(coords_csv)
coords[, sowing_date := ymd(sowing_date)]
coords[, harvest_date := ymd(harvest_date)]
coords <- coords[!is.na(env_ID) & !is.na(sowing_date) & !is.na(harvest_date)]
coords <- coords[harvest_date >= sowing_date]

setorder(coords, env_ID, sowing_date, harvest_date)
coords <- coords[, .SD[1], by = env_ID]
coords[, gp_lat := snap(as.numeric(lat), res_deg)]
coords[, gp_lon := snap(as.numeric(lon), res_deg)]
coords[, gp_id := sprintf("ERA5_%.2f_%.2f", gp_lat, gp_lon)]

cat("Unique env_ID:", nrow(coords), "| Unique grids:", uniqueN(coords$gp_id), "\n")

cache_env <- new.env(parent = emptyenv())
miss_log <- list()
res_list <- vector("list", nrow(coords))
m <- 0L

get_grid_year_dt <- function(gp_lat, gp_lon, gp_id, year_value) {
  key <- sprintf("%s_%d", gp_id, year_value)
  if (exists(key, envir = cache_env, inherits = FALSE)) {
    return(get(key, envir = cache_env, inherits = FALSE))
  }

  path <- make_local_path(gp_lat, gp_lon, year_value, local_dir)
  if (!file.exists(path)) return(NULL)
  dt <- try(read_local_year_hourly(path), silent = TRUE)
  if (inherits(dt, "try-error")) return(NULL)
  assign(key, dt, envir = cache_env)
  dt
}

for (i in seq_len(nrow(coords))) {
  env_id <- coords$env_ID[i]
  sow <- coords$sowing_date[i]
  harv <- coords$harvest_date[i]
  gp_lat <- coords$gp_lat[i]
  gp_lon <- coords$gp_lon[i]
  gp_id <- coords$gp_id[i]

  y1 <- year(sow)
  y2 <- year(harv)
  years_to_read <- if (y1 == y2) c(y1) else c(y1, y2)

  parts <- list()
  miss_years <- integer(0)
  for (yy in years_to_read) {
    dt_year <- get_grid_year_dt(gp_lat, gp_lon, gp_id, yy)
    if (is.null(dt_year)) {
      miss_years <- c(miss_years, yy)
    } else {
      parts[[length(parts) + 1L]] <- dt_year
    }
  }

  if (length(miss_years) > 0) {
    m <- m + 1L
    miss_log[[m]] <- data.table(
      env_ID = env_id,
      gp_id = gp_id,
      gp_lat = gp_lat,
      gp_lon = gp_lon,
      missing_years = paste(miss_years, collapse = ",")
    )
  }

  if (length(parts) == 0) {
    out <- data.table(
      env_ID = env_id,
      bin = 1:100,
      bin_hours = 0L,
      EDH30_sum = 0,
      EDH35_sum = 0,
      EDH30_per_hour = NA_real_,
      EDH35_per_hour = NA_real_,
      totalGDD = NA_real_,
      total_hours = 0L
    )
  } else {
    hourly_all <- rbindlist(parts, use.names = TRUE, fill = TRUE)
    setorder(hourly_all, time_local)
    hourly_all <- unique(hourly_all, by = "time_local")
    out <- summarise_env_to_bins(hourly_all, sow, harv, env_id, thresholds = edh_thr)
  }

  out[, `:=`(
    gp_id = gp_id,
    gp_lat = gp_lat,
    gp_lon = gp_lon,
    sowing_date = sow,
    harvest_date = harv,
    cross_year = (y1 != y2)
  )]
  setcolorder(out, c(
    "env_ID", "gp_id", "gp_lat", "gp_lon", "cross_year",
    "sowing_date", "harvest_date",
    "bin", "bin_hours",
    "EDH30_sum", "EDH35_sum", "EDH30_per_hour", "EDH35_per_hour",
    "totalGDD", "total_hours"
  ))

  res_list[[i]] <- out
  if (i %% 200 == 0) cat("Processed env:", i, "/", nrow(coords), "\n")
}

res_all <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
fwrite(res_all, out_bin)

if (length(miss_log) > 0) {
  fwrite(rbindlist(miss_log), out_miss)
} else {
  fwrite(data.table(), out_miss)
}

cat("Saved:", out_bin, "\n")
cat("Missing log:", out_miss, "\n")
cat("Rows:", nrow(res_all), "\n")
