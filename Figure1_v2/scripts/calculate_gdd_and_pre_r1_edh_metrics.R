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
r1_col <- "R1_date"

t_low <- 10
t_high <- 30
edh_thr <- c(30, 35)
res_deg <- 0.25

out_res <- file.path(output_dir, "env_GDD_EDH_results_from_local.csv")
out_miss <- file.path(output_dir, "missing_local_files.csv")

norm_lon <- function(lon) {
  lon2 <- lon
  idx <- is.finite(lon2) & (lon2 < -180 | lon2 >= 180)
  if (any(idx)) lon2[idx] <- ((lon2[idx] + 180) %% 360) - 180
  lon2
}

snap <- function(x, res = 0.25) res * floor(x / res + 0.5)
clamp_t <- function(x, low, high) pmin(high, pmax(low, x))
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_sum_exceed <- function(x, thr) {
  if (all(is.na(x))) return(NA_real_)
  sum(pmax(0, x - thr), na.rm = TRUE)
}

make_local_path <- function(gp_lat, gp_lon, year, base_dir = local_dir) {
  file.path(base_dir, sprintf("ERA5_local_%.2f_%.2f_%d.csv", gp_lat, gp_lon, year))
}

read_local_year <- function(path) {
  dt <- fread(path)
  if (!all(c("date_local", "t2m_c") %in% names(dt))) {
    stop("Missing required columns in file: ", path)
  }
  dt[, date_local := as.Date(date_local)]
  dt[, t2m_c := as.numeric(t2m_c)]
  if (!("timezone" %in% names(dt))) dt[, timezone := NA_character_]
  dt[, .(date_local, t2m_c, timezone)]
}

hourly_to_daily_indices <- function(hourly_dt, t_low = 10, t_high = 30, edh_thr = c(30, 35)) {
  dt <- copy(hourly_dt)
  dt[, date_local := as.Date(date_local)]
  dt[, t2m_c := as.numeric(t2m_c)]

  daily <- dt[, {
    tmax <- safe_max(t2m_c)
    tmin <- safe_min(t2m_c)
    .(
      Tmax = tmax,
      Tmin = tmin,
      EDH30_d = safe_sum_exceed(t2m_c, edh_thr[1]),
      EDH35_d = safe_sum_exceed(t2m_c, edh_thr[2])
    )
  }, by = .(date_local)]

  daily[, Tmax_star := clamp_t(Tmax, t_low, t_high)]
  daily[, Tmin_star := clamp_t(Tmin, t_low, t_high)]
  daily[, GDD_d := pmax(0, (Tmax_star + Tmin_star) / 2 - t_low)]
  daily[, c("Tmax_star", "Tmin_star") := NULL]
  setorder(daily, date_local)
  daily[]
}

coords <- fread(coords_csv)
need_cols <- c("env_ID", "lat", "lon", "sowing_date", "harvest_date")
missing_cols <- setdiff(need_cols, names(coords))
if (length(missing_cols) > 0) {
  stop("coords.csv is missing required columns: ", paste(missing_cols, collapse = ", "))
}

coords[, lon := norm_lon(as.numeric(lon))]
coords[, lat := as.numeric(lat)]
coords[, sowing_date := ymd(sowing_date)]
coords[, harvest_date := ymd(harvest_date)]

if (r1_col %in% names(coords)) {
  coords[, R1_date := ymd(get(r1_col))]
} else {
  coords[, R1_date := as.Date(NA)]
}

coords <- coords[
  !is.na(env_ID) & !is.na(lat) & !is.na(lon) &
    !is.na(sowing_date) & !is.na(harvest_date)
]
coords <- coords[harvest_date >= sowing_date]
coords[, gp_lat := snap(lat, res_deg)]
coords[, gp_lon := snap(lon, res_deg)]
coords[, gp_id := sprintf("ERA5_%.2f_%.2f", gp_lat, gp_lon)]

cat("Environments:", nrow(coords), "| Unique grid points:", uniqueN(coords$gp_id), "\n")

gp_list <- split(coords, by = "gp_id", keep.by = TRUE)
res_list <- vector("list", length(gp_list))
miss_list <- list()
k <- 0L
m <- 0L

for (gp_name in names(gp_list)) {
  k <- k + 1L
  dgp <- gp_list[[gp_name]]

  gp_lat <- dgp$gp_lat[1]
  gp_lon <- dgp$gp_lon[1]
  years_needed <- sort(unique(unlist(mapply(
    function(a, b) seq(year(a), year(b)),
    dgp$sowing_date,
    dgp$harvest_date
  ))))

  parts <- list()
  missing_years <- integer(0)

  for (yy in years_needed) {
    path <- make_local_path(gp_lat, gp_lon, yy, local_dir)
    if (!file.exists(path)) {
      missing_years <- c(missing_years, yy)
      next
    }
    dt_year <- try(read_local_year(path), silent = TRUE)
    if (inherits(dt_year, "try-error")) {
      missing_years <- c(missing_years, yy)
      next
    }
    parts[[as.character(yy)]] <- dt_year
  }

  if (length(missing_years) > 0) {
    m <- m + 1L
    miss_list[[m]] <- data.table(
      gp_id = gp_name,
      gp_lat = gp_lat,
      gp_lon = gp_lon,
      missing_years = paste(sort(unique(missing_years)), collapse = ",")
    )
  }

  if (length(parts) == 0) {
    res_list[[k]] <- dgp[, .(
      env_ID, gp_id, gp_lat, gp_lon, site, region, sowing_date, harvest_date, R1_date,
      total_days = NA_integer_,
      total_GDD = NA_real_,
      preR1_days = NA_integer_,
      preR1_GDD = NA_real_,
      postR1_days = NA_integer_,
      postR1_GDD = NA_real_,
      total_EDH30 = NA_real_,
      total_EDH35 = NA_real_,
      preR1_EDH30 = NA_real_,
      preR1_EDH35 = NA_real_,
      postR1_EDH30 = NA_real_,
      postR1_EDH35 = NA_real_
    )]
    next
  }

  hourly_all <- rbindlist(parts, use.names = TRUE, fill = TRUE)
  daily_all <- hourly_to_daily_indices(hourly_all, t_low = t_low, t_high = t_high, edh_thr = edh_thr)
  setkey(daily_all, date_local)

  one_gp_res <- dgp[, {
    s <- sowing_date
    h <- harvest_date
    r <- R1_date

    sub_total <- daily_all[date_local >= s & date_local <= h]
    total_days <- nrow(sub_total)

    total_GDD <- if (total_days == 0) NA_real_ else sum(sub_total$GDD_d, na.rm = TRUE)
    total_EDH30 <- if (total_days == 0) NA_real_ else sum(sub_total$EDH30_d, na.rm = TRUE)
    total_EDH35 <- if (total_days == 0) NA_real_ else sum(sub_total$EDH35_d, na.rm = TRUE)

    if (is.na(r) || r < s || r > h) {
      preR1_days <- NA_integer_
      preR1_GDD <- NA_real_
      preR1_EDH30 <- NA_real_
      preR1_EDH35 <- NA_real_
      postR1_days <- NA_integer_
      postR1_GDD <- NA_real_
      postR1_EDH30 <- NA_real_
      postR1_EDH35 <- NA_real_
    } else {
      sub_pre <- daily_all[date_local >= s & date_local <= r]
      preR1_days <- nrow(sub_pre)
      preR1_GDD <- if (preR1_days == 0) 0 else sum(sub_pre$GDD_d, na.rm = TRUE)
      preR1_EDH30 <- if (preR1_days == 0) 0 else sum(sub_pre$EDH30_d, na.rm = TRUE)
      preR1_EDH35 <- if (preR1_days == 0) 0 else sum(sub_pre$EDH35_d, na.rm = TRUE)

      sub_post <- daily_all[date_local > r & date_local <= h]
      postR1_days <- nrow(sub_post)
      postR1_GDD <- if (postR1_days == 0) 0 else sum(sub_post$GDD_d, na.rm = TRUE)
      postR1_EDH30 <- if (postR1_days == 0) 0 else sum(sub_post$EDH30_d, na.rm = TRUE)
      postR1_EDH35 <- if (postR1_days == 0) 0 else sum(sub_post$EDH35_d, na.rm = TRUE)
    }

    .(
      total_days = total_days,
      total_GDD = total_GDD,
      preR1_days = preR1_days,
      preR1_GDD = preR1_GDD,
      postR1_days = postR1_days,
      postR1_GDD = postR1_GDD,
      total_EDH30 = total_EDH30,
      total_EDH35 = total_EDH35,
      preR1_EDH30 = preR1_EDH30,
      preR1_EDH35 = preR1_EDH35,
      postR1_EDH30 = postR1_EDH30,
      postR1_EDH35 = postR1_EDH35
    )
  }, by = .(env_ID, gp_id, gp_lat, gp_lon, site, region, sowing_date, harvest_date, R1_date)]

  res_list[[k]] <- one_gp_res
  if (k %% 20 == 0) cat("Processed grid:", k, "/", length(gp_list), "\n")
}

res_all <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
fwrite(res_all, out_res)

if (length(miss_list) > 0) {
  miss_all <- rbindlist(miss_list, use.names = TRUE, fill = TRUE)
  fwrite(miss_all, out_miss)
} else {
  fwrite(data.table(), out_miss)
}

cat("Saved results:", out_res, "\n")
cat("Saved missing-file log:", out_miss, "\n")
cat("Rows:", nrow(res_all), "\n")
