library(data.table)
library(lubridate)
library(lutz)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

root_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
external_dir <- file.path(root_dir, "external_data")
era_dir <- file.path(external_dir, "ERA5_timeseries_byGPYear_raw_utc")
out_dir <- file.path(external_dir, "ERA5_local_hourly_byLocalYear")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(out_dir, "build_local_dataset_log.csv")

norm_lon <- function(lon) ((lon + 180) %% 360) - 180

detect_time_col <- function(names_vec) {
  candidates <- c("valid_time", "time", "datetime", "date_time", "time_utc")
  direct_hit <- names_vec[tolower(names_vec) %in% candidates]
  if (length(direct_hit) > 0) return(direct_hit[1])
  fuzzy_hit <- names_vec[grepl("time", tolower(names_vec))]
  if (length(fuzzy_hit) > 0) return(fuzzy_hit[1])
  NA_character_
}

make_era_path <- function(lat, lon, year, base_dir = era_dir) {
  file.path(base_dir, sprintf("ERA5_%.2f_%.2f_%d.csv", lat, lon, year))
}

parse_utc_time <- function(x) {
  if (inherits(x, "POSIXct")) return(with_tz(x, "UTC"))
  x <- as.character(x)
  x <- gsub("T", " ", x, fixed = TRUE)
  x <- gsub("Z$", "", x)

  parsed <- suppressWarnings(parse_date_time(
    x,
    orders = c(
      "Ymd HMS", "Ymd HM",
      "Y/m/d H:M:S", "Y/m/d H:M",
      "Y-m-d H:M:S", "Y-m-d H:M"
    ),
    tz = "UTC",
    exact = FALSE
  ))
  as.POSIXct(parsed, tz = "UTC")
}

read_gp_year <- function(path) {
  dt <- fread(path)
  time_col <- detect_time_col(names(dt))
  if (is.na(time_col)) stop("No time column detected in file: ", path)
  if (!("t2m" %in% names(dt))) stop("Missing t2m in file: ", path)
  if (!all(c("latitude", "longitude") %in% names(dt))) {
    stop("Missing latitude/longitude columns in file: ", path)
  }

  dt[, longitude := norm_lon(longitude)]
  dt[, time_utc := parse_utc_time(get(time_col))]
  dt <- dt[!is.na(time_utc)]
  if (nrow(dt) == 0) stop("All time_utc values were parsed as NA in file: ", path)

  if (!("tp" %in% names(dt))) dt[, tp := NA_real_]
  if (!("ssrd" %in% names(dt))) dt[, ssrd := NA_real_]
  dt[, .(time_utc, latitude, longitude, t2m, tp, ssrd)]
}

build_one_gp_local_year <- function(gp_lat, gp_lon, year_value, era_dir, out_dir, tz_cache) {
  gid <- sprintf("ERA5_%.2f_%.2f", gp_lat, gp_lon)
  out_file <- file.path(out_dir, sprintf("ERA5_local_%.2f_%.2f_%d.csv", gp_lat, gp_lon, year_value))

  if (file.exists(out_file) && file.info(out_file)$size > 0) {
    return(list(status = "skip", out = out_file, match_rate = NA_real_, msg = "already exists"))
  }

  if (gid %in% names(tz_cache)) {
    tz_name <- tz_cache[[gid]]
  } else {
    tz_name <- lutz::tz_lookup_coords(lat = gp_lat, lon = gp_lon, method = "fast")
    if (is.na(tz_name) || tz_name == "") tz_name <- "UTC"
    tz_cache[[gid]] <- tz_name
  }

  local_start <- ymd_hms(sprintf("%d-01-01 00:00:00", year_value), tz = tz_name)
  local_end <- ymd_hms(sprintf("%d-12-31 23:00:00", year_value), tz = tz_name)

  full_local <- data.table(time_local = seq(local_start, local_end, by = "1 hour"))
  full_local[, time_utc_key := with_tz(time_local, "UTC")]
  full_local[, time_utc_src := floor_date(time_utc_key, "hour")]
  full_local[, date_local := as.Date(time_local)]

  utc_years <- sort(unique(year(full_local$time_utc_src)))
  parts <- list()
  missing_files <- character(0)

  for (yy in utc_years) {
    path <- make_era_path(gp_lat, gp_lon, yy, era_dir)
    if (!file.exists(path)) {
      missing_files <- c(missing_files, basename(path))
      next
    }
    dt_year <- try(read_gp_year(path), silent = TRUE)
    if (inherits(dt_year, "try-error")) {
      missing_files <- c(missing_files, paste0(basename(path), " [read_error]"))
      next
    }
    parts[[length(parts) + 1L]] <- dt_year
  }

  dt_all <- if (length(parts) > 0) rbindlist(parts, use.names = TRUE, fill = TRUE) else NULL
  if (is.null(dt_all) || nrow(dt_all) == 0) {
    out_na <- full_local[, .(
      time_utc = format(time_utc_key, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      time_local = format(time_local, "%Y-%m-%d %H:%M:%S", tz = tz_name),
      date_local = format(date_local, "%Y-%m-%d"),
      latitude = gp_lat,
      longitude = gp_lon,
      t2m_c = NA_real_,
      tp_mm = NA_real_,
      ssrd_mj = NA_real_,
      timezone = tz_name
    )]
    fwrite(out_na, out_file)
    msg <- if (length(missing_files) > 0) {
      paste("missing:", paste(missing_files, collapse = ";"))
    } else {
      "no source files"
    }
    return(list(status = "ok", out = out_file, match_rate = 0, msg = msg))
  }

  dt_all[, time_utc := floor_date(time_utc, "hour")]
  start_src <- min(full_local$time_utc_src)
  end_src <- max(full_local$time_utc_src)
  dt_all <- dt_all[time_utc >= start_src & time_utc <= end_src]
  setkey(dt_all, time_utc)
  dt_all <- unique(dt_all, by = "time_utc")

  setkey(full_local, time_utc_src)
  out <- dt_all[full_local, on = .(time_utc = time_utc_src)]
  out[is.na(latitude), latitude := gp_lat]
  out[is.na(longitude), longitude := gp_lon]
  out[, t2m_c := t2m - 273.15]
  out[, tp_mm := tp * 1000]
  out[, ssrd_mj := ssrd / 1e6]
  out[, timezone := tz_name]

  match_rate <- mean(!is.na(out$t2m))
  keep <- out[, .(
    time_utc = format(time_utc_key, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    time_local = format(time_local, "%Y-%m-%d %H:%M:%S", tz = tz_name),
    date_local = format(date_local, "%Y-%m-%d"),
    latitude,
    longitude,
    t2m_c,
    tp_mm,
    ssrd_mj,
    timezone
  )]
  fwrite(keep, out_file)

  msg_parts <- character(0)
  if (length(missing_files) > 0) {
    msg_parts <- c(msg_parts, paste("missing:", paste(missing_files, collapse = ";")))
  }
  msg_parts <- c(msg_parts, sprintf("match_rate=%.4f", match_rate))
  list(status = "ok", out = out_file, match_rate = match_rate, msg = paste(msg_parts, collapse = " | "))
}

files <- list.files(
  era_dir,
  pattern = "^ERA5_-?[0-9]+\\.[0-9]+_-?[0-9]+\\.[0-9]+_[0-9]{4}\\.csv$",
  full.names = FALSE
)

parse_one <- function(filename) {
  match <- regexec("^ERA5_(-?[0-9]+\\.[0-9]+)_(-?[0-9]+\\.[0-9]+)_([0-9]{4})\\.csv$", filename)
  groups <- regmatches(filename, match)[[1]]
  if (length(groups) != 4) return(NULL)
  list(
    gp_lat = as.numeric(groups[2]),
    gp_lon = as.numeric(groups[3]),
    year = as.integer(groups[4])
  )
}

task_list <- lapply(files, parse_one)
task_list <- task_list[!sapply(task_list, is.null)]
tasks <- rbindlist(task_list)
tasks <- unique(tasks, by = c("gp_lat", "gp_lon", "year"))
setorder(tasks, gp_lat, gp_lon, year)

cat("Tasks:", nrow(tasks), "\n")

tz_cache <- list()
log_dt <- data.table(
  gp_lat = numeric(),
  gp_lon = numeric(),
  year = integer(),
  status = character(),
  match_rate = numeric(),
  out = character(),
  msg = character()
)

for (i in seq_len(nrow(tasks))) {
  gp_lat <- tasks$gp_lat[i]
  gp_lon <- tasks$gp_lon[i]
  yy <- tasks$year[i]

  res <- try(build_one_gp_local_year(gp_lat, gp_lon, yy, era_dir, out_dir, tz_cache), silent = TRUE)
  if (inherits(res, "try-error")) {
    log_dt <- rbind(log_dt, data.table(
      gp_lat = gp_lat,
      gp_lon = gp_lon,
      year = yy,
      status = "fail",
      match_rate = NA_real_,
      out = NA_character_,
      msg = as.character(res)
    ))
  } else {
    log_dt <- rbind(log_dt, data.table(
      gp_lat = gp_lat,
      gp_lon = gp_lon,
      year = yy,
      status = res$status,
      match_rate = res$match_rate,
      out = res$out,
      msg = res$msg
    ))
  }

  if (i %% 100 == 0 || i == nrow(tasks)) {
    fwrite(log_dt, log_file)
    cat(sprintf("Processed %d / %d\n", i, nrow(tasks)))
  }
}

fwrite(log_dt, log_file)
cat("Saved log:", log_file, "\n")
