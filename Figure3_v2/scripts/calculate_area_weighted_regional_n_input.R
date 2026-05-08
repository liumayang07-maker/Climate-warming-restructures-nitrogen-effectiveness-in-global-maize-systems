library(sf)
library(dplyr)
library(tidyr)

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
external_dir <- file.path(root_dir, "external_data", "regional_n_input_raw")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
setwd(output_dir)

df <- readRDS(file.path(external_dir, "N_application_maize_masked_1990_2020.rds"))
admin1 <- st_read(file.path(external_dir, "ne_10m_admin_1_states_provinces", "ne_10m_admin_1_states_provinces.shp"))

sort(unique(admin1$name[admin1$admin == "China"]))

# Public-release note
regions_mapping <- list(
  Europe = list(admin = c("Portugal", "Spain", "Hungary", "Poland", "Italy", "Croatia", "Germany")),
  North_America_Irrigated = list(admin = "United States of America", "Canada", name = c("Kansas", "Nebraska", "Oklahoma", "Texas", "Colorado", "South Dakota", "Ontario", "Québec")),
  North_china_Plain = list(admin = "China", name = c("Beijing", "Tianjin", "Hebei", "Henan", "Shandong")),
  North_America_Rainfed = list(admin = "United States of America", name = c("Alabama", "Arkansas", "Illinois", "Indiana", "Iowa", "Kentucky",
                                                                  "Louisiana", "Minnesota", "Mississippi", "Missouri", "Ohio",
                                                                  "Virginia", "Wisconsin")),
  Northwest_China = list(admin = "China", name = c("Shaanxi","Shanxi", "Gansu", "Ningxia", "Qinghai", "Xinjiang")),
  Northeast_China = list(admin = "China", name = c("Heilongjiang", "Liaoning", "Jilin", "Inner Mongol")),
  Southwest_China = list(admin = "China", name = c("Sichuan", "Yunnan", "Guizhou", "Chongqing", "Guangxi")),
  South_America = list(admin = c("Argentina", "Brazil", "Colombia")),
  South_Asia = list(admin = c("India", "Iran", "Pakistan", "Nepal", "Bangladesh")),
  East_Africa = list(admin = c("Zambia", "Kenya", "Uganda", "Egypt", "Tanzania")),
  West_Africa = list(admin = c("Nigeria", "Niger", "Ghana"))
)

# Public-release note
get_region_boundary <- function(region, region_name) {
  filter_expr <- admin1$admin %in% region$admin
  if (!is.null(region$name)) filter_expr <- filter_expr & admin1$name %in% region$name
  d <- admin1[filter_expr, ]
  d$region <- region_name
  d
}
region_boundaries_named <- lapply(names(regions_mapping), function(name) {
  get_region_boundary(regions_mapping[[name]], name)
})
highlight_regions_colored <- do.call(rbind, region_boundaries_named)
highlight_regions_colored <- st_make_valid(highlight_regions_colored)
highlight_regions_colored_union <- highlight_regions_colored %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# Public-release note
pts_sf <- st_as_sf(df, coords = c("long", "lat"), crs = 4326, remove = FALSE)
pts_sf_region <- st_join(pts_sf, highlight_regions_colored_union["region"], left = TRUE)
df$region <- pts_sf_region$region
df <- df[!is.na(df$region), ]

# Public-release note
fert_vars <- grep("_(deep|surface)$", names(df), value = TRUE)

df <- df %>%
  filter(year >= 1990 & year <= 2020, Harvested_area > 0) %>%
  mutate(
    across(all_of(fert_vars), ~replace_na(., 0)),  # Public-release note
    N_synthetic = (
      AA_deep + AA_surface +
        AN_deep + AN_surface +
        AP_deep + AP_surface +
        AS_deep + AS_surface +
        CAN_deep + CAN_surface +
        NK_deep + NK_surface +
        NPK_deep + NPK_surface +
        NS_deep + NS_surface +
        ONP_deep + ONP_surface +
        ONS_deep + ONS_surface +
        Urea_deep + Urea_surface
    ) * Harvested_area,
    N_manure = (MA_deep + MA_surface) * Harvested_area,
    N_residue = (CR_deep + CR_surface) * Harvested_area,
    N_total = N_synthetic + N_manure + N_residue
  )

# Public-release note
region_year_nrate <- df %>%
  group_by(region, year) %>%
  summarise(
    Area_sum = sum(Harvested_area, na.rm = TRUE),
    N_total_sum = sum(N_total, na.rm = TRUE),
    N_synth_sum = sum(N_synthetic, na.rm = TRUE),
    N_manure_sum = sum(N_manure, na.rm = TRUE),
    N_residue_sum = sum(N_residue, na.rm = TRUE),
    N_total_rate = N_total_sum / Area_sum,
    N_synth_rate = N_synth_sum / Area_sum,
    N_manure_rate = N_manure_sum / Area_sum,
    N_residue_rate = N_residue_sum / Area_sum,
    .groups = "drop"
  )

# Public-release note
write.csv(region_year_nrate, "11region_year_weighted_N_rate_1990_2020.csv", row.names = FALSE)
cat("[DONE] /（kg N/ha）！\n")
