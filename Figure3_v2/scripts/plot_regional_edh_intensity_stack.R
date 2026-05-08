library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

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
setwd(output_dir)
dat <- read.csv(file.path(input_dir, "meta_data_v2.csv"), stringsAsFactors = FALSE)

# ---- column names ----
col_lnRR       <- "lnRR"
col_env        <- "Env_Type"        # Normal/PFR/GFR
col_ref        <- "reference_ID"
col_sy         <- "site_year_ID"
col_continent  <- "continent"
col_region     <- "region"
col_mean_c     <- "mean_c"
col_mean_t     <- "mean_t"
col_EDH35_sum    <- "EDH35_total"
col_EDH35_label  <- "EDH35_label"
col_N_rate     <- "N_rate"
col_exp_duration <- "Exp_year_duration"
col_env_ID      <- "env_ID"

env_levels <- c("Normal","PFR","GFR")
env_use <- c("PFR", "GFR")
EDH35_levels <- c("None", "PFR_L", "PFR_M", "PFR_H", "GFR_L", "GFR_M", "GFR_H")

# =========================
# Public-release note
# =========================
d_reg <- dat %>%
  mutate(
    region      = as.character(.data[[col_region]]),
    env_ID      = .data[[col_env_ID]],
    EDH35_label = factor(.data[[col_EDH35_label]], levels = EDH35_levels)
  ) %>%
  filter(!is.na(region), region != "", region != "South China",
         !is.na(env_ID),
         !is.na(EDH35_label)) %>%
  dplyr::select(env_ID, region, EDH35_label) %>%
  distinct(env_ID, .keep_all = TRUE)

# =========================
# Public-release note
# =========================
prop_reg <- d_reg %>%
  count(region, EDH35_label, name = "n") %>%
  complete(region, EDH35_label, fill = list(n = 0)) %>%
  group_by(region) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# =========================
# Public-release note
# =========================
region_order <- d_reg %>%
  count(region, name = "n_env") %>%
  arrange(desc(n_env)) %>%
  pull(region)

prop_reg <- prop_reg %>%
  mutate(region = factor(region, levels = region_order))

# =========================
# Public-release note
# Public-release note
# =========================
edd_cols <- c(
  "None"  = "#ffd7c2",
  "PFR_L" = "#e0e8f2",
  "PFR_M" = "#8b9bc1",
  "PFR_H" = "#35568a",
  "GFR_L" = "#cde3d3",
  "GFR_M" = "#6eb169",
  "GFR_H" = "#397c52"
)

# Public-release note
stopifnot(all(levels(prop_reg$EDH35_label) %in% names(edd_cols)))

write.csv(prop_reg, "EDH_intensity_stack_Figure.csv")

reg_n <- d_reg %>% count(region, name = "n_env")

# =========================
# Public-release note
# =========================
# ---- size tuned for 39 x 60 mm ----
base_size <- 6.4
axis_title_fs <- 6.2
tick_fs <- 6.0
ytext_fs <- 5.8
legend_fs <- 5.8

p_reg_stack_flip <- ggplot(prop_reg, aes(y = region, x = prop, fill = EDH35_label)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = edd_cols, breaks = EDH35_levels, drop = FALSE) +
  scale_x_continuous(
    limits = c(0, 1.10),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = percent_format(accuracy = 1)
  ) +
  labs(x = "Proportion (%)", y = NULL, fill = "Total EDH intensity") +
  geom_text(
    data = reg_n %>% mutate(region = factor(region, levels = region_order)),
    aes(y = region, x = 1.08, label = paste0("n=", n_env)),
    inherit.aes = FALSE,
    size = 2.6
  ) +
  theme_bw(base_size = base_size) +
  theme(
    panel.grid = element_blank(),
    
    axis.title.x = element_text(size = axis_title_fs),
    axis.text.x  = element_text(size = tick_fs),
    axis.text.y  = element_text(size = ytext_fs),
    
    # Public-release note
    legend.position = "bottom",
    legend.title = element_text(size = legend_fs),
    legend.text  = element_text(size = legend_fs),
    legend.key.height = unit(2.6, "mm"),
    legend.key.width  = unit(4.0, "mm"),
    legend.spacing.x  = unit(1.2, "mm"),
    
    # Public-release note
    plot.margin = margin(t = 1.5, r = 10, b = 1.5, l = 1.5, unit = "mm")
  ) +
  coord_cartesian(xlim = c(0, 1.10), clip = "off")

p_reg_stack_flip

ggsave(
 "region_EDH_intensity.pdf",
 p_reg_stack_flip, width = 60, height = 60, units = "mm",
  dpi = 300, useDingbats = FALSE
)
