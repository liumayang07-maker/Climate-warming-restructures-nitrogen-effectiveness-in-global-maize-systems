# =========================================================
# 0) Packages & data
# =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(nlme)
  library(emmeans)
  library(multcomp)
  library(multcompView)
  library(splines)
})

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

subset_countries <- c("USA")
output_prefix <- "usa"

out_file <- function(name) {
  file.path(output_dir, paste0(output_prefix, "_", name))
}

# ---- column names ----
col <- list(col_AEN       <- "AEN",
            col_env        <- "Env_Type",        # Normal/PFR/GFR
            col_ref        <- "reference_ID",
            col_sy         <- "site_year_ID",
            col_continent  <- "continent",
            col_country    <- "country",
            col_mean_c     <- "mean_c",
            col_mean_t     <- "mean_t",
            col_EDH_sum    <- "EDH35_total",
            col_EDH_label  <- "EDH35_label",
            col_N_rate     <- "N_rate",
            col_N_rate_cluster   <- "N_rate_cluster",
            col_env_ID     <- "env_ID",
            col_N0         <- "Label_of_control_treatment",
            col_N_strategy <- "N_strategy",
            col_N_S        <- "N_S",
            col_N_V        <- "N_V",
            col_N_R        <- "N_R",
            col_plant_density <- "plant_density_3cluster",
            col_HY_strategy2 <- "HY_strategy2",
            col_region_N_label <- "region_N_label"
)

env_levels <- c("Normal","PFR","GFR")

# ---- trim settings ----
q_lo <- 0.05
q_hi <- 0.95

# ---- plotting palette (consistent) ----
env_cols <- c(Normal = "#1b9e77", PFR = "#d95f02", GFR = "#7570b3")

# =========================================================
# 1) Data prep (row-level) + centering
# =========================================================
d0 <- dat %>%
  mutate(
    AEN         = suppressWarnings(as.numeric(.data[[col_AEN]])),
    env_group    = factor(.data[[col_env]], levels = env_levels),
    reference_ID = as.factor(.data[[col_ref]]),
    site_year_ID = as.factor(.data[[col_sy]]),
    continent    = as.factor(.data[[col_continent]]),
    country      = as.factor(.data[[col_country]]),
    N0           = as.factor(.data[[col_N0]]),
    N_strategy   = as.factor(.data[[col_N_strategy]]),
    plant_density = as.factor(.data[[col_plant_density]]),
    HY_strategy2  = as.factor(.data[[col_HY_strategy2]]),
    N_rate_cluster = as.factor(.data[[col_N_rate_cluster]]),
    EDH_label    = as.factor(.data[[col_EDH_label]]),
    region_N_label    = as.factor(.data[[col_region_N_label]]),
    mean_c       = suppressWarnings(as.numeric(.data[[col_mean_c]])),
    mean_t       = suppressWarnings(as.numeric(.data[[col_mean_t]])),
    EDH_sum      = suppressWarnings(as.numeric(.data[[col_EDH_sum]])),
    N_rate       = suppressWarnings(as.numeric(.data[[col_N_rate]])),
    N_S          = suppressWarnings(as.numeric(.data[[col_N_S]])),
    N_V          = suppressWarnings(as.numeric(.data[[col_N_V]])),
    N_R          = suppressWarnings(as.numeric(.data[[col_N_R]]))
  ) %>%
  filter(
    is.finite(AEN),
    !is.na(env_group),
    !is.na(reference_ID),
    !is.na(site_year_ID),
    !is.na(continent),
    is.finite(N_rate)
  ) %>%
  mutate(
    Nc = N_rate - mean(N_rate, na.rm = TRUE)  # global centering (for linear sensitivity)
  )


d_trim <- d0 %>%
  filter(region_N_label %in% c("RN")) %>%
  filter(country %in% subset_countries)

d_trim <- d_trim %>%
  filter(!is.na(env_group), !is.na(plant_density), is.finite(AEN)) %>%
  mutate(
    N_strategy  = factor(N_strategy),
    env_group  = droplevels(env_group)
  )

dA <- d_trim %>%
  group_by(env_group) %>%
  mutate(
    q05 = quantile(AEN, q_lo, na.rm = TRUE),
    q95 = quantile(AEN, q_hi, na.rm = TRUE),
    keep = AEN >= q05 & AEN <= q95
  ) %>%
  ungroup() %>%
  filter(keep)

# cell counts (for labels)
count_cell <- dA %>%
  group_by(plant_density, env_group) %>%
  summarise(
    k     = dplyr::n(),
    n_ref = n_distinct(reference_ID),
    n_sy  = n_distinct(site_year_ID),
    .groups = "drop"
  ) %>%
  mutate(n_lab = paste0("n=", k))

print(count_cell)


m0 <- nlme::lme(
  fixed  = AEN ~ env_group * plant_density + Nc,
  random = ~ 1 | reference_ID/site_year_ID,
  data   = dA,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

anova(m0)

# =========================================================
# emmeans (within env_group) + Tukey pairwise p-values
# Plot: Y=env_group, X=% effect; show raw scatter distribution
# =========================================================
env_levels <- levels(factor(d_trim$env_group))
den_levels <- levels(factor(d_trim$plant_density))

# -------------------------
# 1) emmeans: within env_group, compare plant_density
# -------------------------
emm1 <- emmeans(m0, ~ plant_density | env_group)

emm_df <- as.data.frame(emm1) %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    plant_density = factor(plant_density, levels = den_levels),
    AEN     = emmean,
    AEN_lwr = lower.CL,
    AEN_upr = upper.CL,
    pct     = AEN,
    pct_lwr = AEN_lwr,
    pct_upr = AEN_upr
  ) %>%
  left_join(count_cell, by = c("env_group","plant_density")) %>%
  mutate(n_lab = paste0("n=", k))
write.csv(emm_df, out_file("aen.csv"))
# -------------------------
# 2) Tukey pairwise within env_group
# -------------------------
pw1 <- as.data.frame(pairs(emm1, adjust = "tukey")) %>%
  tidyr::separate(contrast, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    env_group = factor(env_group, levels = env_levels),
    g1 = factor(g1, levels = den_levels),
    g2 = factor(g2, levels = den_levels),
    p_lab = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )
write.csv(pw1, out_file("aen_tukey_p.csv"))
# =========================
# 3) plotting cfg tuned for 76 x 43 mm
# =========================
cfg <- list(
  dodge_w   = 0.62,
  n_groups  = NULL,
  
  n_y       = -15,     # (pre-flip y) -> becomes horizontal axis after flip
  n_angle   = 0,
  n_size    = 1.9,
  n_vjust   = 0.9,
  
  br_base_gap = 5.0,
  br_step     = 8.0,
  br_tick     = 2.8,
  br_lwd      = 0.25,
  br_x_pad    = 0.00,
  
  p_text_dy   = 2.5,
  p_size      = 1.9,
  
  pt_size     = 1.35,
  pt_stroke   = 0.35,
  err_lwd     = 0.25,
  err_width   = 0.28,
  
  base_size   = 5.0,
  
  ylim_pad_top    = 2.0,
  ylim_pad_bottom = 1.5
)
if (is.null(cfg$n_groups)) cfg$n_groups <- length(den_levels)

# -------------------------
# 4) colors for density
# -------------------------
den_cols <- setNames(c("#F6C45C", "#F25336", "#66A5ED"), den_levels)

# -------------------------
# 5) dodge offsets (align with position_dodge)
# -------------------------
den_offset <- tibble(
  plant_density = factor(den_levels, levels = den_levels),
  idx = seq_along(den_levels)
) %>%
  mutate(
    x_off = (idx - (cfg$n_groups + 1)/2) * (cfg$dodge_w / cfg$n_groups)
  ) %>%
  dplyr::select(plant_density, x_off)

x_map <- tibble(
  env_group = factor(env_levels, levels = env_levels),
  x_base = seq_along(env_levels)
)

# -------------------------
# 6) bracket data (pct scale)
# -------------------------
y_top <- emm_df %>%
  group_by(env_group) %>%
  summarise(ymax = max(pct_upr, na.rm = TRUE), .groups = "drop")

br_df <- pw1 %>%
  left_join(x_map, by = "env_group") %>%
  left_join(den_offset, by = c("g1" = "plant_density")) %>%
  rename(x1_off = x_off) %>%
  left_join(den_offset, by = c("g2" = "plant_density")) %>%
  rename(x2_off = x_off) %>%
  left_join(y_top, by = "env_group") %>%
  group_by(env_group) %>%
  mutate(
    br_rank = row_number(),
    y    = ymax + cfg$br_base_gap + (br_rank - 1) * cfg$br_step,
    x1   = x_base + x1_off - cfg$br_x_pad,
    x2   = x_base + x2_off + cfg$br_x_pad,
    xmid = (x1 + x2) / 2,
    y_txt = y + cfg$p_text_dy
  ) %>%
  ungroup()

# -------------------------
# 7) y-limits (pct scale) for post-flip horizontal axis
# -------------------------
ymin_data <- min(emm_df$pct_lwr, na.rm = TRUE)
ymin_plot <- -5
ymax_plot <- 60

# =========================
# 8) plot (forest only + flip)
# =========================
p_forest_env_density_pct <- ggplot(emm_df, aes(x = env_group, y = pct, color = plant_density)) +
  
  geom_errorbar(
    aes(ymin = pct_lwr, ymax = pct_upr),
    position = position_dodge(width = cfg$dodge_w),
    width = cfg$err_width,
    linewidth = cfg$err_lwd
  ) +
  geom_point(
    position = position_dodge(width = cfg$dodge_w),
    shape = 21, fill = "white",
    size = cfg$pt_size, stroke = cfg$pt_stroke
  ) +
  
  # --- n labels (appear on left after flip) ---
  geom_text(
    data = emm_df,
    aes(x = env_group, y = 0, label = n_lab, group = plant_density),
    position = position_dodge(width = cfg$dodge_w),
    color = "black",
    size = cfg$n_size,
    vjust = cfg$n_vjust,
    inherit.aes = FALSE,
    angle = 45
  ) +
  
  # --- brackets ---
  geom_segment(
    data = br_df,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    color = "black", linewidth = cfg$br_lwd
  ) +
  geom_segment(
    data = br_df,
    aes(x = x1, xend = x1, y = y, yend = y - cfg$br_tick),
    inherit.aes = FALSE,
    color = "black", linewidth = cfg$br_lwd
  ) +
  geom_segment(
    data = br_df,
    aes(x = x2, xend = x2, y = y, yend = y - cfg$br_tick),
    inherit.aes = FALSE,
    color = "black", linewidth = cfg$br_lwd
  ) +
  geom_text(
    data = br_df,
    aes(x = xmid, y = y_txt, label = p_lab),
    inherit.aes = FALSE,
    color = "black", size = cfg$p_size
  ) +
  
  scale_color_manual(values = den_cols, breaks = den_levels) +
  labs(
    x = "Env group",
    y = "Agronomic efficiency of N (kg grain / kg N)",
    color = NULL
  ) +
  theme_bw(base_size = cfg$base_size) +
  theme(
    panel.grid = element_blank(),
    
    legend.position = c(0.70, 0.94),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 5.0),
    legend.key.width  = grid::unit(3.2, "mm"),
    legend.key.height = grid::unit(1.5, "mm"),
    
    axis.line = element_line(color = "black", linewidth = 0.30),
    axis.text.x = element_text(color = "black", size = 6.0),
    axis.text.y = element_text(color = "black", size = 6.0),
    axis.title.x = element_text(size = 6.0),
    axis.title.y = element_text(size = 6.0),
    
    axis.ticks.x = element_line(color = "black", linewidth = 0.30),
    axis.ticks.y = element_line(color = "black", linewidth = 0.30),
    axis.ticks.length = grid::unit(1.2, "mm"),
    
    plot.margin = margin(2, 10, 2, 4)  # room for right p labels after flip
  ) +
  coord_cartesian(ylim = c(ymin_plot, ymax_plot), clip = "off") +
  ## y-axis ticks (fixed)
  scale_y_continuous(breaks = seq(0, 60, by = 20))

print(p_forest_env_density_pct)

ggsave(
  out_file("aen_forest_pct_by_env_group_density_tukeyp_flip.pdf"),
  p_forest_env_density_pct,
  width = 76, height = 38, units = "mm", dpi = 300
)
