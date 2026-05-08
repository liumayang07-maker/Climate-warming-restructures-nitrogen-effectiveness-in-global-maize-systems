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

# ---- column names ----
col <- list(col_lnRR       <- "lnRR",
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
    lnRR         = suppressWarnings(as.numeric(.data[[col_lnRR]])),
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
    is.finite(lnRR),
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
  filter(country %in% c("USA"))

d_trim <- d_trim %>%
  filter(!is.na(env_group), !is.na(plant_density), is.finite(lnRR)) %>%
  mutate(
    N_strategy  = factor(N_strategy),
    env_group  = droplevels(env_group)
  )

dA <- d_trim %>%
  group_by(env_group) %>%
  mutate(
    q05 = quantile(lnRR, q_lo, na.rm = TRUE),
    q95 = quantile(lnRR, q_hi, na.rm = TRUE),
    keep = lnRR >= q05 & lnRR <= q95
  ) %>%
  ungroup() %>%
  filter(keep)

write.csv(dA, "USA_China_Region_N_data_set.csv")

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
  fixed  = lnRR ~ env_group * plant_density + Nc,
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
    lnRR     = emmean,
    lnRR_lwr = lower.CL,
    lnRR_upr = upper.CL,
    pct     = (exp(lnRR)     - 1) * 100,
    pct_lwr = (exp(lnRR_lwr) - 1) * 100,
    pct_upr = (exp(lnRR_upr) - 1) * 100
  ) %>%
  left_join(count_cell, by = c("env_group","plant_density")) %>%
  mutate(n_lab = paste0("n=", k))

write.csv(emm_df, "USA_lnRR.csv")

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
write.csv(pw1, "USA_tukey_P.csv")
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
  
  br_base_gap = 13.0,
  br_step     = 25.0,
  br_tick     = 8.0,
  br_lwd      = 0.25,
  br_x_pad    = 0.00,
  
  p_text_dy   = 8.0,
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
ymin_plot <- min(ymin_data, cfg$n_y) - cfg$ylim_pad_bottom
ymax_plot <- 250

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
    y = "Relative yield change (%)",
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
  scale_y_continuous(breaks = seq(0, 250, by = 60))

print(p_forest_env_density_pct)

ggsave(
  "USA_forest_pct_by_env_group_density_TukeyP_flip.pdf",
  p_forest_env_density_pct,
  width = 76, height = 38, units = "mm", dpi = 300
)
#####################################################################
# -------------------------
# 1) within-env_group pairwise Wilcoxon (density comparisons)
# -------------------------
pairs_den <- list(
  c(den_levels[1], den_levels[2]),
  c(den_levels[1], den_levels[3]),
  c(den_levels[2], den_levels[3])
)

wilcox_res2 <- dA %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    plant_density = factor(plant_density, levels = den_levels)
  ) %>%
  group_by(env_group) %>%
  group_modify(~{
    df <- .x
    out <- lapply(pairs_den, function(pp){
      a <- df$mean_t[df$plant_density == pp[1]]
      b <- df$mean_t[df$plant_density == pp[2]]
      if (length(a) < 2 || length(b) < 2) {
        return(data.frame(group1 = pp[1], group2 = pp[2], p_value = NA_real_))
      }
      data.frame(
        group1 = pp[1],
        group2 = pp[2],
        p_value = wilcox.test(a, b, exact = FALSE)$p.value
      )
    }) %>% bind_rows()
    
    out$p_adj_bh   <- p.adjust(out$p_value, method = "BH")
    out$p_adj_holm <- p.adjust(out$p_value, method = "holm")
    out
  }) %>%
  ungroup()

print(wilcox_res2)
write.csv(wilcox_res2, "USA_yield_wilcox.csv", row.names = FALSE)

# -------------------------
# 2) letters within each env_group (choose BH or Holm)
# -------------------------
pick_p <- "p_adj_bh"

letters_df2 <- wilcox_res2 %>%
  filter(!is.na(.data[[pick_p]])) %>%
  mutate(comp = paste(group1, group2, sep = "-")) %>%
  group_by(env_group) %>%
  group_modify(~{
    pv <- .x[[pick_p]]
    names(pv) <- .x$comp
    L <- multcompView::multcompLetters(pv)$Letters
    data.frame(plant_density = names(L), Letters = unname(L))
  }) %>%
  ungroup() %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    plant_density = factor(plant_density, levels = den_levels)
  )

# -------------------------
# 3) annotation y positions
# -------------------------
y_pos2 <- dA %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    plant_density = factor(plant_density, levels = den_levels)
  ) %>%
  group_by(env_group, plant_density) %>%
  summarise(y = max(mean_t, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y * 1.03)

ann_df2 <- y_pos2 %>%
  left_join(letters_df2, by = c("env_group","plant_density"))

# -------------------------
# 4) plot: X = env_group, fill = plant_density
# -------------------------
# Optional alternative palette:
# den_cols <- setNames(c("#1B9E77", "#D95F02", "#7570B3"), den_levels)

p_box_letters2 <- ggplot(
  dA %>% mutate(
    env_group     = factor(env_group, levels = env_levels),
    plant_density = factor(plant_density, levels = den_levels)
  ),
  aes(x = env_group, y = mean_t, fill = plant_density)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.5,
    outlier.alpha = 0.25,
    outlier.size  = 0.5,
    outlier.shape = 16,
    outlier.stroke = 0.2,
    outlier.colour = "black"
  ) +
  stat_summary(
    aes(group = plant_density),
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.8),
    shape = 23, fill = NA,
    size = 1.0, stroke = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = den_cols, breaks = den_levels) +
  labs(
    x = "Env group",
    y = "Yield under N fertilization (Mg / ha)",
    fill = "Plant density"
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.75, 0.90),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.width  = grid::unit(3, "mm"),
    legend.key.height = grid::unit(4, "mm"),
    
    axis.line = element_line(color = "black", linewidth = 0.35),
    axis.text.x = element_text(color = "black", hjust = 0.5, size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.35),
    axis.ticks.y = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = grid::unit(1.5, "mm")
  ) +
  coord_cartesian(ylim = c(0, 25), clip = "off") +
  scale_y_continuous(breaks = seq(0, 25, by = 5)) +
  geom_text(
    data = ann_df2,
    aes(x = env_group, y = y, label = Letters, group = plant_density),
    position = position_dodge(width = 0.8),
    vjust = 0, size = 2.5, inherit.aes = FALSE
  )

print(p_box_letters2)

ggsave("box_yield_envXdensity_wilcox_letters.pdf", p_box_letters2,
        width = 151, height = 45, units = "mm", dpi = 300)

################################################################################
# -------------------------
# 0) data
# -------------------------
d_plot <- dA %>%
  filter(is.finite(EDH_sum), is.finite(Effect),
         EDH_sum >= 0, EDH_sum <= 500, Effect >= 50) %>%
  mutate(plant_density = factor(plant_density))

den_levels <- levels(d_plot$plant_density)
den_cols <- setNames(c("#F6C45C", "#F25336", "#66A5ED")[seq_along(den_levels)], den_levels)

# -------------------------
# 1) fit lm per density + extract equation, R2, p
# -------------------------
fit_df <- d_plot %>%
  group_by(plant_density) %>%
  group_modify(~{
    m <- lm(Effect ~ EDH_sum, data = .x)
    sm <- summary(m)
    
    a  <- unname(coef(m)[1])
    b  <- unname(coef(m)[2])
    r2 <- unname(sm$r.squared)
    p  <- unname(sm$coefficients["EDH_sum", "Pr(>|t|)"])
    
    tibble(
      a = a, b = b, r2 = r2, p = p,
      lab = paste0(
        "y = ", sprintf("%.3f", a), " + ", sprintf("%.3f", b), "x, ",
        "R2 = ", sprintf("%.3f", r2), ", P = ",
        ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
      )
    )
  }) %>%
  ungroup()
write.csv(fit_df, "USA_liner_lnRR.csv")
# -------------------------
# 2) annotation positions (stack 3 lines of text)
# -------------------------
x0 <- min(d_plot$EDH_sum, na.rm = TRUE) + 0.02 * diff(range(d_plot$EDH_sum, na.rm = TRUE))
y_top <- max(d_plot$Effect, na.rm = TRUE)
y_rng <- diff(range(d_plot$Effect, na.rm = TRUE))

ann_df <- fit_df %>%
  mutate(
    x = x0,
    y = y_top - (as.numeric(plant_density) - 1) * 0.12 * y_rng  # Adjust spacing if labels overlap
  )

# ---- size-tuned params for 48 x 48 mm ----
cfg <- list(
  base_size   = 7,
  pt_size     = 1.10,
  pt_stroke   = 0.22,
  pt_alpha    = 0.28,
  
  line_lwd    = 0.55,
  rib_alpha   = 0.18,
  
  eq_size     = 2.0,
  eq_lh       = 0.92,
  
  leg_text    = 5.8,
  leg_key_w   = 2.6,   # mm
  leg_key_h   = 2.2,   # mm
  
  tick_len_mm = 1.1,
  axis_lwd    = 0.25,
  tick_lwd    = 0.25
)

p_line <- ggplot(d_plot, aes(x = EDH_sum, y = Effect, color = plant_density, fill = plant_density)) +
  geom_point(shape = 21, fill = "white", alpha = cfg$pt_alpha,
             size = cfg$pt_size, stroke = cfg$pt_stroke) +
  
  # SE ribbon uses fill mapping; line uses color mapping
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              linewidth = cfg$line_lwd, alpha = cfg$rib_alpha) +
  
  geom_text(
    data = ann_df,
    aes(x = x, y = y, label = lab, color = plant_density),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = cfg$eq_size, lineheight = cfg$eq_lh
  ) +
  
  scale_color_manual(values = den_cols, breaks = den_levels) +
  scale_fill_manual(values  = den_cols, breaks = den_levels) +
  
  theme_bw(base_size = cfg$base_size) +
  theme(
    panel.grid = element_blank(),
    
    legend.position = "none",
    
    axis.line  = element_line(color = "black", linewidth = cfg$axis_lwd),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = cfg$tick_lwd),
    axis.ticks.length = grid::unit(cfg$tick_len_mm, "mm"),
    
    plot.margin = margin(2, 2, 2, 2)
  ) +
  labs(x = "Total EDH", y = "Relative yield change (%)") +
  coord_cartesian(ylim = c(0, 300), xlim = c(0, 500), clip = "off") +
  scale_y_continuous(breaks = seq(0, 300, by = 80)) +
  scale_x_continuous(breaks = seq(0, 500, by = 125))

print(p_line)

ggsave("USA_lines_EDH_Effect.pdf", p_line,
       width = 54, height = 48, units = "mm", dpi = 300)



# -------------------------
# AEN
# -------------------------
# -------------------------
# 0) data
# -------------------------
d_plot <- dA %>%
  filter(is.finite(EDH_sum), is.finite(AEN),
         EDH_sum >= 0, EDH_sum <= 500, AEN >= 15) %>%
  mutate(plant_density = factor(plant_density))

den_levels <- levels(d_plot$plant_density)
den_cols <- setNames(c("#F6C45C", "#F25336", "#66A5ED")[seq_along(den_levels)], den_levels)

# -------------------------
# 1) fit lm per density + extract equation, R2, p
# -------------------------
fit_df <- d_plot %>%
  group_by(plant_density) %>%
  group_modify(~{
    m <- lm(AEN ~ EDH_sum, data = .x)
    sm <- summary(m)
    
    a  <- unname(coef(m)[1])
    b  <- unname(coef(m)[2])
    r2 <- unname(sm$r.squared)
    p  <- unname(sm$coefficients["EDH_sum", "Pr(>|t|)"])
    
    tibble(
      a = a, b = b, r2 = r2, p = p,
      lab = paste0(
        "y = ", sprintf("%.3f", a), " + ", sprintf("%.3f", b), "x, ",
        "R2 = ", sprintf("%.3f", r2), ", P = ",
        ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
      )
    )
  }) %>%
  ungroup()
write.csv(fit_df, "USA_liner_AEN.csv")
# -------------------------
# 2) annotation positions (stack 3 lines of text)
# -------------------------
x0 <- min(d_plot$EDH_sum, na.rm = TRUE) + 0.02 * diff(range(d_plot$EDH_sum, na.rm = TRUE))
y_top <- max(d_plot$AEN, na.rm = TRUE)
y_rng <- diff(range(d_plot$AEN, na.rm = TRUE))

ann_df <- fit_df %>%
  mutate(
    x = x0,
    y = y_top - (as.numeric(plant_density) - 1) * 0.12 * y_rng  # Adjust spacing if labels overlap
  )

# ---- size-tuned params for 48 x 48 mm ----
cfg <- list(
  base_size   = 7,
  pt_size     = 1.10,
  pt_stroke   = 0.22,
  pt_alpha    = 0.28,
  
  line_lwd    = 0.55,
  rib_alpha   = 0.18,
  
  eq_size     = 2.0,
  eq_lh       = 0.92,
  
  leg_text    = 5.8,
  leg_key_w   = 2.6,   # mm
  leg_key_h   = 2.2,   # mm
  
  tick_len_mm = 1.1,
  axis_lwd    = 0.25,
  tick_lwd    = 0.25
)

p_line <- ggplot(d_plot, aes(x = EDH_sum, y = AEN, color = plant_density, fill = plant_density)) +
  geom_point(shape = 21, fill = "white", alpha = cfg$pt_alpha,
             size = cfg$pt_size, stroke = cfg$pt_stroke) +
  
  # SE ribbon uses fill mapping; line uses color mapping
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              linewidth = cfg$line_lwd, alpha = cfg$rib_alpha) +
  
  geom_text(
    data = ann_df,
    aes(x = x, y = y, label = lab, color = plant_density),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = cfg$eq_size, lineheight = cfg$eq_lh
  ) +
  
  scale_color_manual(values = den_cols, breaks = den_levels) +
  scale_fill_manual(values  = den_cols, breaks = den_levels) +
  
  theme_bw(base_size = cfg$base_size) +
  theme(
    panel.grid = element_blank(),
    
    legend.position = "none",
    
    axis.line  = element_line(color = "black", linewidth = cfg$axis_lwd),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = cfg$tick_lwd),
    axis.ticks.length = grid::unit(cfg$tick_len_mm, "mm"),
    
    plot.margin = margin(2, 2, 2, 2)
  ) +
  labs(x = "Total EDH", y = "Agronomic efficiency of N (kg grain / kg N)") +
  coord_cartesian(ylim = c(0, 60), xlim = c(0, 500), clip = "off") +
  scale_y_continuous(breaks = seq(0, 60, by = 15)) +
  scale_x_continuous(breaks = seq(0, 500, by = 125))

print(p_line)

ggsave("USA_lines_EDH_AEN.pdf", p_line,
       width = 54, height = 48, units = "mm", dpi = 300)



# -------------------------
# Yield under N fertilization
# -------------------------
# -------------------------
# 0) data
# -------------------------
d_plot <- dA %>%
  filter(is.finite(EDH_sum), is.finite(mean_t),
         EDH_sum >= 0, EDH_sum <= 500) %>%
  mutate(plant_density = factor(plant_density))

den_levels <- levels(d_plot$plant_density)
den_cols <- setNames(c("#F6C45C", "#F25336", "#66A5ED")[seq_along(den_levels)], den_levels)

# -------------------------
# 1) fit lm per density + extract equation, R2, p
# -------------------------
fit_df <- d_plot %>%
  group_by(plant_density) %>%
  group_modify(~{
    m <- lm(mean_t ~ EDH_sum, data = .x)
    sm <- summary(m)
    
    a  <- unname(coef(m)[1])
    b  <- unname(coef(m)[2])
    r2 <- unname(sm$r.squared)
    p  <- unname(sm$coefficients["EDH_sum", "Pr(>|t|)"])
    
    tibble(
      a = a, b = b, r2 = r2, p = p,
      lab = paste0(
        "y = ", sprintf("%.3f", a), " + ", sprintf("%.3f", b), "x, ",
        "R2 = ", sprintf("%.3f", r2), ", P = ",
        ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
      )
    )
  }) %>%
  ungroup()
write.csv(fit_df, "USA_liner_yield.csv")
# -------------------------
# 2) annotation positions (stack 3 lines of text)
# -------------------------
x0 <- min(d_plot$EDH_sum, na.rm = TRUE) + 0.02 * diff(range(d_plot$EDH_sum, na.rm = TRUE))
y_top <- max(d_plot$mean_t, na.rm = TRUE)
y_rng <- diff(range(d_plot$mean_t, na.rm = TRUE))

ann_df <- fit_df %>%
  mutate(
    x = x0,
    y = y_top - (as.numeric(plant_density) - 1) * 0.12 * y_rng  # Adjust spacing if labels overlap
  )

# ---- size-tuned params for 48 x 48 mm ----
cfg <- list(
  base_size   = 7,
  pt_size     = 1.10,
  pt_stroke   = 0.22,
  pt_alpha    = 0.28,
  
  line_lwd    = 0.55,
  rib_alpha   = 0.18,
  
  eq_size     = 2.0,
  eq_lh       = 0.92,
  
  leg_text    = 5.8,
  leg_key_w   = 2.6,   # mm
  leg_key_h   = 2.2,   # mm
  
  tick_len_mm = 1.1,
  axis_lwd    = 0.25,
  tick_lwd    = 0.25
)

p_line <- ggplot(d_plot, aes(x = EDH_sum, y = mean_t, color = plant_density, fill = plant_density)) +
  geom_point(shape = 21, fill = "white", alpha = cfg$pt_alpha,
             size = cfg$pt_size, stroke = cfg$pt_stroke) +
  
  # SE ribbon uses fill mapping; line uses color mapping
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              linewidth = cfg$line_lwd, alpha = cfg$rib_alpha) +
  
  geom_text(
    data = ann_df,
    aes(x = x, y = y, label = lab, color = plant_density),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = cfg$eq_size, lineheight = cfg$eq_lh
  ) +
  
  scale_color_manual(values = den_cols, breaks = den_levels) +
  scale_fill_manual(values  = den_cols, breaks = den_levels) +
  
  theme_bw(base_size = cfg$base_size) +
  theme(
    panel.grid = element_blank(),
    
    legend.position = "none",
    
    axis.line  = element_line(color = "black", linewidth = cfg$axis_lwd),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = cfg$tick_lwd),
    axis.ticks.length = grid::unit(cfg$tick_len_mm, "mm"),
    
    plot.margin = margin(2, 2, 2, 2)
  ) +
  labs(x = "Total EDH", y = "Yield under N fertilization (Mg / ha)") +
  coord_cartesian(ylim = c(0, 20), xlim = c(0, 500), clip = "off") +
  scale_y_continuous(breaks = seq(0, 20, by = 5)) +
  scale_x_continuous(breaks = seq(0, 500, by = 125))

print(p_line)

ggsave("USA_lines_EDH_mean_t.pdf", p_line,
       width = 54, height = 48, units = "mm", dpi = 300)


