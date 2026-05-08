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
dat <- read.csv(file.path(input_dir, "region_meta_data_v2.csv"), stringsAsFactors = FALSE)

# ---- column names ----
col_AEN       <- "AEN"
col_env        <- "Env_Type"        # Normal/PFR/GFR
col_ref        <- "reference_ID"
col_sy         <- "site_year_ID"
col_continent  <- "continent"
col_country    <- "country"
col_mean_c     <- "mean_c"
col_mean_t     <- "mean_t"
col_EDH_sum    <- "EDH35_total"
col_EDH_label  <- "EDH35_label"
col_N_rate     <- "N_rate"
col_env_ID     <- "env_ID"
col_N0         <- "Label_of_control_treatment"

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
    mean_c       = suppressWarnings(as.numeric(.data[[col_mean_c]])),
    mean_t       = suppressWarnings(as.numeric(.data[[col_mean_t]])),
    EDH_sum      = suppressWarnings(as.numeric(.data[[col_EDH_sum]])),
    N_rate       = suppressWarnings(as.numeric(.data[[col_N_rate]]))
  ) %>%
  filter(
    is.finite(AEN),
    !is.na(env_group),
    !is.na(reference_ID),
    !is.na(site_year_ID),
    !is.na(continent),
    is.finite(N_rate),
  ) %>%
  mutate(
    Nc = N_rate - mean(N_rate, na.rm = TRUE)  # global centering (for linear sensitivity)
  )

d_trim <- d0 %>%
   group_by(env_group) %>%
   mutate(
    q05 = quantile(AEN, q_lo, na.rm = TRUE),
    q95 = quantile(AEN, q_hi, na.rm = TRUE),
    keep = AEN >= q05 & AEN <= q95
   ) %>%
   ungroup() %>%
   filter(keep)


# cell counts (for labels)
count_cell <- d_trim %>%
  group_by(continent, env_group) %>%
  summarise(
    k     = dplyr::n(),
    n_ref = n_distinct(reference_ID),
    n_sy  = n_distinct(site_year_ID),
    .groups = "drop"
  ) %>%
  mutate(n_lab = paste0("n=", k))

print(count_cell)

# =========================================================
# 4) Sensitivity: linear vs spline control for N_rate (ML)
# =========================================================
m_lin <- nlme::lme(
  fixed  = AEN ~ env_group * continent + Nc,
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d_trim,
  method = "ML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

m_ns3 <- nlme::lme(
  fixed  = AEN ~ env_group * continent + ns(N_rate, df = 3),
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d_trim,
  method = "ML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

m_ns4 <- nlme::lme(
  fixed  = AEN ~ env_group * continent + ns(N_rate, df = 4),
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d_trim,
  method = "ML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

print(AIC(m_lin, m_ns3, m_ns4))
print(anova(m_lin, m_ns3, m_ns4))  # Public-release note


# =========================================================
# 5) Main-text model choice: m0 (shared spline curve)
#    Main text focus: Env differences + continent differences
# =========================================================
m0 <- nlme::lme(
  fixed  = AEN ~ env_group * continent + ns(N_rate, df = 3),
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d_trim,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

print(anova(m0))

# =========================================================
# emmeans (within continent) + Tukey pairwise p-values
# =========================================================
# emmeans: within each continent, 3 env means
emm0 <- emmeans(m0, ~ env_group | continent)

# Public-release note
emm_df <- as.data.frame(emm0) %>%
  mutate(
    pct     = emmean,
    pct_lwr = lower.CL,
    pct_upr = upper.CL,
    continent = factor(continent, levels = levels(d_trim$continent)),
    env_group = factor(env_group, levels = env_levels)
  ) %>%
  left_join(count_cell, by = c("continent","env_group")) %>%
  mutate(
    n_lab = paste0("n=", k)
  )
write.csv(emm_df, "Data_Supplement Figure3-2B.csv")
# Public-release note
n_offset <- 0.06 * (max(emm_df$pct_upr, na.rm = TRUE) - min(emm_df$pct_lwr, na.rm = TRUE))
emm_df <- emm_df %>% mutate(y_n = pct_lwr - n_offset)

# Public-release note
pw0 <- as.data.frame(pairs(emm0, adjust = "tukey")) %>%
  tidyr::separate(contrast, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    continent = factor(continent, levels = levels(d_trim$continent)),
    g1 = factor(g1, levels = env_levels),
    g2 = factor(g2, levels = env_levels),
    p_lab = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )
write.csv(pw0, "Data_Supplement Figure3-2B-1.csv")
# =========================
# Public-release note
# =========================
cfg <- list(
  # Public-release note
  dodge_w   = 0.72,  # Public-release note
  n_groups  = NULL,  # Public-release note
  
  # Public-release note
  n_y       = 10,  # Public-release note
  n_angle   = 45,
  n_size    = 2.6,
  n_vjust   = 1,
  
  # Public-release note
  br_base_gap = 5.0,  # Public-release note
  br_step     = 6.0,  # Public-release note
  br_tick     = 1.5,  # Public-release note
  br_lwd      = 0.25,  # Public-release note
  
  # Public-release note
  br_x_pad    = 0.00,  # Public-release note
  
  # Public-release note
  p_text_dy   = 2.0,  # Public-release note
  p_size      = 2.0,
  
  # Public-release note
  pt_size     = 1.5,
  pt_stroke   = 0.4,
  err_lwd     = 0.25,
  err_width   = 0.32,  # Public-release note
  
  # Public-release note
  base_size   = 8,
  axis_x_angle = 0,
  
  # Public-release note
  ylim_pad_top    = 3.5,
  ylim_pad_bottom = 2.0
)

# =========================
# Public-release note
# =========================
if (is.null(cfg$n_groups)) cfg$n_groups <- length(env_levels)

# Public-release note
env_offset <- tibble(
  env_group = factor(env_levels, levels = env_levels),
  idx = seq_along(env_levels)
) %>%
  mutate(
    x_off = (idx - (cfg$n_groups + 1)/2) * (cfg$dodge_w / cfg$n_groups)
  ) %>%
  dplyr::select(env_group, x_off)

# Public-release note
cont_lv <- levels(factor(emm_df$continent))
cont_map <- tibble(
  continent = factor(cont_lv, levels = cont_lv),
  x_base = seq_along(cont_lv)
)

# =========================
# Public-release note
# =========================
y_top <- emm_df %>%
  group_by(continent) %>%
  summarise(ymax = max(pct_upr, na.rm = TRUE), .groups = "drop")

br_df <- pw0 %>%
  left_join(cont_map, by = "continent") %>%
  left_join(env_offset, by = c("g1" = "env_group")) %>%
  rename(x1_off = x_off) %>%
  left_join(env_offset, by = c("g2" = "env_group")) %>%
  rename(x2_off = x_off) %>%
  left_join(y_top, by = "continent") %>%
  group_by(continent) %>%
  mutate(
    br_rank = row_number(),  # Public-release note
    y   = ymax + cfg$br_base_gap + (br_rank - 1) * cfg$br_step,
    x1  = x_base + x1_off - cfg$br_x_pad,
    x2  = x_base + x2_off + cfg$br_x_pad,
    xmid = (x1 + x2) / 2,
    y_txt = y + cfg$p_text_dy
  ) %>%
  ungroup()

# =========================
# Public-release note
# =========================
ymin_data <- min(emm_df$pct_lwr, na.rm = TRUE)
ymin_n    <- cfg$n_y  # Public-release note
ymin_plot <- min(ymin_data, ymin_n) - cfg$ylim_pad_bottom

ymax_plot <- 60

# =========================
# Public-release note
# =========================
p_forest_continent <- ggplot(emm_df, aes(x = continent, y = pct, color = env_group)) +
  
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
  
  # --- n labels ---
  geom_text(
    data = emm_df,
    aes(x = continent, y = 5, label = n_lab, group = env_group, angle = cfg$n_angle),
    position = position_dodge(width = cfg$dodge_w),
    color = "black",
    size = cfg$n_size,
    vjust = cfg$n_vjust,
    inherit.aes = FALSE
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
  
  scale_color_manual(values = env_cols, breaks = env_levels) +
  labs(
    x = "continent",
    y = "Agronomic efficiency of N (kg grain / kg N)",
    color = "Env"
  ) +
  theme_bw(base_size = cfg$base_size) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.70, 0.90),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.line = element_line(color = "black", linewidth = 0.35),
    # Public-release note
    axis.text.x = element_text(color = "black", hjust = 0.5, size = 8),  # Public-release note
    axis.text.y = element_text(color = "black", size = 8),
    
    # Public-release note
    axis.ticks.x = element_line(color = "black", linewidth = 0.35),
    axis.ticks.y = element_line(color = "black", linewidth = 0.35),
    
    # Public-release note
    axis.ticks.length = grid::unit(1.5, "mm")
    
  ) +
  coord_cartesian(ylim = c(ymin_plot, ymax_plot), clip = "off") +
## y-axis ticks (fixed)
  scale_y_continuous(breaks = seq(0, 60, by = 15))

print(p_forest_continent)

ggsave(
  "AEN_region_N_forest_by_continent_env_tukeyP.pdf",
  p_forest_continent,
  width = 110, height = 55, units = "mm", dpi = 300
)


##### Public-release note
# --- within-continent pairwise Wilcoxon ---
pairs_env <- list(c("Normal","PFR"), c("Normal","GFR"), c("PFR","GFR"))

wilcox_res <- d_trim %>%
  group_by(continent) %>%
  group_modify(~{
    df <- .x
    out <- lapply(pairs_env, function(pp){
      a <- df$mean_t[df$env_group == pp[1]]
      b <- df$mean_t[df$env_group == pp[2]]
      if (length(a) < 2 || length(b) < 2) {
        return(data.frame(group1=pp[1], group2=pp[2], p_value=NA_real_))
      }
      data.frame(group1=pp[1], group2=pp[2],
                 p_value = wilcox.test(a, b, exact = FALSE)$p.value)
    }) %>% bind_rows()
    
    out$p_adj_bh   <- p.adjust(out$p_value, method = "BH")
    out
  }) %>%
  ungroup()

print(wilcox_res)

# --- letters (choose BH or Holm) ---
pick_p <- "p_adj_bh"

letters_df <- wilcox_res %>%
  filter(!is.na(.data[[pick_p]])) %>%
  mutate(comp = paste(group1, group2, sep = "-")) %>%
  group_by(continent) %>%
  group_modify(~{
    pv <- .x[[pick_p]]
    names(pv) <- .x$comp
    L <- multcompView::multcompLetters(pv)$Letters
    data.frame(env_group = names(L), Letters = unname(L))
  }) %>%
  ungroup() %>%
  mutate(env_group = factor(env_group, levels = env_levels))

y_pos <- d_trim %>%
  group_by(continent, env_group) %>%
  summarise(y = max(mean_t, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y * 1.03)

ann_df <- y_pos %>%
  left_join(letters_df, by = c("continent","env_group"))

p_box_letters <- ggplot(d_trim, aes(x = continent, y = mean_t, fill = env_group)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.5,
    outlier.alpha = 0.25,
    outlier.size  = 0.5,  # Public-release note
    outlier.shape = 16,  # Public-release note
    outlier.stroke = 0.2,  # Public-release note
    outlier.colour = "black"
  ) +
  stat_summary(
    aes(group = env_group),  # Public-release note
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.8),
    shape = 23, fill = NA,
    size = 1.0, stroke = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = env_cols, breaks = env_levels) +
  labs(x = "continent", y = "Yield under N fertilization (Mg / ha)", fill = "Env") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.80, 0.90),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.width  = grid::unit(3, "mm"),
        legend.key.height = grid::unit(4, "mm"),
        axis.line = element_line(color = "black", linewidth = 0.35),
        # Public-release note
        axis.text.x = element_text(color = "black", hjust = 0.5, size = 8),  # Public-release note
        axis.text.y = element_text(color = "black", size = 8),
        
        # Public-release note
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.y = element_line(color = "black", linewidth = 0.35),
        
        # Public-release note
        axis.ticks.length = grid::unit(1.5, "mm")
  ) +
  coord_cartesian(ylim = c(0, 25), clip = "off") +
  ## y-axis ticks (fixed)
  scale_y_continuous(breaks = seq(0, 25, by = 5)) +
  geom_text(
    data = ann_df,
    aes(x = continent, y = y, label = Letters, group = env_group),
    position = position_dodge(width = 0.8),
    vjust = 0, size = 2.5, inherit.aes = FALSE
  )

print(p_box_letters)

ggsave("region_speical_N_box_mean_t_by_continent_env_wilcox_letters.pdf", p_box_letters,
       width = 110, height = 55, units = "mm", dpi = 300)

################# Public-release note
# --- within-continent pairwise Wilcoxon ---
pairs_env <- list(c("Normal","PFR"), c("Normal","GFR"), c("PFR","GFR"))

wilcox_res <- d_trim %>%
  group_by(continent) %>%
  group_modify(~{
    df <- .x
    out <- lapply(pairs_env, function(pp){
      a <- df$N_rate[df$env_group == pp[1]]
      b <- df$N_rate[df$env_group == pp[2]]
      if (length(a) < 2 || length(b) < 2) {
        return(data.frame(group1=pp[1], group2=pp[2], p_value=NA_real_))
      }
      data.frame(group1=pp[1], group2=pp[2],
                 p_value = wilcox.test(a, b, exact = FALSE)$p.value)
    }) %>% bind_rows()
    
    out$p_adj_bh   <- p.adjust(out$p_value, method = "BH")
    out
  }) %>%
  ungroup()

print(wilcox_res)

# --- letters (choose BH or Holm) ---
pick_p <- "p_adj_bh"

letters_df <- wilcox_res %>%
  filter(!is.na(.data[[pick_p]])) %>%
  mutate(comp = paste(group1, group2, sep = "-")) %>%
  group_by(continent) %>%
  group_modify(~{
    pv <- .x[[pick_p]]
    names(pv) <- .x$comp
    L <- multcompView::multcompLetters(pv)$Letters
    data.frame(env_group = names(L), Letters = unname(L))
  }) %>%
  ungroup() %>%
  mutate(env_group = factor(env_group, levels = env_levels))

y_pos <- d_trim %>%
  group_by(continent, env_group) %>%
  summarise(y = max(N_rate, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y * 1.03)

ann_df <- y_pos %>%
  left_join(letters_df, by = c("continent","env_group"))

p_box_letters <- ggplot(d_trim, aes(x = continent, y = N_rate, fill = env_group)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.5,
    outlier.alpha = 0.25,
    outlier.size  = 0.5,  # Public-release note
    outlier.shape = 16,  # Public-release note
    outlier.stroke = 0.2,  # Public-release note
    outlier.colour = "black"
  ) +
  stat_summary(
    aes(group = env_group),  # Public-release note
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.8),
    shape = 23, fill = NA,
    size = 1.0, stroke = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = env_cols, breaks = env_levels) +
  labs(x = "continent", y = "N application rate (kg / ha)", fill = "Env") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.80, 0.90),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.width  = grid::unit(3, "mm"),
        legend.key.height = grid::unit(4, "mm"),
        axis.line = element_line(color = "black", linewidth = 0.35),
        # Public-release note
        axis.text.x = element_text(color = "black", hjust = 0.5, size = 8),  # Public-release note
        axis.text.y = element_text(color = "black", size = 8),
        
        # Public-release note
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.y = element_line(color = "black", linewidth = 0.35),
        
        # Public-release note
        axis.ticks.length = grid::unit(1.5, "mm")
  ) +
  coord_cartesian(ylim = c(0, 600), clip = "off") +
  ## y-axis ticks (fixed)
  scale_y_continuous(breaks = seq(0, 600, by = 125)) 
print(p_box_letters)

ggsave("region_speical_N_box_N_rate_by_continent_env_wilcox_letters.pdf", p_box_letters,
       width = 110, height = 55, units = "mm", dpi = 300)

##### Public-release note
# --- within-continent pairwise Wilcoxon ---
pairs_env <- list(c("Normal","PFR"), c("Normal","GFR"), c("PFR","GFR"))

d_trim <- d_trim %>%
  mutate(.row_id = row_number()) %>%
  group_by(env_ID) %>%
  arrange(.row_id, .by_group = TRUE) %>%
  filter(is.na(lag(mean_c)) | mean_c != lag(mean_c)) %>%
  ungroup() %>%
  dplyr::select(-.row_id)


wilcox_res <- d_trim %>%
  group_by(continent) %>%
  group_modify(~{
    df <- .x
    out <- lapply(pairs_env, function(pp){
      a <- df$mean_c[df$env_group == pp[1]]
      b <- df$mean_c[df$env_group == pp[2]]
      if (length(a) < 2 || length(b) < 2) {
        return(data.frame(group1=pp[1], group2=pp[2], p_value=NA_real_))
      }
      data.frame(group1=pp[1], group2=pp[2],
                 p_value = wilcox.test(a, b, exact = FALSE)$p.value)
    }) %>% bind_rows()
    
    out$p_adj_bh   <- p.adjust(out$p_value, method = "BH")
    out
  }) %>%
  ungroup()

print(wilcox_res)

# --- letters (choose BH or Holm) ---
pick_p <- "p_adj_bh"

letters_df <- wilcox_res %>%
  filter(!is.na(.data[[pick_p]])) %>%
  mutate(comp = paste(group1, group2, sep = "-")) %>%
  group_by(continent) %>%
  group_modify(~{
    pv <- .x[[pick_p]]
    names(pv) <- .x$comp
    L <- multcompView::multcompLetters(pv)$Letters
    data.frame(env_group = names(L), Letters = unname(L))
  }) %>%
  ungroup() %>%
  mutate(env_group = factor(env_group, levels = env_levels))

y_pos <- d_trim %>%
  group_by(continent, env_group) %>%
  summarise(y = max(mean_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y * 1.03)

ann_df <- y_pos %>%
  left_join(letters_df, by = c("continent","env_group"))

p_box_letters <- ggplot(d_trim, aes(x = continent, y = mean_c, fill = env_group)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.5,
    outlier.alpha = 0.25,
    outlier.size  = 0.5,  # Public-release note
    outlier.shape = 16,  # Public-release note
    outlier.stroke = 0.2,  # Public-release note
    outlier.colour = "black"
  ) +
  stat_summary(
    aes(group = env_group),  # Public-release note
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.8),
    shape = 23, fill = NA,
    size = 1.0, stroke = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = env_cols, breaks = env_levels) +
  labs(x = "continent", y = "Yield under N0 level (Mg / ha)", fill = "Env") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.80, 0.90),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.width  = grid::unit(3, "mm"),
        legend.key.height = grid::unit(4, "mm"),
        axis.line = element_line(color = "black", linewidth = 0.35),
        # Public-release note
        axis.text.x = element_text(color = "black", hjust = 0.5, size = 8),  # Public-release note
        axis.text.y = element_text(color = "black", size = 8),
        
        # Public-release note
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.y = element_line(color = "black", linewidth = 0.35),
        
        # Public-release note
        axis.ticks.length = grid::unit(1.5, "mm")
  ) +
  coord_cartesian(ylim = c(0, 25), clip = "off") +
  ## y-axis ticks (fixed)
  scale_y_continuous(breaks = seq(0, 25, by = 5)) +
  geom_text(
    data = ann_df,
    aes(x = continent, y = y, label = Letters, group = env_group),
    position = position_dodge(width = 0.8),
    vjust = 0, size = 2.5, inherit.aes = FALSE
  )

print(p_box_letters)

ggsave("Global_N0_box_mean_c_by_continent_env_wilcox_letters.pdf", p_box_letters,
       width = 110, height = 55, units = "mm", dpi = 300)





###############China&USA comparation##############

d_trim <- d_trim %>%
  dplyr::filter(country %in% c("China", "USA"))

# cell counts (for labels)
count_cell <- d_trim %>%
  group_by(country, env_group) %>%
  summarise(
    k     = dplyr::n(),
    n_ref = n_distinct(reference_ID),
    n_sy  = n_distinct(site_year_ID),
    .groups = "drop"
  ) %>%
  mutate(n_lab = paste0("n=", k))

print(count_cell)

# =========================================================
# 5) Main-text model choice: m0 (shared spline curve)
#    Main text focus: Env differences + country differences
# =========================================================
m0 <- nlme::lme(
  fixed  = AEN ~ env_group * country + ns(N_rate, df = 3),
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d_trim,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

print(anova(m0))

# =========================================================
# emmeans (within country) + Tukey pairwise p-values
# =========================================================
# emmeans: within each country, 3 env means
emm0 <- emmeans(m0, ~ env_group | country)

# Public-release note
emm_df <- as.data.frame(emm0) %>%
  mutate(
    pct     = emmean,
    pct_lwr = lower.CL,
    pct_upr = upper.CL,
    country = factor(country, levels = levels(d_trim$country)),
    env_group = factor(env_group, levels = env_levels)
  ) %>%
  left_join(count_cell, by = c("country","env_group")) %>%
  mutate(
    n_lab = paste0("n=", k)
  )
write.csv(emm_df,"Supplement Figure 3-3_B-1.csv")
# Public-release note
n_offset <- 0.06 * (max(emm_df$pct_upr, na.rm = TRUE) - min(emm_df$pct_lwr, na.rm = TRUE))
emm_df <- emm_df %>% mutate(y_n = pct_lwr - n_offset)

# Public-release note
pw0 <- as.data.frame(pairs(emm0, adjust = "tukey")) %>%
  tidyr::separate(contrast, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    country = factor(country, levels = levels(d_trim$country)),
    g1 = factor(g1, levels = env_levels),
    g2 = factor(g2, levels = env_levels),
    p_lab = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  )
write.csv(pw0,"Supplement Figure 3-3_B-2.csv")
# =========================
# Public-release note
# =========================
cfg <- list(
  # Public-release note
  dodge_w   = 0.72,  # Public-release note
  n_groups  = NULL,  # Public-release note
  
  # Public-release note
  n_y       = 10,  # Public-release note
  n_angle   = 45,
  n_size    = 2.6,
  n_vjust   = 1,
  
  # Public-release note
  br_base_gap = 2.0,  # Public-release note
  br_step     = 4.0,  # Public-release note
  br_tick     = 1.5,  # Public-release note
  br_lwd      = 0.25,  # Public-release note
  
  # Public-release note
  br_x_pad    = 0.00,  # Public-release note
  
  # Public-release note
  p_text_dy   = 1.5,  # Public-release note
  p_size      = 2.0,
  
  # Public-release note
  pt_size     = 1.5,
  pt_stroke   = 0.4,
  err_lwd     = 0.25,
  err_width   = 0.32,  # Public-release note
  
  # Public-release note
  base_size   = 8,
  axis_x_angle = 0,
  
  # Public-release note
  ylim_pad_top    = 3.5,
  ylim_pad_bottom = 2.0
)
# =========================
# Public-release note
# =========================
if (is.null(cfg$n_groups)) cfg$n_groups <- length(env_levels)

# Public-release note
env_offset <- tibble(
  env_group = factor(env_levels, levels = env_levels),
  idx = seq_along(env_levels)
) %>%
  mutate(
    x_off = (idx - (cfg$n_groups + 1)/2) * (cfg$dodge_w / cfg$n_groups)
  ) %>%
  dplyr::select(env_group, x_off)

# Public-release note
cont_lv <- levels(factor(emm_df$country))
cont_map <- tibble(
  country = factor(cont_lv, levels = cont_lv),
  x_base = seq_along(cont_lv)
)

# =========================
# Public-release note
# =========================
y_top <- emm_df %>%
  group_by(country) %>%
  summarise(ymax = max(pct_upr, na.rm = TRUE), .groups = "drop")

br_df <- pw0 %>%
  left_join(cont_map, by = "country") %>%
  left_join(env_offset, by = c("g1" = "env_group")) %>%
  rename(x1_off = x_off) %>%
  left_join(env_offset, by = c("g2" = "env_group")) %>%
  rename(x2_off = x_off) %>%
  left_join(y_top, by = "country") %>%
  group_by(country) %>%
  mutate(
    br_rank = row_number(),  # Public-release note
    y   = ymax + cfg$br_base_gap + (br_rank - 1) * cfg$br_step,
    x1  = x_base + x1_off - cfg$br_x_pad,
    x2  = x_base + x2_off + cfg$br_x_pad,
    xmid = (x1 + x2) / 2,
    y_txt = y + cfg$p_text_dy
  ) %>%
  ungroup()

# =========================
# Public-release note
# =========================
ymin_data <- min(emm_df$pct_lwr, na.rm = TRUE)
ymin_n    <- cfg$n_y  # Public-release note
ymin_plot <- min(ymin_data, ymin_n) - cfg$ylim_pad_bottom

ymax_plot <- 50

# =========================
# Public-release note
# =========================
p_forest_country <- ggplot(emm_df, aes(x = country, y = pct, color = env_group)) +
  
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
  
  # --- n labels ---
  geom_text(
    data = emm_df,
    aes(x = country, y = 10, label = n_lab, group = env_group, angle = cfg$n_angle),
    position = position_dodge(width = cfg$dodge_w),
    color = "black",
    size = cfg$n_size,
    vjust = cfg$n_vjust,
    inherit.aes = FALSE
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
  
  scale_color_manual(values = env_cols, breaks = env_levels) +
  labs(
    x = "country",
    y = "Agronomic efficiency of N (kg grain / kg N)",
    color = "Env"
  ) +
  theme_bw(base_size = cfg$base_size) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.70, 0.90),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.line = element_line(color = "black", linewidth = 0.35),
    # Public-release note
    axis.text.x = element_text(color = "black", hjust = 0.5, size = 8),  # Public-release note
    axis.text.y = element_text(color = "black", size = 8),
    
    # Public-release note
    axis.ticks.x = element_line(color = "black", linewidth = 0.35),
    axis.ticks.y = element_line(color = "black", linewidth = 0.35),
    
    # Public-release note
    axis.ticks.length = grid::unit(1.5, "mm")
    
  ) +
  coord_cartesian(ylim = c(ymin_plot, ymax_plot), clip = "off") +
  ## y-axis ticks (fixed)
  scale_y_continuous(breaks = seq(0, 50, by = 10))

print(p_forest_country)

ggsave(
  "China_USA_speical_N_forest_AEN.pdf",
  p_forest_country,
  width = 55, height = 55, units = "mm", dpi = 300
)

# --- within-country pairwise Wilcoxon ---
pairs_env <- list(c("Normal","PFR"), c("Normal","GFR"), c("PFR","GFR"))

wilcox_res <- d_trim %>%
  group_by(country) %>%
  group_modify(~{
    df <- .x
    out <- lapply(pairs_env, function(pp){
      a <- df$mean_t[df$env_group == pp[1]]
      b <- df$mean_t[df$env_group == pp[2]]
      if (length(a) < 2 || length(b) < 2) {
        return(data.frame(group1=pp[1], group2=pp[2], p_value=NA_real_))
      }
      data.frame(group1=pp[1], group2=pp[2],
                 p_value = wilcox.test(a, b, exact = FALSE)$p.value)
    }) %>% bind_rows()
    
    out$p_adj_bh   <- p.adjust(out$p_value, method = "BH")
    out
  }) %>%
  ungroup()

print(wilcox_res)

# --- letters (choose BH or Holm) ---
pick_p <- "p_adj_bh"

letters_df <- wilcox_res %>%
  filter(!is.na(.data[[pick_p]])) %>%
  mutate(comp = paste(group1, group2, sep = "-")) %>%
  group_by(country) %>%
  group_modify(~{
    pv <- .x[[pick_p]]
    names(pv) <- .x$comp
    L <- multcompView::multcompLetters(pv)$Letters
    data.frame(env_group = names(L), Letters = unname(L))
  }) %>%
  ungroup() %>%
  mutate(env_group = factor(env_group, levels = env_levels))

y_pos <- d_trim %>%
  group_by(country, env_group) %>%
  summarise(y = max(mean_t, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y * 1.03)

ann_df <- y_pos %>%
  left_join(letters_df, by = c("country","env_group"))

p_box_letters <- ggplot(d_trim, aes(x = country, y = mean_t, fill = env_group)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.5,
    outlier.alpha = 0.25,
    outlier.size  = 0.5,  # Public-release note
    outlier.shape = 16,  # Public-release note
    outlier.stroke = 0.2,  # Public-release note
    outlier.colour = "black"
  ) +
  stat_summary(
    aes(group = env_group),  # Public-release note
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.8),
    shape = 23, fill = NA,
    size = 1.0, stroke = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = env_cols, breaks = env_levels) +
  labs(x = "country", y = "Yield under N fertilization (Mg / ha)", fill = "Env") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.80, 0.90),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.width  = grid::unit(3, "mm"),
        legend.key.height = grid::unit(4, "mm"),
        axis.line = element_line(color = "black", linewidth = 0.35),
        # Public-release note
        axis.text.x = element_text(color = "black", hjust = 0.5, size = 8),  # Public-release note
        axis.text.y = element_text(color = "black", size = 8),
        
        # Public-release note
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.y = element_line(color = "black", linewidth = 0.35),
        
        # Public-release note
        axis.ticks.length = grid::unit(1.5, "mm")
  ) +
  coord_cartesian(ylim = c(0, 25), clip = "off") +
  ## y-axis ticks (fixed)
  scale_y_continuous(breaks = seq(0, 25, by = 5)) +
  geom_text(
    data = ann_df,
    aes(x = country, y = y, label = Letters, group = env_group),
    position = position_dodge(width = 0.8),
    vjust = 0, size = 2.5, inherit.aes = FALSE
  )

print(p_box_letters)

ggsave("China_USA_speical_N_box_mean_t_by_country_env_wilcox_letters.pdf", p_box_letters,
       width = 70, height = 55, units = "mm", dpi = 300)


# --- within-country pairwise Wilcoxon ---
pairs_env <- list(c("Normal","PFR"), c("Normal","GFR"), c("PFR","GFR"))

d_trim <- d_trim %>%
  mutate(.row_id = row_number()) %>%
  group_by(env_ID) %>%
  arrange(.row_id, .by_group = TRUE) %>%
  filter(is.na(lag(mean_c)) | mean_c != lag(mean_c)) %>%
  ungroup() %>%
  dplyr::select(-.row_id)


wilcox_res <- d_trim %>%
  group_by(country) %>%
  group_modify(~{
    df <- .x
    out <- lapply(pairs_env, function(pp){
      a <- df$mean_c[df$env_group == pp[1]]
      b <- df$mean_c[df$env_group == pp[2]]
      if (length(a) < 2 || length(b) < 2) {
        return(data.frame(group1=pp[1], group2=pp[2], p_value=NA_real_))
      }
      data.frame(group1=pp[1], group2=pp[2],
                 p_value = wilcox.test(a, b, exact = FALSE)$p.value)
    }) %>% bind_rows()
    
    out$p_adj_bh   <- p.adjust(out$p_value, method = "BH")
    out
  }) %>%
  ungroup()

print(wilcox_res)

# --- letters (choose BH or Holm) ---
pick_p <- "p_adj_bh"

letters_df <- wilcox_res %>%
  filter(!is.na(.data[[pick_p]])) %>%
  mutate(comp = paste(group1, group2, sep = "-")) %>%
  group_by(country) %>%
  group_modify(~{
    pv <- .x[[pick_p]]
    names(pv) <- .x$comp
    L <- multcompView::multcompLetters(pv)$Letters
    data.frame(env_group = names(L), Letters = unname(L))
  }) %>%
  ungroup() %>%
  mutate(env_group = factor(env_group, levels = env_levels))

y_pos <- d_trim %>%
  group_by(country, env_group) %>%
  summarise(y = max(mean_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y * 1.03)

ann_df <- y_pos %>%
  left_join(letters_df, by = c("country","env_group"))

p_box_letters <- ggplot(d_trim, aes(x = country, y = mean_c, fill = env_group)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.5,
    outlier.alpha = 0.25,
    outlier.size  = 0.5,  # Public-release note
    outlier.shape = 16,  # Public-release note
    outlier.stroke = 0.2,  # Public-release note
    outlier.colour = "black"
  ) +
  stat_summary(
    aes(group = env_group),  # Public-release note
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.8),
    shape = 23, fill = NA,
    size = 1.0, stroke = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = env_cols, breaks = env_levels) +
  labs(x = "country", y = "Yield under N0  (Mg / ha)", fill = "Env") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.80, 0.90),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.width  = grid::unit(3, "mm"),
        legend.key.height = grid::unit(4, "mm"),
        axis.line = element_line(color = "black", linewidth = 0.35),
        # Public-release note
        axis.text.x = element_text(color = "black", hjust = 0.5, size = 8),  # Public-release note
        axis.text.y = element_text(color = "black", size = 8),
        
        # Public-release note
        axis.ticks.x = element_line(color = "black", linewidth = 0.35),
        axis.ticks.y = element_line(color = "black", linewidth = 0.35),
        
        # Public-release note
        axis.ticks.length = grid::unit(1.5, "mm")
  ) +
  coord_cartesian(ylim = c(0, 25), clip = "off") +
  ## y-axis ticks (fixed)
  scale_y_continuous(breaks = seq(0, 25, by = 5)) +
  geom_text(
    data = ann_df,
    aes(x = country, y = y, label = Letters, group = env_group),
    position = position_dodge(width = 0.8),
    vjust = 0, size = 2.5, inherit.aes = FALSE
  )

print(p_box_letters)

ggsave("China_USA_N0_box_mean_t_by_country_env_wilcox_letters.pdf", p_box_letters,
       width = 55, height = 55, units = "mm", dpi = 300)
