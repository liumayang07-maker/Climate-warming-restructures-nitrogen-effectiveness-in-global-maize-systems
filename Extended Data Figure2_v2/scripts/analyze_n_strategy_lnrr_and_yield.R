# =========================================================
# 0) Packages & data
# =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(nlme)
  library(emmeans)
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
col_lnRR       <- "lnRR"
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
col_N_rate_cluster <- "N_rate_cluster"
col_env_ID     <- "env_ID"
col_N0         <- "Label_of_control_treatment"
col_N_strategy <- "N_strategy"
col_N_S        <- "N_S"
col_N_V        <- "N_V"
col_N_R        <- "N_R"
col_plant_density <- "plant_density_cluster"
col_region_N_label <- "region_N_label"

env_levels <- c("Normal","PFR","GFR")

# ---- trim settings ----
q_lo <- 0.05
q_hi <- 0.95

summarise_yield_box <- function(d, panel, subset_label) {
  d %>%
    filter(is.finite(mean_t)) %>%
    mutate(
      env_group = factor(env_group, levels = env_levels),
      N_strategy = droplevels(factor(N_strategy))
    ) %>%
    group_by(env_group, N_strategy) %>%
    summarise(
      panel = panel,
      subset = subset_label,
      n = dplyr::n(),
      n_ref = n_distinct(reference_ID),
      n_site_year = n_distinct(site_year_ID),
      mean = mean(mean_t, na.rm = TRUE),
      sd = sd(mean_t, na.rm = TRUE),
      se = sd / sqrt(n),
      median = median(mean_t, na.rm = TRUE),
      q25 = quantile(mean_t, 0.25, na.rm = TRUE),
      q75 = quantile(mean_t, 0.75, na.rm = TRUE),
      iqr = IQR(mean_t, na.rm = TRUE),
      min = min(mean_t, na.rm = TRUE),
      max = max(mean_t, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    relocate(panel, subset, env_group, N_strategy)
}

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
    N_rate_cluster = as.factor(.data[[col_N_rate_cluster]]),
    EDH_label    = as.factor(.data[[col_EDH_label]]),
    region_N_label = as.factor(.data[[col_region_N_label]]),
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
    Nc = N_rate - mean(N_rate, na.rm = TRUE)  # global centering
  )

# ---- subset: RN + selected countries + four N-management strategies ----
d_trim <- d0 %>%
  filter(region_N_label %in% c("RN")) %>%
  filter(country %in% c("USA", "China")) %>%  
  filter(!is.na(env_group), !is.na(N_strategy), is.finite(lnRR)) %>%
  mutate(
    env_group  = droplevels(env_group),
    N_strategy = droplevels(factor(N_strategy))
  ) %>%
  filter(N_strategy %in% c("S","SV","SR","SVR"))

write.csv(d_trim, "USA_China_N_strategy_data_set.csv")

# ---- within-env trimming of lnRR (5%-95%) ----
dA <- d_trim %>%
  group_by(env_group) %>%
  mutate(
    q05  = quantile(lnRR, q_lo, na.rm = TRUE),
    q95  = quantile(lnRR, q_hi, na.rm = TRUE),
    keep = lnRR >= q05 & lnRR <= q95
  ) %>%
  ungroup() %>%
  filter(keep)

# =========================================================
# 2) Cell counts (for labels)
# =========================================================
count_cell <- dA %>%
  group_by(N_strategy, env_group) %>%
  summarise(
    k     = dplyr::n(),
    n_ref = n_distinct(reference_ID),
    n_sy  = n_distinct(site_year_ID),
    .groups = "drop"
  ) %>%
  mutate(n_lab = paste0("n=", k))

print(count_cell)

# =========================================================
# 3) Mixed model
# =========================================================
m0 <- nlme::lme(
  fixed  = lnRR ~ env_group * N_strategy + Nc,
  random = ~ 1 | reference_ID/site_year_ID,
  data   = dA,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

print(anova(m0))

# =========================================================
# 4) emmeans + omnibus test within each env_group
# =========================================================
env_levels2 <- levels(factor(dA$env_group))
N_levels2   <- c("S","SV","SR","SVR")
emm1 <- emmeans(m0, ~ N_strategy | env_group)

emm_df <- as.data.frame(emm1) %>%
  mutate(
    env_group  = factor(env_group, levels = env_levels2),
    N_strategy = factor(N_strategy, levels = N_levels2),
    lnRR     = emmean,
    lnRR_lwr = lower.CL,
    lnRR_upr = upper.CL,
    pct     = (exp(lnRR)     - 1) * 100,
    pct_lwr = (exp(lnRR_lwr) - 1) * 100,
    pct_upr = (exp(lnRR_upr) - 1) * 100
  ) %>%
  left_join(count_cell, by = c("env_group","N_strategy")) %>%
  mutate(n_lab = paste0("n=", k))

write.csv(emm_df, "Extended Data Figure 2a_data1.csv")

# ---- Omnibus per env_group: H0 (within env): all N_strategy means equal ----
jt0 <- as.data.frame(test(pairs(emm1), joint = TRUE))

# robustly standardize the grouping column to "env_group"
grp_col <- setdiff(names(jt0), c("df1","df2","F.ratio","p.value","statistic","p.value.adj","p.adj"))
if (length(grp_col) == 0) grp_col <- names(jt0)[1]   # fallback
jt_by_env <- jt0 %>%
  rename(env_group = all_of(grp_col[1])) %>%
  mutate(
    env_group = factor(env_group, levels = env_levels2),
    p_lab = case_when(
      p.value < 0.001 ~ "P<0.001",
      TRUE ~ paste0("P=", sprintf("%.3f", p.value))
    )
  )

print(jt_by_env)

write.csv(jt_by_env, "Extended Data Figure 2a_data2.csv")

# ---- position p labels above the highest CI in each env ----
y_top <- emm_df %>%
  group_by(env_group) %>%
  summarise(ymax = max(pct_upr, na.rm = TRUE), .groups = "drop")

p_df <- jt_by_env %>%
  left_join(y_top, by = "env_group") %>%
  mutate(y = ymax + 8)   # Adjust vertical spacing if labels overlap

# =========================================================
# 5) Plot (means+CI + n labels + ONE omnibus p per env)
# =========================================================
cfg <- list(
  dodge_w   = 0.62,
  n_size    = 1.9,
  n_vjust   = 0.9,
  p_size    = 2.0,
  pt_size   = 1.35,
  pt_stroke = 0.35,
  err_lwd   = 0.25,
  err_width = 0.28,
  base_size = 5.0
)

N_cols <- setNames(c("#83A0A4", "#B2A1C0", "#F2C666", "#989568"), N_levels2)
ymin_plot <- -5
ymax_plot <- 120

p_forest_env_strategy_pct <- ggplot(emm_df, aes(x = env_group, y = pct, color = N_strategy)) +
  
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
  
  # n labels (aligned with dodge)
  geom_text(
    data = emm_df,
    aes(x = env_group, y = 10, label = n_lab, group = N_strategy),
    position = position_dodge(width = cfg$dodge_w),
    color = "black",
    size = cfg$n_size,
    vjust = cfg$n_vjust,
    inherit.aes = FALSE,
    angle = 45
  ) +
  
  # one omnibus p per env_group
  geom_text(
    data = p_df,
    aes(x = env_group, y = y, label = p_lab),
    inherit.aes = FALSE,
    color = "black",
    size = cfg$p_size,
    vjust = 0
  ) +
  
  scale_color_manual(values = N_cols, breaks = N_levels2) +
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
    
    plot.margin = margin(2, 6, 2, 4)
  ) +
  coord_cartesian(ylim = c(ymin_plot, ymax_plot), clip = "off") +
  scale_y_continuous(breaks = seq(0, 120, by = 30))

print(p_forest_env_strategy_pct)

ggsave(
  "USA_China_lnRR_forest_OmnibusP_N_strategy.pdf",
  p_forest_env_strategy_pct,
  width = 76, height = 38, units = "mm", dpi = 300
)


#####################################################################
# -------------------------
# 1) within-env_group pairwise Wilcoxon (density comparisons)
# -------------------------
pairs_N <- list(
  c(N_levels2[1], N_levels2[2]),
  c(N_levels2[1], N_levels2[3]),
  c(N_levels2[1], N_levels2[4]),
  c(N_levels2[2], N_levels2[3]),
  c(N_levels2[2], N_levels2[4]),
  c(N_levels2[3], N_levels2[4])
)

wilcox_res2 <- dA %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  ) %>%
  group_by(env_group) %>%
  group_modify(~{
    df <- .x
    out <- lapply(pairs_N, function(pp){
      a <- df$mean_t[df$N_strategy == pp[1]]
      b <- df$mean_t[df$N_strategy == pp[2]]
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
write.csv(wilcox_res2, "Extended Data Figure 2g_data.csv", row.names = FALSE)
yield_summary2 <- summarise_yield_box(dA, "Extended Data Figure 2g", "USA_China")
write.csv(yield_summary2, "Extended Data Figure 2g_yield_summary.csv", row.names = FALSE)

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
    data.frame(N_strategy = names(L), Letters = unname(L))
  }) %>%
  ungroup() %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  )

# -------------------------
# 3) annotation y positions
# -------------------------
y_pos2 <- dA %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  ) %>%
  group_by(env_group, N_strategy) %>%
  summarise(y = max(mean_t, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y * 1.03)

ann_df2 <- y_pos2 %>%
  left_join(letters_df2, by = c("env_group","N_strategy"))

# -------------------------
# 4) plot: X = env_group, fill = N_strategy
# -------------------------
# Optional alternative palette:
# Optional alternative palette:

p_box_letters2 <- ggplot(
  dA %>% mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  ),
  aes(x = env_group, y = mean_t, fill = N_strategy)
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
    aes(group = N_strategy),
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.8),
    shape = 23, fill = NA,
    size = 1.0, stroke = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = N_cols, breaks = N_levels2) +
  labs(
    x = "Env group",
    y = "Yield under N fertilization (Mg / ha)",
    fill = "Nitrogen strategy"
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
    aes(x = env_group, y = y, label = Letters, group = N_strategy),
    position = position_dodge(width = 0.8),
    vjust = 0, size = 2.5, inherit.aes = FALSE
  )

print(p_box_letters2)

ggsave("USA_China_box_yield_N_strategy_wilcox_letters.pdf", p_box_letters2,
       width = 153, height = 45, units = "mm", dpi = 300)

# China subset


# ---- subset: RN + selected countries + four N-management strategies ----
d_trim <- d0 %>%
  filter(region_N_label %in% c("RN")) %>%
  filter(country %in% c("China")) %>%  
  filter(!is.na(env_group), !is.na(N_strategy), is.finite(lnRR)) %>%
  mutate(
    env_group  = droplevels(env_group),
    N_strategy = droplevels(factor(N_strategy))
  ) %>%
  filter(N_strategy %in% c("S","SV","SR","SVR"))

# ---- within-env trimming of lnRR (5%-95%) ----
dA <- d_trim %>%
  group_by(env_group) %>%
  mutate(
    q05  = quantile(lnRR, q_lo, na.rm = TRUE),
    q95  = quantile(lnRR, q_hi, na.rm = TRUE),
    keep = lnRR >= q05 & lnRR <= q95
  ) %>%
  ungroup() %>%
  filter(keep)

# =========================================================
# 2) Cell counts (for labels)
# =========================================================
count_cell <- dA %>%
  group_by(N_strategy, env_group) %>%
  summarise(
    k     = dplyr::n(),
    n_ref = n_distinct(reference_ID),
    n_sy  = n_distinct(site_year_ID),
    .groups = "drop"
  ) %>%
  mutate(n_lab = paste0("n=", k))

print(count_cell)

# =========================================================
# 3) Mixed model
# =========================================================
m0 <- nlme::lme(
  fixed  = lnRR ~ env_group * N_strategy + Nc,
  random = ~ 1 | reference_ID/site_year_ID,
  data   = dA,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

print(anova(m0))

# =========================================================
# 4) emmeans + omnibus test within each env_group
# =========================================================
env_levels2 <- levels(factor(dA$env_group))
N_levels2   <- c("S","SV","SR","SVR")
emm1 <- emmeans(m0, ~ N_strategy | env_group)

emm_df <- as.data.frame(emm1) %>%
  mutate(
    env_group  = factor(env_group, levels = env_levels2),
    N_strategy = factor(N_strategy, levels = N_levels2),
    lnRR     = emmean,
    lnRR_lwr = lower.CL,
    lnRR_upr = upper.CL,
    pct     = (exp(lnRR)     - 1) * 100,
    pct_lwr = (exp(lnRR_lwr) - 1) * 100,
    pct_upr = (exp(lnRR_upr) - 1) * 100
  ) %>%
  left_join(count_cell, by = c("env_group","N_strategy")) %>%
  mutate(n_lab = paste0("n=", k))
write.csv(emm_df, "Extended Data Figure 2c_data1.csv")

# Optional grouped contrast: (SR + SVR) - (S + SV)
# Coefficient order must match levels(N_strategy) = c("S", "SV", "SR", "SVR")
L <- list("SR+SVR vs S+SV" = c(-0.5, -0.5, 0.5, 0.5))

ct <- contrast(emm1, method = L) %>%
  test(adjust = "none") %>%     # One contrast per environment; no further adjustment applied here
  as.data.frame()

# Optional: apply BH correction across the three environment-specific contrasts
ct <- ct %>%
  mutate(p_BH = p.adjust(p.value, method = "BH"))

ct

# ---- Omnibus per env_group: H0 (within env): all N_strategy means equal ----
jt0 <- as.data.frame(test(pairs(emm1), joint = TRUE))

# robustly standardize the grouping column to "env_group"
grp_col <- setdiff(names(jt0), c("df1","df2","F.ratio","p.value","statistic","p.value.adj","p.adj"))
if (length(grp_col) == 0) grp_col <- names(jt0)[1]   # fallback
jt_by_env <- jt0 %>%
  rename(env_group = all_of(grp_col[1])) %>%
  mutate(
    env_group = factor(env_group, levels = env_levels2),
    p_lab = case_when(
      p.value < 0.001 ~ "P<0.001",
      TRUE ~ paste0("P=", sprintf("%.3f", p.value))
    )
  )

print(jt_by_env)
write.csv(jt_by_env, "Extended Data Figure 2c_data2.csv")
# ---- position p labels above the highest CI in each env ----
y_top <- emm_df %>%
  group_by(env_group) %>%
  summarise(ymax = max(pct_upr, na.rm = TRUE), .groups = "drop")

p_df <- jt_by_env %>%
  left_join(y_top, by = "env_group") %>%
  mutate(y = ymax + 8)   # Adjust vertical spacing if labels overlap

# =========================================================
# 5) Plot (means+CI + n labels + ONE omnibus p per env)
# =========================================================
cfg <- list(
  dodge_w   = 0.62,
  n_size    = 1.9,
  n_vjust   = 0.9,
  p_size    = 2.0,
  pt_size   = 1.35,
  pt_stroke = 0.35,
  err_lwd   = 0.25,
  err_width = 0.28,
  base_size = 5.0
)

N_cols <- setNames(c("#83A0A4", "#B2A1C0", "#F2C666", "#989568"), N_levels2)
ymin_plot <- -5
ymax_plot <- 120

p_forest_env_strategy_pct <- ggplot(emm_df, aes(x = env_group, y = pct, color = N_strategy)) +
  
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
  
  # n labels (aligned with dodge)
  geom_text(
    data = emm_df,
    aes(x = env_group, y = 10, label = n_lab, group = N_strategy),
    position = position_dodge(width = cfg$dodge_w),
    color = "black",
    size = cfg$n_size,
    vjust = cfg$n_vjust,
    inherit.aes = FALSE,
    angle = 45
  ) +
  
  # one omnibus p per env_group
  geom_text(
    data = p_df,
    aes(x = env_group, y = y, label = p_lab),
    inherit.aes = FALSE,
    color = "black",
    size = cfg$p_size,
    vjust = 0
  ) +
  
  scale_color_manual(values = N_cols, breaks = N_levels2) +
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
    
    plot.margin = margin(2, 6, 2, 4)
  ) +
  coord_cartesian(ylim = c(ymin_plot, ymax_plot), clip = "off") +
  scale_y_continuous(breaks = seq(0, 120, by = 30))

print(p_forest_env_strategy_pct)

ggsave(
  "China_lnRR_forest_OmnibusP_N_strategy.pdf",
  p_forest_env_strategy_pct,
  width = 76, height = 38, units = "mm", dpi = 300
)
####################################################################################
# -------------------------
# 1) within-env_group pairwise Wilcoxon (density comparisons)
# -------------------------
pairs_N <- list(
  c(N_levels2[1], N_levels2[2]),
  c(N_levels2[1], N_levels2[3]),
  c(N_levels2[1], N_levels2[4]),
  c(N_levels2[2], N_levels2[3]),
  c(N_levels2[2], N_levels2[4]),
  c(N_levels2[3], N_levels2[4])
)

wilcox_res2 <- dA %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  ) %>%
  group_by(env_group) %>%
  group_modify(~{
    df <- .x
    out <- lapply(pairs_N, function(pp){
      a <- df$mean_t[df$N_strategy == pp[1]]
      b <- df$mean_t[df$N_strategy == pp[2]]
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
write.csv(wilcox_res2, "Extended Data Figure 2h_china_yield_wilcox.csv", row.names = FALSE)
yield_summary2 <- summarise_yield_box(dA, "Extended Data Figure 2h", "China")
write.csv(yield_summary2, "Extended Data Figure 2h_china_yield_summary.csv", row.names = FALSE)

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
    data.frame(N_strategy = names(L), Letters = unname(L))
  }) %>%
  ungroup() %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  )

# -------------------------
# 3) annotation y positions
# -------------------------
y_pos2 <- dA %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  ) %>%
  group_by(env_group, N_strategy) %>%
  summarise(y = max(mean_t, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y * 1.03)

ann_df2 <- y_pos2 %>%
  left_join(letters_df2, by = c("env_group","N_strategy"))

# -------------------------
# 4) plot: X = env_group, fill = N_strategy
# -------------------------
# Optional alternative palette:
# Optional alternative palette:

p_box_letters2 <- ggplot(
  dA %>% mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  ),
  aes(x = env_group, y = mean_t, fill = N_strategy)
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
    aes(group = N_strategy),
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.8),
    shape = 23, fill = NA,
    size = 1.0, stroke = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = N_cols, breaks = N_levels2) +
  labs(
    x = "Env group",
    y = "Yield under N fertilization (Mg / ha)",
    fill = "Nitrogen strategy"
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
    aes(x = env_group, y = y, label = Letters, group = N_strategy),
    position = position_dodge(width = 0.8),
    vjust = 0, size = 2.5, inherit.aes = FALSE
  )

print(p_box_letters2)

ggsave("China_box_yield_N_strategy_wilcox_letters.pdf", p_box_letters2,
       width = 100, height = 45, units = "mm", dpi = 300)

# USA subset
# ---- subset: RN + selected countries + four N-management strategies ----
d_trim <- d0 %>%
  filter(region_N_label %in% c("RN"),
         country %in% c("USA")) %>%
  filter(!is.na(env_group), !is.na(N_strategy), is.finite(lnRR)) %>%
  mutate(
    env_group  = droplevels(env_group),
    N_strategy = droplevels(factor(N_strategy))
  ) %>%
  filter(N_strategy %in% c("S","SV","SR","SVR"))

d_trim <- d_trim %>%
  mutate(
    N_strategy = case_when(
      N_strategy == "S" ~ "S",
      N_strategy %in% c("SV","SR","SVR") ~ "SS",
      TRUE ~ NA_character_
    ),
    N_strategy = factor(N_strategy, levels = c("S","SS"))
  ) %>%
  filter(!is.na(N_strategy))

# ---- within-env trimming of lnRR (5%-95%) ----
dA <- d_trim %>%
  group_by(env_group) %>%
  mutate(
    q05  = quantile(lnRR, q_lo, na.rm = TRUE),
    q95  = quantile(lnRR, q_hi, na.rm = TRUE),
    keep = lnRR >= q05 & lnRR <= q95
  ) %>%
  ungroup() %>%
  filter(keep)

# =========================================================
# 2) Cell counts (for labels)
# =========================================================
count_cell <- dA %>%
  group_by(N_strategy, env_group) %>%
  summarise(
    k     = dplyr::n(),
    n_ref = n_distinct(reference_ID),
    n_sy  = n_distinct(site_year_ID),
    .groups = "drop"
  ) %>%
  mutate(n_lab = paste0("n=", k))

print(count_cell)

# =========================================================
# 3) Mixed model
# =========================================================
m0 <- nlme::lme(
  fixed  = lnRR ~ env_group * N_strategy + Nc,
  random = ~ 1 | reference_ID/site_year_ID,
  data   = dA,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

print(anova(m0))

# =========================================================
# 4) emmeans + omnibus test within each env_group
# =========================================================
env_levels2 <- levels(factor(dA$env_group))
N_levels2   <- c("S","SS")
emm1 <- emmeans(m0, ~ N_strategy | env_group)

emm_df <- as.data.frame(emm1) %>%
  mutate(
    env_group  = factor(env_group, levels = env_levels2),
    N_strategy = factor(N_strategy, levels = N_levels2),
    lnRR     = emmean,
    lnRR_lwr = lower.CL,
    lnRR_upr = upper.CL,
    pct     = (exp(lnRR)     - 1) * 100,
    pct_lwr = (exp(lnRR_lwr) - 1) * 100,
    pct_upr = (exp(lnRR_upr) - 1) * 100
  ) %>%
  left_join(count_cell, by = c("env_group","N_strategy")) %>%
  mutate(n_lab = paste0("n=", k))
write.csv(emm_df, "Extended Data Figure 2e_data1.csv")
# ---- Omnibus per env_group: H0 (within env): all N_strategy means equal ----
jt0 <- as.data.frame(test(pairs(emm1), joint = TRUE))

# robustly standardize the grouping column to "env_group"
grp_col <- setdiff(names(jt0), c("df1","df2","F.ratio","p.value","statistic","p.value.adj","p.adj"))
if (length(grp_col) == 0) grp_col <- names(jt0)[1]   # fallback
jt_by_env <- jt0 %>%
  rename(env_group = all_of(grp_col[1])) %>%
  mutate(
    env_group = factor(env_group, levels = env_levels2),
    p_lab = case_when(
      p.value < 0.001 ~ "P<0.001",
      TRUE ~ paste0("P=", sprintf("%.3f", p.value))
    )
  )

print(jt_by_env)
write.csv(jt_by_env, "Extended Data Figure 2e_data2.csv")
# ---- position p labels above the highest CI in each env ----
y_top <- emm_df %>%
  group_by(env_group) %>%
  summarise(ymax = max(pct_upr, na.rm = TRUE), .groups = "drop")

p_df <- jt_by_env %>%
  left_join(y_top, by = "env_group") %>%
  mutate(y = ymax + 8)   # Adjust vertical spacing if labels overlap

# =========================================================
# 5) Plot (means+CI + n labels + ONE omnibus p per env)
# =========================================================
cfg <- list(
  dodge_w   = 0.62,
  n_size    = 1.9,
  n_vjust   = 0.9,
  p_size    = 2.0,
  pt_size   = 1.35,
  pt_stroke = 0.35,
  err_lwd   = 0.25,
  err_width = 0.28,
  base_size = 5.0
)

N_cols <- setNames(c("#83A0A4", "#F2C666"), N_levels2)
ymin_plot <- 20
ymax_plot <- 180

p_forest_env_strategy_pct <- ggplot(emm_df, aes(x = env_group, y = pct, color = N_strategy)) +
  
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
  
  # n labels (aligned with dodge)
  geom_text(
    data = emm_df,
    aes(x = env_group, y = 35, label = n_lab, group = N_strategy),
    position = position_dodge(width = cfg$dodge_w),
    color = "black",
    size = cfg$n_size,
    vjust = cfg$n_vjust,
    inherit.aes = FALSE,
    angle = 45
  ) +
  
  # one omnibus p per env_group
  geom_text(
    data = p_df,
    aes(x = env_group, y = y, label = p_lab),
    inherit.aes = FALSE,
    color = "black",
    size = cfg$p_size,
    vjust = 0
  ) +
  
  scale_color_manual(values = N_cols, breaks = N_levels2) +
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
    
    plot.margin = margin(2, 6, 2, 4)
  ) +
  coord_cartesian(ylim = c(ymin_plot, ymax_plot), clip = "off") +
  scale_y_continuous(breaks = seq(20, 180, by = 40))

print(p_forest_env_strategy_pct)

ggsave(
  "USA_lnRR_forest_OmnibusP_N_strategy.pdf",
  p_forest_env_strategy_pct,
  width = 76, height = 38, units = "mm", dpi = 300
)

# -------------------------
# 1) within-env_group pairwise Wilcoxon (density comparisons)
# -------------------------

N_levels2   <- c("S","SS")
pairs_N <- list(
  c(N_levels2[1], N_levels2[2])
)

wilcox_res2 <- dA %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  ) %>%
  group_by(env_group) %>%
  group_modify(~{
    df <- .x
    out <- lapply(pairs_N, function(pp){
      a <- df$mean_t[df$N_strategy == pp[1]]
      b <- df$mean_t[df$N_strategy == pp[2]]
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
write.csv(wilcox_res2, "Extended Data Figure 2i_USA_yield_wilcox.csv", row.names = FALSE)
yield_summary2 <- summarise_yield_box(dA, "Extended Data Figure 2i", "USA")
write.csv(yield_summary2, "Extended Data Figure 2i_USA_yield_summary.csv", row.names = FALSE)

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
    data.frame(N_strategy = names(L), Letters = unname(L))
  }) %>%
  ungroup() %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  )

# -------------------------
# 3) annotation y positions
# -------------------------
y_pos2 <- dA %>%
  mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  ) %>%
  group_by(env_group, N_strategy) %>%
  summarise(y = max(mean_t, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y * 1.03)

ann_df2 <- y_pos2 %>%
  left_join(letters_df2, by = c("env_group","N_strategy"))

# -------------------------
# 4) plot: X = env_group, fill = N_strategy
# -------------------------
# Optional alternative palette:
# Optional alternative palette:

p_box_letters2 <- ggplot(
  dA %>% mutate(
    env_group     = factor(env_group, levels = env_levels),
    N_strategy = factor(N_strategy, levels = N_levels2)
  ),
  aes(x = env_group, y = mean_t, fill = N_strategy)
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
    aes(group = N_strategy),
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.8),
    shape = 23, fill = NA,
    size = 1.0, stroke = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = N_cols, breaks = N_levels2) +
  labs(
    x = "Env group",
    y = "Yield under N fertilization (Mg / ha)",
    fill = "Nitrogen strategy"
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
    aes(x = env_group, y = y, label = Letters, group = N_strategy),
    position = position_dodge(width = 0.8),
    vjust = 0, size = 2.5, inherit.aes = FALSE
  )

print(p_box_letters2)

ggsave("USA_box_yield_N_strategy_wilcox_letters.pdf", p_box_letters2,
       width = 53, height = 45, units = "mm", dpi = 300)
