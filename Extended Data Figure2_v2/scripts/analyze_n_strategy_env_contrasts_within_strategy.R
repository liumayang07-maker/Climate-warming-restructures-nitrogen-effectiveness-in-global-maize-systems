# Supplementary contrasts for Extended Data Figure 2.
# Purpose: compare heat-exposure environments within each N-management strategy.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(nlme)
  library(emmeans)
  library(multcompView)
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

env_levels <- c("Normal", "PFR", "GFR")
q_lo <- 0.05
q_hi <- 0.95

dat <- read.csv(file.path(input_dir, "meta_data_v2.csv"), stringsAsFactors = FALSE)

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

prep_metric <- function(metric_col, subset_countries, collapse_usa_strategy = FALSE) {
  dat %>%
    mutate(
      response = suppressWarnings(as.numeric(.data[[metric_col]])),
      env_group = factor(.data[["Env_Type"]], levels = env_levels),
      reference_ID = as.factor(.data[["reference_ID"]]),
      site_year_ID = as.factor(.data[["site_year_ID"]]),
      continent = as.factor(.data[["continent"]]),
      country = as.factor(.data[["country"]]),
      N_strategy = as.factor(.data[["N_strategy"]]),
      region_N_label = as.factor(.data[["region_N_label"]]),
      N_rate = suppressWarnings(as.numeric(.data[["N_rate"]])),
      mean_t = suppressWarnings(as.numeric(.data[["mean_t"]]))
    ) %>%
    filter(
      is.finite(response),
      !is.na(env_group),
      !is.na(reference_ID),
      !is.na(site_year_ID),
      !is.na(continent),
      is.finite(N_rate)
    ) %>%
    filter(region_N_label %in% "RN", country %in% subset_countries) %>%
    filter(!is.na(N_strategy), N_strategy %in% c("S", "SV", "SR", "SVR")) %>%
    mutate(
      N_strategy = if (collapse_usa_strategy) {
        ifelse(as.character(N_strategy) == "S", "S", "SS")
      } else {
        as.character(N_strategy)
      },
      N_strategy = factor(N_strategy, levels = if (collapse_usa_strategy) c("S", "SS") else c("S", "SV", "SR", "SVR")),
      Nc = N_rate - mean(N_rate, na.rm = TRUE)
    ) %>%
    group_by(env_group) %>%
    mutate(
      q05 = quantile(response, q_lo, na.rm = TRUE),
      q95 = quantile(response, q_hi, na.rm = TRUE),
      keep = response >= q05 & response <= q95
    ) %>%
    ungroup() %>%
    filter(keep) %>%
    mutate(
      env_group = droplevels(env_group),
      N_strategy = droplevels(N_strategy),
      reference_ID = droplevels(reference_ID),
      site_year_ID = droplevels(site_year_ID)
    )
}

model_env_within_strategy <- function(d, metric_name, response_transform = c("identity", "lnrr_pct")) {
  response_transform <- match.arg(response_transform)

  model <- nlme::lme(
    fixed = response ~ env_group * N_strategy + Nc,
    random = ~ 1 | reference_ID/site_year_ID,
    data = d,
    method = "REML",
    na.action = na.omit,
    control = lmeControl(opt = "optim")
  )

  emm <- emmeans(model, ~ env_group | N_strategy)

  estimates <- as.data.frame(emm) %>%
    rename(estimate_model = emmean) %>%
    mutate(
      metric = metric_name,
      response_estimate = if (response_transform == "lnrr_pct") (exp(estimate_model) - 1) * 100 else estimate_model,
      response_lwr = if (response_transform == "lnrr_pct") (exp(lower.CL) - 1) * 100 else lower.CL,
      response_upr = if (response_transform == "lnrr_pct") (exp(upper.CL) - 1) * 100 else upper.CL
    )

  count_cell <- d %>%
    group_by(N_strategy, env_group) %>%
    summarise(
      k = dplyr::n(),
      n_ref = n_distinct(reference_ID),
      n_sy = n_distinct(site_year_ID),
      .groups = "drop"
    )

  estimates <- estimates %>%
    left_join(count_cell, by = c("N_strategy", "env_group"))

  contrasts <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    tidyr::separate(contrast, into = c("g1", "g2"), sep = " - ", remove = FALSE) %>%
    mutate(
      metric = metric_name,
      p_lab = format_p(p.value)
    ) %>%
    left_join(
      estimates %>%
        select(N_strategy, env_group, response_estimate, response_lwr, response_upr, k) %>%
        rename(
          g1 = env_group,
          g1_estimate = response_estimate,
          g1_lwr = response_lwr,
          g1_upr = response_upr,
          g1_k = k
        ),
      by = c("N_strategy", "g1")
    ) %>%
    left_join(
      estimates %>%
        select(N_strategy, env_group, response_estimate, response_lwr, response_upr, k) %>%
        rename(
          g2 = env_group,
          g2_estimate = response_estimate,
          g2_lwr = response_lwr,
          g2_upr = response_upr,
          g2_k = k
        ),
      by = c("N_strategy", "g2")
    ) %>%
    mutate(response_difference = g1_estimate - g2_estimate)

  jt0 <- as.data.frame(test(pairs(emm), joint = TRUE))
  grp_col <- setdiff(names(jt0), c("df1", "df2", "F.ratio", "p.value", "statistic", "p.value.adj", "p.adj"))
  if (length(grp_col) == 0) grp_col <- names(jt0)[1]
  omnibus <- jt0 %>%
    rename(N_strategy = all_of(grp_col[1])) %>%
    mutate(
      N_strategy = factor(N_strategy, levels = levels(droplevels(d$N_strategy))),
      metric = metric_name,
      p_lab = case_when(
        p.value < 0.001 ~ "P<0.001",
        TRUE ~ paste0("P=", sprintf("%.3f", p.value))
      )
    )

  list(model = model, estimates = estimates, contrasts = contrasts, omnibus = omnibus)
}

yield_wilcox_within_strategy <- function(d) {
  pairs_env <- list(c("Normal", "PFR"), c("Normal", "GFR"), c("PFR", "GFR"))

  d_yield <- d %>%
    filter(is.finite(mean_t)) %>%
    mutate(
      env_group = factor(env_group, levels = env_levels),
      N_strategy = droplevels(N_strategy)
    )

  summary <- d_yield %>%
    group_by(N_strategy, env_group) %>%
    summarise(
      k = dplyr::n(),
      mean = mean(mean_t, na.rm = TRUE),
      median = median(mean_t, na.rm = TRUE),
      q25 = quantile(mean_t, 0.25, na.rm = TRUE),
      q75 = quantile(mean_t, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  contrasts <- d_yield %>%
    group_by(N_strategy) %>%
    group_modify(~{
      df <- .x
      out <- lapply(pairs_env, function(pp) {
        a <- df$mean_t[df$env_group == pp[1]]
        b <- df$mean_t[df$env_group == pp[2]]
        if (length(a) < 2 || length(b) < 2) {
          return(data.frame(g1 = pp[1], g2 = pp[2], p_value = NA_real_))
        }
        data.frame(
          g1 = pp[1],
          g2 = pp[2],
          p_value = wilcox.test(a, b, exact = FALSE)$p.value
        )
      }) %>% bind_rows()
      out$p_adj_bh <- p.adjust(out$p_value, method = "BH")
      out$p_adj_holm <- p.adjust(out$p_value, method = "holm")
      out
    }) %>%
    ungroup() %>%
    mutate(
      metric = "yield_under_N",
      p_lab_bh = format_p(p_adj_bh)
    ) %>%
    left_join(
      summary %>%
        select(N_strategy, env_group, mean, median, k) %>%
        rename(g1 = env_group, g1_mean = mean, g1_median = median, g1_k = k),
      by = c("N_strategy", "g1")
    ) %>%
    left_join(
      summary %>%
        select(N_strategy, env_group, mean, median, k) %>%
        rename(g2 = env_group, g2_mean = mean, g2_median = median, g2_k = k),
      by = c("N_strategy", "g2")
    ) %>%
    mutate(mean_difference = g1_mean - g2_mean)

  list(summary = summary, contrasts = contrasts)
}

plot_model_forest <- function(estimates, contrasts, y_lab, out_name, y_breaks, y_limits) {
  strategy_levels <- levels(factor(estimates$N_strategy))
  estimates <- estimates %>%
    mutate(
      N_strategy = factor(N_strategy, levels = strategy_levels),
      env_group = factor(env_group, levels = env_levels),
      n_lab = paste0("n=", k)
    )

  cfg <- list(
    dodge_w = 0.62,
    n_groups = length(env_levels),
    n_y = y_limits[1] + 0.05 * diff(y_limits),
    br_base_gap = 0.08 * diff(y_limits),
    br_step = 0.12 * diff(y_limits),
    br_tick = 0.035 * diff(y_limits),
    p_text_dy = 0.035 * diff(y_limits),
    br_lwd = 0.25,
    p_size = 1.9,
    n_size = 1.9,
    pt_size = 1.35,
    pt_stroke = 0.35,
    err_lwd = 0.25,
    err_width = 0.28,
    base_size = 5.0
  )

  env_offset <- tibble(
    env_group = factor(env_levels, levels = env_levels),
    idx = seq_along(env_levels)
  ) %>%
    mutate(x_off = (idx - (cfg$n_groups + 1) / 2) * (cfg$dodge_w / cfg$n_groups)) %>%
    select(env_group, x_off)

  x_map <- tibble(
    N_strategy = factor(strategy_levels, levels = strategy_levels),
    x_base = seq_along(strategy_levels)
  )

  y_top <- estimates %>%
    group_by(N_strategy) %>%
    summarise(ymax = max(response_upr, na.rm = TRUE), .groups = "drop")

  br_df <- contrasts %>%
    mutate(
      N_strategy = factor(N_strategy, levels = strategy_levels),
      g1 = factor(g1, levels = env_levels),
      g2 = factor(g2, levels = env_levels)
    ) %>%
    left_join(x_map, by = "N_strategy") %>%
    left_join(env_offset, by = c("g1" = "env_group")) %>%
    rename(x1_off = x_off) %>%
    left_join(env_offset, by = c("g2" = "env_group")) %>%
    rename(x2_off = x_off) %>%
    left_join(y_top, by = "N_strategy") %>%
    group_by(N_strategy) %>%
    mutate(
      br_rank = row_number(),
      y = ymax + cfg$br_base_gap + (br_rank - 1) * cfg$br_step,
      x1 = x_base + x1_off,
      x2 = x_base + x2_off,
      xmid = (x1 + x2) / 2,
      y_txt = y + cfg$p_text_dy
    ) %>%
    ungroup()

  p <- ggplot(estimates, aes(x = N_strategy, y = response_estimate, color = env_group)) +
    geom_errorbar(
      aes(ymin = response_lwr, ymax = response_upr),
      position = position_dodge(width = cfg$dodge_w),
      width = cfg$err_width,
      linewidth = cfg$err_lwd
    ) +
    geom_point(
      position = position_dodge(width = cfg$dodge_w),
      shape = 21,
      fill = "white",
      size = cfg$pt_size,
      stroke = cfg$pt_stroke
    ) +
    geom_text(
      aes(x = N_strategy, y = cfg$n_y, label = n_lab, group = env_group),
      position = position_dodge(width = cfg$dodge_w),
      color = "black",
      size = cfg$n_size,
      vjust = 0.9,
      inherit.aes = FALSE,
      angle = 45
    ) +
    geom_segment(
      data = br_df,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      color = "black",
      linewidth = cfg$br_lwd
    ) +
    geom_segment(
      data = br_df,
      aes(x = x1, xend = x1, y = y, yend = y - cfg$br_tick),
      inherit.aes = FALSE,
      color = "black",
      linewidth = cfg$br_lwd
    ) +
    geom_segment(
      data = br_df,
      aes(x = x2, xend = x2, y = y, yend = y - cfg$br_tick),
      inherit.aes = FALSE,
      color = "black",
      linewidth = cfg$br_lwd
    ) +
    geom_text(
      data = br_df,
      aes(x = xmid, y = y_txt, label = p_lab),
      inherit.aes = FALSE,
      color = "black",
      size = cfg$p_size
    ) +
    scale_color_manual(values = c(Normal = "#1b9e77", PFR = "#d95f02", GFR = "#7570b3"), breaks = env_levels) +
    labs(x = "N management strategy", y = y_lab, color = NULL) +
    theme_bw(base_size = cfg$base_size) +
    theme(
      panel.grid = element_blank(),
      legend.position = c(0.70, 0.94),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 5.0),
      legend.key.width = grid::unit(3.2, "mm"),
      legend.key.height = grid::unit(1.5, "mm"),
      axis.line = element_line(color = "black", linewidth = 0.30),
      axis.text.x = element_text(color = "black", size = 6.0),
      axis.text.y = element_text(color = "black", size = 6.0),
      axis.title.x = element_text(size = 6.0),
      axis.title.y = element_text(size = 6.0),
      axis.ticks.x = element_line(color = "black", linewidth = 0.30),
      axis.ticks.y = element_line(color = "black", linewidth = 0.30),
      axis.ticks.length = grid::unit(1.2, "mm"),
      plot.margin = margin(2, 10, 2, 4)
    ) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    scale_y_continuous(breaks = y_breaks)

  ggsave(file.path(output_dir, out_name), p, width = 76, height = 38, units = "mm", dpi = 300)
  invisible(p)
}

plot_model_forest_omnibus <- function(estimates, omnibus, y_lab, out_name, y_breaks, y_limits) {
  strategy_levels <- levels(factor(estimates$N_strategy))
  estimates <- estimates %>%
    mutate(
      N_strategy = factor(N_strategy, levels = strategy_levels),
      env_group = factor(env_group, levels = env_levels),
      n_lab = paste0("n=", k)
    )

  cfg <- list(
    dodge_w = 0.62,
    n_y = y_limits[1] + 0.05 * diff(y_limits),
    p_gap = 0.08 * diff(y_limits),
    p_size = 2.0,
    n_size = 1.9,
    pt_size = 1.35,
    pt_stroke = 0.35,
    err_lwd = 0.25,
    err_width = 0.28,
    base_size = 5.0
  )

  p_df <- estimates %>%
    group_by(N_strategy) %>%
    summarise(ymax = max(response_upr, na.rm = TRUE), .groups = "drop") %>%
    left_join(omnibus, by = "N_strategy") %>%
    mutate(y = ymax + cfg$p_gap)

  p <- ggplot(estimates, aes(x = N_strategy, y = response_estimate, color = env_group)) +
    geom_errorbar(
      aes(ymin = response_lwr, ymax = response_upr),
      position = position_dodge(width = cfg$dodge_w),
      width = cfg$err_width,
      linewidth = cfg$err_lwd
    ) +
    geom_point(
      position = position_dodge(width = cfg$dodge_w),
      shape = 21,
      fill = "white",
      size = cfg$pt_size,
      stroke = cfg$pt_stroke
    ) +
    geom_text(
      aes(x = N_strategy, y = cfg$n_y, label = n_lab, group = env_group),
      position = position_dodge(width = cfg$dodge_w),
      color = "black",
      size = cfg$n_size,
      vjust = 0.9,
      inherit.aes = FALSE,
      angle = 45
    ) +
    geom_text(
      data = p_df,
      aes(x = N_strategy, y = y, label = p_lab),
      inherit.aes = FALSE,
      color = "black",
      size = cfg$p_size
    ) +
    scale_color_manual(values = c(Normal = "#1b9e77", PFR = "#d95f02", GFR = "#7570b3"), breaks = env_levels) +
    labs(x = "N management strategy", y = y_lab, color = NULL) +
    theme_bw(base_size = cfg$base_size) +
    theme(
      panel.grid = element_blank(),
      legend.position = c(0.70, 0.94),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 5.0),
      legend.key.width = grid::unit(3.2, "mm"),
      legend.key.height = grid::unit(1.5, "mm"),
      axis.line = element_line(color = "black", linewidth = 0.30),
      axis.text.x = element_text(color = "black", size = 6.0),
      axis.text.y = element_text(color = "black", size = 6.0),
      axis.title.x = element_text(size = 6.0),
      axis.title.y = element_text(size = 6.0),
      axis.ticks.x = element_line(color = "black", linewidth = 0.30),
      axis.ticks.y = element_line(color = "black", linewidth = 0.30),
      axis.ticks.length = grid::unit(1.2, "mm"),
      plot.margin = margin(2, 10, 2, 4)
    ) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    scale_y_continuous(breaks = y_breaks)

  ggsave(file.path(output_dir, out_name), p, width = 76, height = 38, units = "mm", dpi = 300)
  invisible(p)
}

plot_yield_box <- function(d, yield_res, out_name) {
  strategy_levels <- levels(factor(d$N_strategy))
  d_plot <- d %>%
    filter(is.finite(mean_t)) %>%
    mutate(
      N_strategy = factor(N_strategy, levels = strategy_levels),
      env_group = factor(env_group, levels = env_levels)
    )

  letters_df <- yield_res$contrasts %>%
    filter(!is.na(p_adj_bh)) %>%
    mutate(comp = paste(g1, g2, sep = "-")) %>%
    group_by(N_strategy) %>%
    group_modify(~{
      pv <- .x$p_adj_bh
      names(pv) <- .x$comp
      letters <- multcompView::multcompLetters(pv)$Letters
      data.frame(env_group = names(letters), Letters = unname(letters))
    }) %>%
    ungroup() %>%
    mutate(
      N_strategy = factor(N_strategy, levels = strategy_levels),
      env_group = factor(env_group, levels = env_levels)
    )

  y_pos <- d_plot %>%
    group_by(N_strategy, env_group) %>%
    summarise(y = max(mean_t, na.rm = TRUE) * 1.03, .groups = "drop")

  ann_df <- y_pos %>%
    left_join(letters_df, by = c("N_strategy", "env_group"))

  p <- ggplot(d_plot, aes(x = N_strategy, y = mean_t, fill = env_group)) +
    geom_boxplot(
      position = position_dodge(width = 0.72),
      width = 0.58,
      outlier.shape = NA,
      linewidth = 0.30
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      size = 1.4,
      fill = "white",
      color = "black",
      position = position_dodge(width = 0.72),
      stroke = 0.25
    ) +
    geom_text(
      data = ann_df,
      aes(x = N_strategy, y = y, label = Letters, group = env_group),
      position = position_dodge(width = 0.72),
      inherit.aes = FALSE,
      size = 2.0,
      color = "black"
    ) +
    scale_fill_manual(values = c(Normal = "#1b9e77", PFR = "#d95f02", GFR = "#7570b3"), breaks = env_levels) +
    labs(x = "N management strategy", y = expression("Grain yield under N (Mg ha"^-1*")"), fill = NULL) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      legend.position = c(0.75, 0.90),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.direction = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.key.width = grid::unit(3, "mm"),
      legend.key.height = grid::unit(4, "mm"),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.text.x = element_text(color = "black", hjust = 0.5, size = 8),
      axis.text.y = element_text(color = "black", size = 8),
      axis.ticks.x = element_line(color = "black", linewidth = 0.35),
      axis.ticks.y = element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = grid::unit(1.5, "mm")
    ) +
    coord_cartesian(ylim = c(0, 25), clip = "off") +
    scale_y_continuous(breaks = seq(0, 25, by = 5))

  ggsave(file.path(output_dir, out_name), p, width = 151, height = 45, units = "mm", dpi = 300)
  invisible(p)
}

run_subset <- function(output_prefix, subset_countries, collapse_usa_strategy = FALSE) {
  d_lnrr <- prep_metric("lnRR", subset_countries, collapse_usa_strategy)
  d_aen <- prep_metric("AEN", subset_countries, collapse_usa_strategy)

  write.csv(d_lnrr, file.path(output_dir, paste0(output_prefix, "_N_strategy_within_strategy_data_set.csv")), row.names = FALSE)

  lnrr <- model_env_within_strategy(d_lnrr, "lnRR", "lnrr_pct")
  aen <- model_env_within_strategy(d_aen, "AEN", "identity")
  yield <- yield_wilcox_within_strategy(d_lnrr)

  write.csv(lnrr$estimates, file.path(output_dir, paste0(output_prefix, "_env_estimates_within_N_strategy_lnRR.csv")), row.names = FALSE)
  write.csv(lnrr$contrasts, file.path(output_dir, paste0(output_prefix, "_env_contrasts_within_N_strategy_lnRR.csv")), row.names = FALSE)
  write.csv(lnrr$omnibus, file.path(output_dir, paste0(output_prefix, "_env_omnibus_within_N_strategy_lnRR.csv")), row.names = FALSE)
  write.csv(aen$estimates, file.path(output_dir, paste0(output_prefix, "_env_estimates_within_N_strategy_AEN.csv")), row.names = FALSE)
  write.csv(aen$contrasts, file.path(output_dir, paste0(output_prefix, "_env_contrasts_within_N_strategy_AEN.csv")), row.names = FALSE)
  write.csv(aen$omnibus, file.path(output_dir, paste0(output_prefix, "_env_omnibus_within_N_strategy_AEN.csv")), row.names = FALSE)
  write.csv(yield$summary, file.path(output_dir, paste0(output_prefix, "_env_summary_within_N_strategy_yield.csv")), row.names = FALSE)
  write.csv(yield$contrasts, file.path(output_dir, paste0(output_prefix, "_env_contrasts_within_N_strategy_yield_wilcox.csv")), row.names = FALSE)

  plot_model_forest(
    lnrr$estimates,
    lnrr$contrasts,
    "Relative yield change (%)",
    paste0(output_prefix, "_lnRR_forest_by_N_strategy_env_TukeyP_flip.pdf"),
    seq(0, 180, by = 40),
    c(-10, 180)
  )

  plot_model_forest_omnibus(
    lnrr$estimates,
    lnrr$omnibus,
    "Relative yield change (%)",
    paste0(output_prefix, "_lnRR_forest_by_N_strategy_env_OmnibusP.pdf"),
    seq(0, 180, by = 40),
    c(-10, 180)
  )

  plot_model_forest(
    aen$estimates,
    aen$contrasts,
    "Agronomic efficiency of N (kg grain / kg N)",
    paste0(output_prefix, "_AEN_forest_by_N_strategy_env_TukeyP_flip.pdf"),
    seq(0, 60, by = 20),
    c(-5, 60)
  )

  plot_model_forest_omnibus(
    aen$estimates,
    aen$omnibus,
    "Agronomic efficiency of N (kg grain / kg N)",
    paste0(output_prefix, "_AEN_forest_by_N_strategy_env_OmnibusP.pdf"),
    seq(0, 60, by = 20),
    c(-5, 60)
  )

  plot_yield_box(
    d_lnrr,
    yield,
    paste0(output_prefix, "_box_yield_N_strategyXenv_wilcox_letters.pdf")
  )
}

run_subset("USA_China", c("USA", "China"), collapse_usa_strategy = FALSE)
run_subset("China", "China", collapse_usa_strategy = FALSE)
run_subset("USA", "USA", collapse_usa_strategy = TRUE)

cat("Wrote Extended Data Figure 2 within-N-strategy environment contrasts and plots.\n")
