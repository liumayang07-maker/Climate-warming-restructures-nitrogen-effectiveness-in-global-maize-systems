# Supplementary contrasts for Figure 4.
# Purpose: compare heat-exposure environments within each planting-density class.

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
script_name <- if (length(script_arg) > 0) {
  basename(sub("^--file=", "", script_arg[1]))
} else {
  "analyze_density_env_contrasts_china_usa.R"
}

root_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
input_dir <- file.path(root_dir, "input_data")
output_dir <- file.path(root_dir, "output_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

env_levels <- c("Normal", "PFR", "GFR")
if (grepl("china_usa", script_name, ignore.case = TRUE)) {
  subset_countries <- c("China", "USA")
  output_prefix <- "USA_China"
} else if (grepl("china", script_name, ignore.case = TRUE)) {
  subset_countries <- "China"
  output_prefix <- "China"
} else if (grepl("usa", script_name, ignore.case = TRUE)) {
  subset_countries <- "USA"
  output_prefix <- "USA"
} else {
  subset_countries <- c("China", "USA")
  output_prefix <- "USA_China"
}
q_lo <- 0.05
q_hi <- 0.95

dat <- read.csv(file.path(input_dir, "Figure4_density_model_input.csv"), stringsAsFactors = FALSE)
if (!"plant_density_cluster" %in% names(dat) && "plant_density_3cluster" %in% names(dat)) {
  dat$plant_density_cluster <- dat$plant_density_3cluster
}
if (!"plant_density_3cluster" %in% names(dat) && "plant_density_cluster" %in% names(dat)) {
  dat$plant_density_3cluster <- dat$plant_density_cluster
}


pick_col <- function(candidates) {
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0) {
    stop("Missing required column. Tried: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

col_plant_density <- pick_col(c("plant_density_cluster", "plant_density_3cluster"))

prep_metric <- function(metric_col) {
  dat %>%
    mutate(
      response = suppressWarnings(as.numeric(.data[[metric_col]])),
      env_group = factor(.data[["Env_Type"]], levels = env_levels),
      reference_ID = as.factor(.data[["reference_ID"]]),
      site_year_ID = as.factor(.data[["site_year_ID"]]),
      continent = as.factor(.data[["continent"]]),
      country = as.factor(.data[["country"]]),
      plant_density = as.factor(.data[[col_plant_density]]),
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
    mutate(Nc = N_rate - mean(N_rate, na.rm = TRUE)) %>%
    filter(region_N_label %in% "RN", country %in% subset_countries) %>%
    filter(!is.na(plant_density), is.finite(response)) %>%
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
      plant_density = droplevels(plant_density)
    )
}

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

env_cols <- c(Normal = "#1b9e77", PFR = "#d95f02", GFR = "#7570b3")

model_contrasts <- function(d, metric_name, response_transform = c("identity", "lnrr_pct")) {
  response_transform <- match.arg(response_transform)
  model <- nlme::lme(
    fixed = response ~ env_group * plant_density + Nc,
    random = ~ 1 | reference_ID/site_year_ID,
    data = d,
    method = "REML",
    na.action = na.omit,
    control = lmeControl(opt = "optim")
  )

  model_summary <- as.data.frame(anova(model)) %>%
    tibble::rownames_to_column("term") %>%
    mutate(metric = metric_name) %>%
    mutate(p_lab = format_p(`p-value`)) %>%
    relocate(metric, term)

  emm <- emmeans(model, ~ env_group | plant_density)
  emm_all <- emmeans(model, ~ env_group * plant_density)
  estimates <- as.data.frame(emm) %>%
    rename(estimate_model = emmean) %>%
    mutate(
      metric = metric_name,
      response_estimate = if (response_transform == "lnrr_pct") {
        (exp(estimate_model) - 1) * 100
      } else {
        estimate_model
      },
      response_lwr = if (response_transform == "lnrr_pct") {
        (exp(lower.CL) - 1) * 100
      } else {
        lower.CL
      },
      response_upr = if (response_transform == "lnrr_pct") {
        (exp(upper.CL) - 1) * 100
      } else {
        upper.CL
      }
    )

  count_cell <- d %>%
    group_by(plant_density, env_group) %>%
    summarise(
      k = dplyr::n(),
      n_ref = n_distinct(reference_ID),
      n_sy = n_distinct(site_year_ID),
      .groups = "drop"
    )

  estimates <- estimates %>%
    left_join(count_cell, by = c("plant_density", "env_group"))

  contrasts <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    tidyr::separate(contrast, into = c("g1", "g2"), sep = " - ", remove = FALSE) %>%
    mutate(
      metric = metric_name,
      p_lab = format_p(p.value)
    ) %>%
    left_join(
      estimates %>%
        select(plant_density, env_group, response_estimate, response_lwr, response_upr, k) %>%
        rename(g1 = env_group, g1_estimate = response_estimate, g1_lwr = response_lwr,
               g1_upr = response_upr, g1_k = k),
      by = c("plant_density", "g1")
    ) %>%
    left_join(
      estimates %>%
        select(plant_density, env_group, response_estimate, response_lwr, response_upr, k) %>%
        rename(g2 = env_group, g2_estimate = response_estimate, g2_lwr = response_lwr,
               g2_upr = response_upr, g2_k = k),
      by = c("plant_density", "g2")
    ) %>%
    mutate(response_difference = g1_estimate - g2_estimate)

  penalty_contrasts <- contrast(
    emm,
    method = list(
      Normal_minus_PFR = c(1, -1, 0),
      Normal_minus_GFR = c(1, 0, -1)
    ),
    by = "plant_density",
    adjust = "none"
  )

  penalty_df <- as.data.frame(penalty_contrasts) %>%
    mutate(
      metric = metric_name,
      p_lab = format_p(p.value)
    )

  penalty_diff <- as.data.frame(pairs(penalty_contrasts, by = "contrast", adjust = "tukey")) %>%
    rename(
      density_contrast = contrast1,
      heat_penalty = contrast
    ) %>%
    mutate(metric = metric_name) %>%
    tidyr::separate(density_contrast, into = c("density_1", "density_2"), sep = " - ", remove = FALSE) %>%
    left_join(
      penalty_df %>%
        select(plant_density, contrast, estimate, SE, df, t.ratio, p.value) %>%
        rename(
          density_1 = plant_density,
          heat_penalty = contrast,
          density_1_penalty = estimate,
          density_1_SE = SE,
          density_1_df = df,
          density_1_t_ratio = t.ratio,
          density_1_p_value = p.value
        ),
      by = c("density_1", "heat_penalty")
    ) %>%
    left_join(
      penalty_df %>%
        select(plant_density, contrast, estimate, SE, df, t.ratio, p.value) %>%
        rename(
          density_2 = plant_density,
          heat_penalty = contrast,
          density_2_penalty = estimate,
          density_2_SE = SE,
          density_2_df = df,
          density_2_t_ratio = t.ratio,
          density_2_p_value = p.value
        ),
      by = c("density_2", "heat_penalty")
    ) %>%
    mutate(
      p_lab = format_p(p.value),
      interpretation = ifelse(
        estimate < 0,
        paste0(density_1, " penalty smaller than ", density_2),
        paste0(density_1, " penalty larger than ", density_2)
      )
    )

  list(
    model = model,
    estimates = estimates,
    contrasts = contrasts,
    heat_penalties = penalty_df,
    penalty_diff = penalty_diff,
    model_summary = model_summary
  )
}

yield_wilcox_within_density <- function(d) {
  pairs_env <- list(c("Normal", "PFR"), c("Normal", "GFR"), c("PFR", "GFR"))

  d_yield <- d %>%
    filter(is.finite(mean_t)) %>%
    mutate(
      env_group = factor(env_group, levels = env_levels),
      plant_density = droplevels(plant_density)
    )

  summary <- d_yield %>%
    group_by(plant_density, env_group) %>%
    summarise(
      k = dplyr::n(),
      mean = mean(mean_t, na.rm = TRUE),
      median = median(mean_t, na.rm = TRUE),
      q25 = quantile(mean_t, 0.25, na.rm = TRUE),
      q75 = quantile(mean_t, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  contrasts <- d_yield %>%
    group_by(plant_density) %>%
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
        select(plant_density, env_group, mean, median, k) %>%
        rename(g1 = env_group, g1_mean = mean, g1_median = median, g1_k = k),
      by = c("plant_density", "g1")
    ) %>%
    left_join(
      summary %>%
        select(plant_density, env_group, mean, median, k) %>%
        rename(g2 = env_group, g2_mean = mean, g2_median = median, g2_k = k),
      by = c("plant_density", "g2")
    ) %>%
    mutate(mean_difference = g1_mean - g2_mean)

  list(summary = summary, contrasts = contrasts)
}

plot_model_forest <- function(estimates, contrasts, metric_name, y_lab, out_name, y_breaks, y_limits) {
  density_levels <- levels(factor(estimates$plant_density))
  estimates <- estimates %>%
    mutate(
      plant_density = factor(plant_density, levels = density_levels),
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
    plant_density = factor(density_levels, levels = density_levels),
    x_base = seq_along(density_levels)
  )

  y_top <- estimates %>%
    group_by(plant_density) %>%
    summarise(ymax = max(response_upr, na.rm = TRUE), .groups = "drop")

  br_df <- contrasts %>%
    mutate(
      plant_density = factor(plant_density, levels = density_levels),
      g1 = factor(g1, levels = env_levels),
      g2 = factor(g2, levels = env_levels)
    ) %>%
    left_join(x_map, by = "plant_density") %>%
    left_join(env_offset, by = c("g1" = "env_group")) %>%
    rename(x1_off = x_off) %>%
    left_join(env_offset, by = c("g2" = "env_group")) %>%
    rename(x2_off = x_off) %>%
    left_join(y_top, by = "plant_density") %>%
    group_by(plant_density) %>%
    mutate(
      br_rank = row_number(),
      y = ymax + cfg$br_base_gap + (br_rank - 1) * cfg$br_step,
      x1 = x_base + x1_off,
      x2 = x_base + x2_off,
      xmid = (x1 + x2) / 2,
      y_txt = y + cfg$p_text_dy
    ) %>%
    ungroup()

  p <- ggplot(estimates, aes(x = plant_density, y = response_estimate, color = env_group)) +
    geom_errorbar(
      aes(ymin = response_lwr, ymax = response_upr),
      position = position_dodge(width = cfg$dodge_w),
      width = cfg$err_width,
      linewidth = cfg$err_lwd
    ) +
    geom_point(
      position = position_dodge(width = cfg$dodge_w),
      shape = 21, fill = "white",
      size = cfg$pt_size, stroke = cfg$pt_stroke
    ) +
    geom_text(
      aes(x = plant_density, y = cfg$n_y, label = n_lab, group = env_group),
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
    labs(x = "Planting density", y = y_lab, color = NULL) +
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

plot_penalty_diff_forest <- function(penalty_diff, out_name, x_lab, x_breaks, x_limits) {
  df <- penalty_diff %>%
    mutate(
      heat_penalty = factor(
        heat_penalty,
        levels = c("Normal_minus_PFR", "Normal_minus_GFR"),
        labels = c("Normal - PFR", "Normal - GFR")
      ),
      density_contrast = factor(
        density_contrast,
        levels = c("High - Low", "Conventional - Low", "Conventional - High")
      ),
      lower.CL = estimate - qt(0.975, df) * SE,
      upper.CL = estimate + qt(0.975, df) * SE,
      p_text = paste0("P=", p_lab)
    )

  cfg <- list(
    base_size = 5.0,
    pt_size = 1.35,
    pt_stroke = 0.35,
    err_lwd = 0.25,
    zero_lwd = 0.25,
    p_size = 1.9
  )

  p_x <- x_limits[2] - 0.02 * diff(x_limits)

  p <- ggplot(df, aes(y = density_contrast, x = estimate)) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = cfg$zero_lwd, linetype = "dashed") +
    geom_segment(
      aes(x = lower.CL, xend = upper.CL, yend = density_contrast),
      color = "black",
      linewidth = cfg$err_lwd
    ) +
    geom_point(shape = 21, fill = "white", color = "black", size = cfg$pt_size, stroke = cfg$pt_stroke) +
    geom_text(
      aes(x = p_x, label = p_text),
      hjust = 1,
      color = "black",
      size = cfg$p_size
    ) +
    facet_wrap(~heat_penalty, nrow = 1) +
    labs(x = x_lab, y = NULL) +
    theme_bw(base_size = cfg$base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.0),
      axis.line = element_line(color = "black", linewidth = 0.30),
      axis.text.x = element_text(color = "black", size = 6.0),
      axis.text.y = element_text(color = "black", size = 6.0),
      axis.title.x = element_text(size = 6.0),
      axis.ticks.x = element_line(color = "black", linewidth = 0.30),
      axis.ticks.y = element_line(color = "black", linewidth = 0.30),
      axis.ticks.length = grid::unit(1.2, "mm"),
      plot.margin = margin(2, 8, 2, 4)
    ) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    scale_x_continuous(breaks = x_breaks)

  ggsave(file.path(output_dir, out_name), p, width = 76, height = 34, units = "mm", dpi = 300)
  invisible(p)
}

plot_yield_box <- function(d, yield_res, out_name) {
  density_levels <- levels(factor(d$plant_density))
  d_plot <- d %>%
    filter(is.finite(mean_t)) %>%
    mutate(
      plant_density = factor(plant_density, levels = density_levels),
      env_group = factor(env_group, levels = env_levels)
    )

  letters_df <- yield_res$contrasts %>%
    filter(!is.na(p_adj_bh)) %>%
    mutate(comp = paste(g1, g2, sep = "-")) %>%
    group_by(plant_density) %>%
    group_modify(~{
      pv <- .x$p_adj_bh
      names(pv) <- .x$comp
      letters <- multcompView::multcompLetters(pv)$Letters
      data.frame(env_group = names(letters), Letters = unname(letters))
    }) %>%
    ungroup() %>%
    mutate(
      plant_density = factor(plant_density, levels = density_levels),
      env_group = factor(env_group, levels = env_levels)
    )

  y_pos <- d_plot %>%
    group_by(plant_density, env_group) %>%
    summarise(y = max(mean_t, na.rm = TRUE) * 1.03, .groups = "drop")

  ann_df <- y_pos %>%
    left_join(letters_df, by = c("plant_density", "env_group"))

  p <- ggplot(d_plot, aes(x = plant_density, y = mean_t, fill = env_group)) +
    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = 0.5,
      outlier.alpha = 0.25,
      outlier.size = 0.5,
      outlier.shape = 16,
      outlier.stroke = 0.2,
      outlier.colour = "black"
    ) +
    stat_summary(
      aes(group = env_group),
      fun = mean,
      geom = "point",
      position = position_dodge(width = 0.8),
      shape = 23, fill = NA,
      size = 1.0, stroke = 0.4,
      color = "black"
    ) +
    geom_text(
      data = ann_df,
      aes(x = plant_density, y = y, label = Letters, group = env_group),
      position = position_dodge(width = 0.8),
      vjust = 0,
      size = 2.5,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = env_cols, breaks = env_levels) +
    labs(
      x = "Planting density",
      y = "Yield under N fertilization (Mg / ha)",
      fill = NULL
    ) +
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

build_input_qc <- function(d, metric_name) {
  overall <- tibble::tibble(
    metric = metric_name,
    n_rows = nrow(d),
    n_ref = dplyr::n_distinct(d$reference_ID),
    n_sy = dplyr::n_distinct(d$site_year_ID),
    n_country = dplyr::n_distinct(d$country),
    countries = paste(sort(unique(as.character(d$country))), collapse = "|")
  )

  by_cell <- d %>%
    dplyr::count(country, env_group, plant_density, name = "n") %>%
    dplyr::mutate(metric = metric_name) %>%
    dplyr::relocate(metric)

  list(overall = overall, by_cell = by_cell)
}

d_lnrr <- prep_metric("lnRR")
d_aen <- prep_metric("AEN")

qc_lnrr <- build_input_qc(d_lnrr, "lnRR")
qc_aen <- build_input_qc(d_aen, "AEN")
write.csv(
  dplyr::bind_rows(qc_lnrr$overall, qc_aen$overall),
  file.path(output_dir, paste0(output_prefix, "_model_input_qc_overall.csv")),
  row.names = FALSE
)
write.csv(
  dplyr::bind_rows(qc_lnrr$by_cell, qc_aen$by_cell),
  file.path(output_dir, paste0(output_prefix, "_model_input_qc_cell_counts.csv")),
  row.names = FALSE
)

lnrr <- model_contrasts(d_lnrr, "lnRR", "lnrr_pct")
aen <- model_contrasts(d_aen, "AEN", "identity")
yield <- yield_wilcox_within_density(d_lnrr)

write.csv(lnrr$estimates, file.path(output_dir, paste0(output_prefix, "_env_estimates_within_density_lnRR.csv")), row.names = FALSE)
write.csv(lnrr$contrasts, file.path(output_dir, paste0(output_prefix, "_env_contrasts_within_density_lnRR.csv")), row.names = FALSE)
write.csv(lnrr$heat_penalties, file.path(output_dir, paste0(output_prefix, "_heat_penalties_by_density_lnRR.csv")), row.names = FALSE)
write.csv(lnrr$penalty_diff, file.path(output_dir, paste0(output_prefix, "_heat_penalty_contrast_of_contrasts_lnRR.csv")), row.names = FALSE)
write.csv(lnrr$model_summary, file.path(output_dir, paste0(output_prefix, "_model_anova_lnRR.csv")), row.names = FALSE)
write.csv(aen$estimates, file.path(output_dir, paste0(output_prefix, "_env_estimates_within_density_AEN.csv")), row.names = FALSE)
write.csv(aen$contrasts, file.path(output_dir, paste0(output_prefix, "_env_contrasts_within_density_AEN.csv")), row.names = FALSE)
write.csv(aen$heat_penalties, file.path(output_dir, paste0(output_prefix, "_heat_penalties_by_density_AEN.csv")), row.names = FALSE)
write.csv(aen$penalty_diff, file.path(output_dir, paste0(output_prefix, "_heat_penalty_contrast_of_contrasts_AEN.csv")), row.names = FALSE)
write.csv(aen$model_summary, file.path(output_dir, paste0(output_prefix, "_model_anova_AEN.csv")), row.names = FALSE)
write.csv(dplyr::bind_rows(lnrr$model_summary, aen$model_summary), file.path(output_dir, paste0(output_prefix, "_model_anova_summary.csv")), row.names = FALSE)
write.csv(yield$summary, file.path(output_dir, paste0(output_prefix, "_env_summary_within_density_yield.csv")), row.names = FALSE)
write.csv(yield$contrasts, file.path(output_dir, paste0(output_prefix, "_env_contrasts_within_density_yield_wilcox.csv")), row.names = FALSE)

plot_model_forest(
  lnrr$estimates,
  lnrr$contrasts,
  "lnRR",
  "Relative yield change (%)",
  paste0(output_prefix, "_forest_pct_by_density_env_TukeyP_flip.pdf"),
  seq(0, 250, by = 60),
  c(-15, 250)
)

plot_model_forest(
  aen$estimates,
  aen$contrasts,
  "AEN",
  "Agronomic efficiency of N (kg grain / kg N)",
  paste0(output_prefix, "_AEN_forest_by_density_env_TukeyP_flip.pdf"),
  seq(0, 60, by = 20),
  c(-5, 60)
)

plot_penalty_diff_forest(
  lnrr$penalty_diff,
  paste0(output_prefix, "_lnRR_heat_penalty_contrast_of_contrasts_forest.pdf"),
  "Difference in heat-risk penalty (lnRR)",
  seq(-0.3, 0.3, by = 0.15),
  c(-0.36, 0.36)
)

plot_yield_box(
  d_lnrr,
  yield,
  paste0(output_prefix, "_box_yield_densityXenv_wilcox_letters.pdf")
)

cat("Wrote ", output_prefix, " within-density environment contrasts, contrast-of-contrasts and plots for lnRR, AEN and yield under N.\n", sep = "")
cat("QC saved: ", output_prefix, "_model_input_qc_overall.csv and ", output_prefix, "_model_input_qc_cell_counts.csv\n", sep = "")
