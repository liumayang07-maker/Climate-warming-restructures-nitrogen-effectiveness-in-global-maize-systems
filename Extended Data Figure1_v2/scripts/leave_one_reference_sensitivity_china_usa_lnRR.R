# Leave-one-reference sensitivity for Extended Data Figure 1 interaction contrasts.
# Purpose: test whether density-specific heat-penalty contrasts are driven by individual references.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(nlme)
  library(emmeans)
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
subset_countries <- c("China", "USA")
q_lo <- 0.05
q_hi <- 0.95

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

dat <- read.csv(file.path(input_dir, "Extended_Data_Figure1_density_model_input.csv"), stringsAsFactors = FALSE)
if (!"plant_density_cluster" %in% names(dat) && "plant_density_3cluster" %in% names(dat)) {
  dat$plant_density_cluster <- dat$plant_density_3cluster
}
if (!"plant_density_3cluster" %in% names(dat) && "plant_density_cluster" %in% names(dat)) {
  dat$plant_density_3cluster <- dat$plant_density_cluster
}

col_plant_density <- if ("plant_density_cluster" %in% names(dat)) {
  "plant_density_cluster"
} else {
  "plant_density_3cluster"
}

prep_lnrr <- function(dat) {
  dat %>%
    mutate(
      response = suppressWarnings(as.numeric(.data[["lnRR"]])),
      env_group = factor(.data[["Env_Type"]], levels = env_levels),
      reference_ID = as.factor(.data[["reference_ID"]]),
      site_year_ID = as.factor(.data[["site_year_ID"]]),
      country = as.factor(.data[["country"]]),
      plant_density = as.factor(.data[[col_plant_density]]),
      region_N_label = as.factor(.data[["region_N_label"]]),
      N_rate = suppressWarnings(as.numeric(.data[["N_rate"]]))
    ) %>%
    filter(
      is.finite(response),
      !is.na(env_group),
      !is.na(reference_ID),
      !is.na(site_year_ID),
      is.finite(N_rate),
      region_N_label %in% "RN",
      country %in% subset_countries,
      !is.na(plant_density)
    ) %>%
    mutate(Nc = N_rate - mean(N_rate, na.rm = TRUE)) %>%
    group_by(env_group) %>%
    mutate(
      q05 = quantile(response, q_lo, na.rm = TRUE),
      q95 = quantile(response, q_hi, na.rm = TRUE),
      keep = response >= q05 & response <= q95
    ) %>%
    ungroup() %>%
    filter(keep) %>%
    mutate(
      env_group = factor(env_group, levels = env_levels),
      plant_density = droplevels(plant_density),
      reference_ID = droplevels(reference_ID),
      site_year_ID = droplevels(site_year_ID)
    )
}

has_complete_cells <- function(d) {
  cells <- d %>%
    count(plant_density, env_group, name = "n") %>%
    tidyr::complete(plant_density, env_group, fill = list(n = 0))
  all(cells$n > 0)
}

fit_penalty_diff <- function(d) {
  if (!has_complete_cells(d)) {
    stop("Missing at least one planting-density x environment cell")
  }

  model <- nlme::lme(
    fixed = response ~ env_group * plant_density + Nc,
    random = ~ 1 | reference_ID/site_year_ID,
    data = d,
    method = "REML",
    na.action = na.omit,
    control = lmeControl(opt = "optim")
  )

  emm <- emmeans(model, ~ env_group | plant_density)
  heat_penalties <- contrast(
    emm,
    method = list(
      Normal_minus_PFR = c(1, -1, 0),
      Normal_minus_GFR = c(1, 0, -1)
    ),
    by = "plant_density",
    adjust = "none"
  )

  as.data.frame(pairs(heat_penalties, by = "contrast", adjust = "tukey")) %>%
    rename(
      density_contrast = contrast1,
      heat_penalty = contrast
    ) %>%
    tidyr::separate(density_contrast, into = c("density_1", "density_2"), sep = " - ", remove = FALSE) %>%
    mutate(
      p_lab = format_p(p.value),
      lower.CL = estimate - qt(0.975, df) * SE,
      upper.CL = estimate + qt(0.975, df) * SE
    )
}

d0 <- prep_lnrr(dat)
refs <- sort(unique(as.character(d0$reference_ID)))

baseline <- fit_penalty_diff(d0) %>%
  mutate(
    run_type = "baseline",
    omitted_reference_ID = NA_character_,
    omitted_rows = 0L,
    remaining_rows = nrow(d0),
    remaining_refs = n_distinct(d0$reference_ID),
    status = "ok"
  )

leave_one <- lapply(refs, function(ref_id) {
  d_i <- d0 %>%
    filter(as.character(reference_ID) != ref_id) %>%
    mutate(
      reference_ID = droplevels(reference_ID),
      site_year_ID = droplevels(site_year_ID),
      plant_density = droplevels(plant_density)
    )

  res <- tryCatch(
    fit_penalty_diff(d_i),
    error = function(e) {
      tibble(
        density_contrast = NA_character_,
        density_1 = NA_character_,
        density_2 = NA_character_,
        heat_penalty = NA_character_,
        estimate = NA_real_,
        SE = NA_real_,
        df = NA_real_,
        t.ratio = NA_real_,
        p.value = NA_real_,
        p_lab = NA_character_,
        lower.CL = NA_real_,
        upper.CL = NA_real_,
        status = paste("failed:", conditionMessage(e))
      )
    }
  )

  if (!"status" %in% names(res)) {
    res$status <- "ok"
  }

  res %>%
    mutate(
      run_type = "leave_one_reference",
      omitted_reference_ID = ref_id,
      omitted_rows = nrow(d0) - nrow(d_i),
      remaining_rows = nrow(d_i),
      remaining_refs = n_distinct(d_i$reference_ID)
    )
}) %>%
  bind_rows()

all_runs <- bind_rows(baseline, leave_one) %>%
  mutate(metric = "lnRR")

baseline_key <- baseline %>%
  select(heat_penalty, density_contrast, baseline_estimate = estimate, baseline_p = p.value)

summary_tbl <- leave_one %>%
  filter(status == "ok") %>%
  left_join(baseline_key, by = c("heat_penalty", "density_contrast")) %>%
  group_by(heat_penalty, density_contrast) %>%
  summarise(
    baseline_estimate = first(baseline_estimate),
    baseline_p = first(baseline_p),
    n_success = dplyr::n(),
    estimate_min = min(estimate, na.rm = TRUE),
    estimate_q25 = quantile(estimate, 0.25, na.rm = TRUE),
    estimate_median = median(estimate, na.rm = TRUE),
    estimate_q75 = quantile(estimate, 0.75, na.rm = TRUE),
    estimate_max = max(estimate, na.rm = TRUE),
    sign_consistency = mean(sign(estimate) == sign(first(baseline_estimate)), na.rm = TRUE),
    p_lt_0.05_fraction = mean(p.value < 0.05, na.rm = TRUE),
    p_lt_0.10_fraction = mean(p.value < 0.10, na.rm = TRUE),
    max_abs_shift = max(abs(estimate - first(baseline_estimate)), na.rm = TRUE),
    max_abs_shift_reference_ID = omitted_reference_ID[which.max(abs(estimate - first(baseline_estimate)))],
    .groups = "drop"
  ) %>%
  mutate(
    baseline_p_lab = format_p(baseline_p),
    p05_lab = sprintf("%.2f", p_lt_0.05_fraction),
    sign_lab = sprintf("%.2f", sign_consistency)
  )

write.csv(all_runs, file.path(output_dir, "USA_China_lnRR_leave_one_reference_heat_penalty_coc_all_runs.csv"), row.names = FALSE)
write.csv(summary_tbl, file.path(output_dir, "USA_China_lnRR_leave_one_reference_heat_penalty_coc_summary.csv"), row.names = FALSE)

plot_df <- summary_tbl %>%
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
    label = paste0("sign=", sign_lab, "; P<0.05=", p05_lab)
  )

p <- ggplot(plot_df, aes(y = density_contrast)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.25, linetype = "dashed") +
  geom_segment(aes(x = estimate_min, xend = estimate_max, yend = density_contrast),
               color = "grey65", linewidth = 0.55) +
  geom_segment(aes(x = estimate_q25, xend = estimate_q75, yend = density_contrast),
               color = "black", linewidth = 0.75) +
  geom_point(aes(x = estimate_median), shape = 21, fill = "white", size = 1.4, stroke = 0.35) +
  geom_point(aes(x = baseline_estimate), shape = 23, fill = "black", color = "black", size = 1.5, stroke = 0.35) +
  geom_text(aes(x = 0.36, label = label), hjust = 1, size = 1.9, color = "black") +
  facet_wrap(~heat_penalty, nrow = 1) +
  labs(x = "Leave-one-reference estimate range (lnRR)", y = NULL) +
  coord_cartesian(xlim = c(-0.38, 0.38), clip = "off") +
  scale_x_continuous(breaks = seq(-0.3, 0.3, by = 0.15)) +
  theme_bw(base_size = 5) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 6),
    axis.line = element_line(color = "black", linewidth = 0.30),
    axis.text.x = element_text(color = "black", size = 6),
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.x = element_text(size = 6),
    axis.ticks = element_line(color = "black", linewidth = 0.30),
    axis.ticks.length = grid::unit(1.2, "mm"),
    plot.margin = margin(2, 10, 2, 4)
  )

ggsave(
  file.path(output_dir, "USA_China_lnRR_leave_one_reference_heat_penalty_coc_summary_forest.pdf"),
  p,
  width = 88,
  height = 34,
  units = "mm",
  dpi = 300
)

cat("Wrote leave-one-reference sensitivity for lnRR heat-penalty contrast-of-contrasts.\n")
