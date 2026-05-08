suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)
  library(multcompView)
  library(ggpubr)
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

df <- read.csv(file.path(input_dir, "adaptation_strategy_dataset.csv"), stringsAsFactors = FALSE)

plot_strategy_distribution <- function(
  data,
  keep_cols,
  reference_group,
  fill_cols,
  distribution_prefix,
  mean_difference_file,
  pdf_name,
  n_y,
  letter_y,
  ylim_range,
  y_breaks
) {
  d_long <- data %>%
    select(all_of(keep_cols)) %>%
    pivot_longer(cols = everything(), names_to = "group", values_to = "value") %>%
    mutate(
      group = factor(group, levels = keep_cols),
      value = as.numeric(value)
    ) %>%
    filter(is.finite(value))

  gh <- d_long %>% games_howell_test(value ~ group)
  write.csv(gh, file.path(output_dir, mean_difference_file), row.names = FALSE)

  gh_sig <- gh %>%
    filter(group1 == reference_group) %>%
    mutate(
      label = paste0(sprintf("%+.2f", estimate), sprintf(" P=%.2f", p.adj))
    ) %>%
    add_xy_position(x = "group", fun = "max", step.increase = 0.8)

  dunn <- d_long %>% dunn_test(value ~ group, p.adjust.method = "BH")
  write.csv(dunn, file.path(output_dir, paste0(distribution_prefix, "_dunn_bh_tests_from_script.csv")), row.names = FALSE)

  p_named <- dunn$p.adj
  names(p_named) <- paste(dunn$group1, dunn$group2, sep = "-")
  letters <- multcompView::multcompLetters(p_named)$Letters

  letters_df <- tibble(
    group = factor(names(letters), levels = keep_cols),
    letter = unname(letters)
  )
  write.csv(letters_df, file.path(output_dir, paste0(distribution_prefix, "_dunn_bh_letters_from_script.csv")), row.names = FALSE)

  letter_pos <- d_long %>%
    group_by(group) %>%
    summarise(y = max(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(y = y + 0.05 * diff(range(d_long$value))) %>%
    left_join(letters_df, by = "group")

  n_df <- d_long %>%
    group_by(group) %>%
    summarise(n = sum(is.finite(value)), .groups = "drop") %>%
    mutate(
      y = n_y,
      n_lab = paste0("n=", n)
    )

  mean_sd_df <- d_long %>%
    group_by(group) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = mean - sd,
      ymax = mean + sd
    )
  write.csv(mean_sd_df, file.path(output_dir, paste0(distribution_prefix, "_yield_summary_from_script.csv")), row.names = FALSE)

  p <- ggplot(d_long, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.50, color = "grey15", linewidth = 0.25) +
    geom_boxplot(
      width = 0.18,
      outlier.alpha = 0.80,
      color = "grey15",
      linewidth = 0.25,
      outlier.size = 0.3
    ) +
    scale_fill_manual(values = fill_cols, breaks = keep_cols) +
    geom_errorbar(
      data = mean_sd_df,
      aes(x = group, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      width = 0.14,
      linewidth = 0.25,
      color = "grey10"
    ) +
    geom_point(
      data = mean_sd_df,
      aes(x = group, y = mean),
      inherit.aes = FALSE,
      shape = 23,
      size = 1.5,
      stroke = 0.35,
      fill = "white",
      color = "grey10"
    ) +
    geom_text(
      data = n_df,
      aes(x = group, y = y, label = n_lab),
      inherit.aes = FALSE,
      size = 2,
      vjust = 0
    ) +
    geom_text(
      data = letter_pos,
      aes(x = group, y = letter_y, label = letter),
      inherit.aes = FALSE,
      vjust = 0,
      size = 3
    ) +
    ggpubr::stat_pvalue_manual(
      gh_sig,
      label = "label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.03,
      bracket.size = 0.25,
      size = 2,
      bracket.nudge.y = -10
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.line = element_line(color = "black", linewidth = 0.25),
      axis.text.x = element_text(color = "black", hjust = 1, size = 8, angle = 45),
      axis.text.y = element_text(color = "black", size = 8),
      axis.ticks.x = element_line(color = "black", linewidth = 0.25),
      axis.ticks.y = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = grid::unit(1.5, "mm")
    ) +
    labs(x = NULL, y = "Grain yield (Mg / ha)") +
    coord_cartesian(ylim = ylim_range, clip = "off") +
    scale_y_continuous(breaks = y_breaks)

  ggsave(file.path(output_dir, pdf_name), p, width = 113, height = 64, units = "mm", dpi = 300)
}

plot_strategy_distribution(
  data = df,
  keep_cols = c("RCP_Normal", "RCP_PFR", "WAP1_PFR", "WAP1P_PFR", "WAP2_PFR", "WAP2P_PFR"),
  reference_group = "RCP_PFR",
  fill_cols = c(
    RCP_Normal = "#B2B8FF",
    RCP_PFR = "#83A0A4",
    WAP1_PFR = "#B2A1C0",
    WAP1P_PFR = "#88A9B5",
    WAP2_PFR = "#F2C666",
    WAP2P_PFR = "#989568"
  ),
  distribution_prefix = "Figure4b",
  mean_difference_file = "Figure4c_games_howell_mean_differences_from_script.csv",
  pdf_name = "Figure4b_pfr_yield_distribution.pdf",
  n_y = 0,
  letter_y = 2,
  ylim_range = c(0, 28),
  y_breaks = seq(0, 28, by = 7)
)

plot_strategy_distribution(
  data = df,
  keep_cols = c("RCP_Normal", "RCP_GFR", "WAP1_GFR", "WAP1P_GFR", "WAP2_GFR", "WAP2P_GFR"),
  reference_group = "RCP_GFR",
  fill_cols = c(
    RCP_Normal = "#B2B8FF",
    RCP_GFR = "#83A0A4",
    WAP1_GFR = "#B2A1C0",
    WAP1P_GFR = "#88A9B5",
    WAP2_GFR = "#F2C666",
    WAP2P_GFR = "#989568"
  ),
  distribution_prefix = "Figure4d",
  mean_difference_file = "Figure4e_games_howell_mean_differences_from_script.csv",
  pdf_name = "Figure4d_gfr_yield_distribution.pdf",
  n_y = -3.8,
  letter_y = -2,
  ylim_range = c(-4, 32),
  y_breaks = seq(0, 32, by = 8)
)
