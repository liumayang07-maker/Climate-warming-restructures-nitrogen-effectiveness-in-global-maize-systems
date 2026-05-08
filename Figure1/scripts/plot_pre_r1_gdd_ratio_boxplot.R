library(dplyr)
library(ggplot2)
library(forcats)
library(grid)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

root_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
output_dir <- file.path(root_dir, "output_data")
figure_dir <- file.path(root_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(file.path(output_dir, "env_GDD_EDH_results_from_local.csv"), stringsAsFactors = FALSE)

plot_df <- df %>%
  select(region, total_GDD, preR1_GDD) %>%
  filter(!is.na(preR1_GDD)) %>%
  filter(!is.na(total_GDD), total_GDD > 0) %>%
  mutate(preR1_ratio = preR1_GDD / total_GDD) %>%
  filter(!is.na(region), region != "") %>%
  mutate(region = fct_reorder(region, preR1_ratio, .fun = median, .na_rm = TRUE))

plot_df <- bind_rows(
  plot_df,
  plot_df %>% mutate(region = "Total")
) %>%
  mutate(region = fct_relevel(region, "Total", after = Inf))

p <- ggplot(plot_df, aes(x = region, y = preR1_ratio)) +
  geom_boxplot(outlier.alpha = 0.35, linewidth = 0.25) +
  labs(x = "Region", y = "Pre-R1 GDD / total GDD") +
  theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 9,
      angle = 45,
      hjust = 1,
      vjust = 1,
      margin = margin(t = 2)
    ),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 9, margin = margin(t = 4)),
    axis.title.y = element_text(size = 9, margin = margin(r = 4)),
    plot.margin = margin(t = 2, r = 2, b = 8, l = 2, unit = "mm"),
    axis.ticks.length = unit(1.2, "mm"),
    panel.border = element_rect(linewidth = 0.3)
  )

print(p)

ggsave(
  filename = file.path(figure_dir, "preR1_ratio_by_region_boxplot.pdf"),
  plot = p,
  width = 160,
  height = 90,
  units = "mm",
  dpi = 300
)
