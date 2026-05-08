library(dplyr)
library(ggplot2)
library(grid)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

root_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
input_dir <- file.path(root_dir, "input_data")
figure_dir <- file.path(root_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

meta_path <- file.path(input_dir, "meta_data_v2.csv")
df <- read.csv(meta_path, stringsAsFactors = FALSE)

x <- df$lnRR
x <- x[is.finite(x)]
n <- length(x)

mean_lnrr <- mean(x)
mean_pct <- 100 * (exp(mean_lnrr) - 1)
hist_info <- hist(x, breaks = 60, plot = FALSE)
ymax <- max(hist_info$counts)

p_overall <- ggplot(data.frame(lnRR = x), aes(x = lnRR)) +
  geom_histogram(bins = 60, fill = "#DCDDDD", colour = "#727171", linewidth = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = mean_lnrr, colour = "red", linewidth = 0.25) +
  annotate(
    "text",
    x = 1.35,
    y = 0.95 * ymax,
    hjust = 0,
    label = sprintf("Mean = %s%%\n(n = %d)", round(mean_pct), n)
  ) +
  coord_cartesian(xlim = c(-0.2, 2.0)) +
  labs(x = "Overall effect (lnRR)", y = "Number of observations") +
  theme_classic(base_size = 6)

print(p_overall)
ggsave(
  filename = file.path(figure_dir, "lnrr_histogram.pdf"),
  plot = p_overall,
  width = 70,
  height = 70,
  units = "mm"
)

sum_env <- df %>%
  group_by(Env_Type) %>%
  summarise(
    n = n(),
    mean_lnRR = mean(lnRR, na.rm = TRUE),
    mean_pct = (exp(mean_lnRR) - 1) * 100,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("Mean = %d%%\n(n = %d)", round(mean_pct), n))

bin_width <- 0.08
base_size <- 7.0
line_width_axis <- 0.3
line_width_bar <- 0.18
line_width_vline <- 0.35
tick_len <- unit(1.0, "mm")

palette_env <- c("Normal" = "#E9D9CA", "PFR" = "#2B586A", "GFR" = "#505592")

p_env <- ggplot(df, aes(x = lnRR)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_histogram(
    aes(fill = Env_Type),
    binwidth = bin_width,
    boundary = 0,
    colour = "grey35",
    linewidth = line_width_bar
  ) +
  geom_vline(
    data = sum_env,
    aes(xintercept = mean_lnRR),
    colour = "#d62728",
    linewidth = line_width_vline
  ) +
  geom_text(
    data = sum_env,
    aes(x = mean_lnRR, y = Inf, label = label),
    vjust = 1.15,
    hjust = 0,
    angle = 45,
    size = 2.2,
    lineheight = 1.0,
    colour = "black"
  ) +
  facet_wrap(~ Env_Type, ncol = 1, scales = "fixed") +
  scale_fill_manual(values = palette_env, guide = "none") +
  labs(x = "Overall effect (lnRR)", y = "Number of observations") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = base_size) +
  theme(
    axis.line = element_line(linewidth = line_width_axis, colour = "black"),
    axis.ticks.length = tick_len,
    axis.ticks = element_line(linewidth = 0.25, colour = "black"),
    axis.title.x = element_text(size = base_size + 0.5, margin = margin(t = 3)),
    axis.title.y = element_text(size = base_size + 0.5, margin = margin(r = 3)),
    axis.text = element_text(size = base_size - 0.2, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = base_size + 0.5, face = "bold"),
    panel.spacing = unit(2.0, "mm"),
    plot.margin = margin(3, 6, 3, 3, unit = "mm")
  )

print(p_env)
ggsave(
  filename = file.path(figure_dir, "lnrr_environment_histogram.pdf"),
  plot = p_env,
  width = 62,
  height = 68,
  units = "mm",
  dpi = 300,
  useDingbats = FALSE
)
