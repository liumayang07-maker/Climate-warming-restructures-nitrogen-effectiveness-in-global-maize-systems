library(data.table)
library(ggplot2)
library(patchwork)
library(grid)
library(ggbreak)
library(scales)

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

dt <- fread(file.path(output_dir, "env_EDH_by_relGDD_100bins.csv"))
value_col <- "EDH35_per_hour"
stopifnot(value_col %in% names(dt))

dt[, val := as.numeric(get(value_col))]
dt[is.na(val), val := 0]
dt[, bin := as.integer(bin)]

min_bin <- min(dt$bin, na.rm = TRUE)
max_bin <- max(dt$bin, na.rm = TRUE)
if (min_bin == 0L && max_bin == 99L) {
  dt[, bin := bin + 1L]
  min_bin <- 1L
  max_bin <- 100L
}
if (!(min_bin == 1L && max_bin == 100L)) {
  stop(sprintf("Unexpected bin range: [%d, %d]. Expected 1..100 or 0..99.", min_bin, max_bin))
}

mat_dt <- dcast(dt, env_ID ~ bin, value.var = "val", fill = 0)
bin_cols <- as.character(1:100)
missing_bin_cols <- setdiff(bin_cols, names(mat_dt))
if (length(missing_bin_cols) > 0) {
  for (cc in missing_bin_cols) mat_dt[, (cc) := 0]
}
setcolorder(mat_dt, c("env_ID", bin_cols))

envs <- mat_dt$env_ID
mat <- as.matrix(mat_dt[, ..bin_cols])
rownames(mat) <- envs

row_sum <- rowSums(mat)
zero_env <- names(row_sum[row_sum == 0])
nz_env <- names(row_sum[row_sum > 0])

cat("Total environments:", length(envs), "\n")
cat("All-zero environments:", length(zero_env), "\n")
cat("Non-zero environments:", length(nz_env), "\n")

dt_nz <- dt[env_ID %in% nz_env]
mat_nz <- mat[nz_env, , drop = FALSE]

eps <- 1e-4
mat0 <- mat_nz
mat0[mat0 == 0] <- eps

sds <- apply(mat0, 1, sd)
const_idx <- which(!is.finite(sds) | sds == 0)
if (length(const_idx) > 0) {
  ramp <- seq_len(ncol(mat0)) * 1e-10
  mat0[const_idx, ] <- mat0[const_idx, ] + rep(ramp, each = length(const_idx))
}

cor_mat <- cor(t(mat0), use = "pairwise.complete.obs", method = "pearson")
cor_mat[!is.finite(cor_mat)] <- 0
diag(cor_mat) <- 1

hc <- hclust(as.dist(1 - cor_mat), method = "average")
ord <- hc$labels[hc$order]
dt_nz[, env_ID := factor(env_ID, levels = ord)]

cap <- as.numeric(quantile(dt_nz$val, 0.99, na.rm = TRUE))
dt_nz[, val_cap := pmin(val, cap)]

k <- 2
cl_nz <- cutree(hc, k = k)
dt_nz[, cluster := factor(cl_nz[as.character(env_ID)], levels = sort(unique(cl_nz)))]

cluster_colors <- c("1" = "#2B5A6C", "2" = "#53589A")
anno_dt <- data.table(
  env_ID = factor(ord, levels = ord),
  x = 1L,
  cluster = factor(cl_nz[ord], levels = levels(dt_nz$cluster))
)

base_sz <- 7.5
axis_sz <- 7.0
leg_txt <- 7.0
leg_title <- 7.5
leg_key_mm <- 3.0
leg_key_h_mm <- 3.8
m_strip <- margin(1.2, 0.0, 1.2, 1.2, unit = "mm")
m_heat <- margin(1.2, 1.2, 1.2, 0.0, unit = "mm")

p_strip <- ggplot(anno_dt, aes(x = x, y = env_ID, fill = cluster)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = cluster_colors, name = sprintf("Cluster (k=%d)", k)) +
  theme_void(base_size = base_sz) +
  theme(
    legend.position = "left",
    legend.title = element_text(size = leg_title),
    legend.text = element_text(size = leg_txt),
    legend.key.size = unit(leg_key_mm, "mm"),
    legend.spacing.y = unit(0.8, "mm"),
    legend.box.spacing = unit(0.6, "mm"),
    plot.margin = m_strip
  )

p_heat <- ggplot(dt_nz, aes(x = bin, y = env_ID, fill = val)) +
  geom_tile(width = 1, height = 1) +
  scale_x_continuous(breaks = seq(10, 100, 10), expand = c(0, 0)) +
  geom_vline(xintercept = 48, linetype = "dashed", linewidth = 0.4, color = "black") +
  geom_vline(xintercept = 58, linetype = "dashed", linewidth = 0.4, color = "black") +
  coord_cartesian(xlim = c(0.5, 100.5), clip = "on") +
  scale_fill_gradient(
    low = "white",
    high = "firebrick3",
    limits = c(0, 1),
    oob = scales::squish,
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    name = value_col,
    na.value = "grey90"
  ) +
  labs(x = "Relative thermal progress (%)", y = NULL) +
  theme_bw(base_size = base_sz) +
  theme(
    axis.text.x = element_text(size = axis_sz),
    axis.title.x = element_text(size = axis_sz),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.5),
    legend.title = element_text(size = leg_title),
    legend.text = element_text(size = leg_txt),
    legend.key.width = unit(leg_key_mm, "mm"),
    legend.key.height = unit(leg_key_h_mm, "mm"),
    legend.spacing.y = unit(0.8, "mm"),
    plot.margin = m_heat
  )

p_heatmap <- p_strip + p_heat +
  plot_layout(widths = c(0.05, 0.95), guides = "collect") &
  theme(
    legend.box = "vertical",
    legend.margin = margin(0, 0, 0, 0, unit = "mm")
  )

print(p_heatmap)
ggsave(
  filename = file.path(figure_dir, "EDH_heatmap_cluster_strip.pdf"),
  plot = p_heatmap,
  width = 83,
  height = 85,
  units = "mm",
  device = cairo_pdf
)

fwrite(data.table(env_ID = ord), file.path(output_dir, "env_nonzero_cor_average_order.csv"))
fwrite(data.table(env_ID = zero_env), file.path(output_dir, "env_allzero_list.csv"))

sum_dt <- dt_nz[, .(
  mean_val = mean(val, na.rm = TRUE),
  sd_val = sd(val, na.rm = TRUE),
  n_env = uniqueN(env_ID)
), by = .(cluster, bin)]
sum_dt[, se_val := sd_val / sqrt(pmax(n_env, 1))]

p_curve <- ggplot(sum_dt, aes(x = bin, y = mean_val, color = cluster, group = cluster)) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val, fill = cluster), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.45) +
  scale_color_manual(values = cluster_colors, name = "Cluster") +
  scale_fill_manual(values = cluster_colors, name = "Cluster") +
  scale_x_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = "Relative thermal progress (bins)", y = sprintf("Mean %s", value_col)) +
  theme_bw(base_size = 6.6) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.4),
    axis.title = element_text(size = 6.6),
    axis.text = element_text(size = 6.6),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 6.4),
    legend.text = element_text(size = 6.0),
    legend.key.size = unit(2.6, "mm"),
    legend.spacing.y = unit(0.6, "mm"),
    legend.box.spacing = unit(0.6, "mm"),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

print(p_curve)
ggsave(
  filename = file.path(figure_dir, "EDH_curve.pdf"),
  plot = p_curve,
  width = 65,
  height = 40,
  units = "mm",
  device = cairo_pdf
)

all_env <- unique(dt$env_ID)
env_grp <- data.table(env_ID = all_env)
env_grp[, group := fifelse(
  env_ID %in% zero_env,
  "All-zero",
  fifelse(env_ID %in% names(cl_nz), paste0("Cluster ", cl_nz[env_ID]), NA_character_)
)]

if (any(is.na(env_grp$group))) {
  warning(sprintf(
    "%d environments were not assigned to any of the three groups. Check zero_env and cl_nz coverage.",
    sum(is.na(env_grp$group))
  ))
}

env_grp[, group := factor(group, levels = c("All-zero", "Cluster 1", "Cluster 2"))]
fwrite(env_grp, file.path(output_dir, "env_Type.csv"))

pie_dt <- env_grp[!is.na(group), .(n = .N), by = group][order(group)]
pie_dt[, pct := 100 * n / sum(n)]
pie_dt[, label := sprintf("%.1f%%\n(n=%d)", pct, n)]

pie_cols <- c("All-zero" = "#EBDBCB", "Cluster 1" = "#2B5A6C", "Cluster 2" = "#53589A")

p_pie <- ggplot(pie_dt, aes(x = "", y = n, fill = group)) +
  geom_col(width = 1, color = "black", linewidth = 0.35) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 2.2) +
  scale_fill_manual(values = pie_cols, name = "Environment group") +
  theme_void(base_size = 6.5) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6),
    legend.key.size = unit(2.6, "mm"),
    legend.key.height = unit(2.6, "mm"),
    legend.key.width = unit(2.6, "mm"),
    legend.spacing.y = unit(0.8, "mm"),
    legend.box.spacing = unit(0.5, "mm"),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

ggsave(
  filename = file.path(figure_dir, "environment_group_pie.pdf"),
  plot = p_pie,
  width = 42,
  height = 39,
  units = "mm",
  device = cairo_pdf
)

theme_aux <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.5),
      plot.title = element_text(size = base_size, face = "plain"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      plot.margin = margin(3, 3, 3, 3, unit = "mm")
    )
}

cat("\n========== Global summary: val ==========\n")
print(summary(dt_nz$val))
cat(
  "\nN total:", nrow(dt_nz),
  "\nN NA:", sum(is.na(dt_nz$val)),
  "\nN zero:", sum(dt_nz$val == 0, na.rm = TRUE),
  "\nN >0:", sum(dt_nz$val > 0, na.rm = TRUE), "\n"
)

qs <- c(0, .01, .05, .10, .25, .50, .75, .90, .95, .99, 1)
q_val <- as.data.table(t(quantile(dt_nz$val, probs = qs, na.rm = TRUE)))
q_cap <- as.data.table(t(quantile(dt_nz$val_cap, probs = qs, na.rm = TRUE)))
setnames(q_val, paste0("q", qs * 100))
q_val[, metric := "val"]
setnames(q_cap, paste0("q", qs * 100))
q_cap[, metric := "val_cap"]
q_dt <- rbindlist(list(q_val, q_cap), fill = TRUE)
fwrite(q_dt, file.path(output_dir, "heatmap_value_quantiles.csv"))

p_hist <- ggplot(dt_nz, aes(x = val)) +
  geom_histogram(bins = 60, color = "black", linewidth = 0.25) +
  scale_y_break(c(1500, 43000), scales = 0.25) +
  labs(x = "EDH intensity per hour in each GDD bin (°C)", y = "Count") +
  theme_aux(base_size = 10)

ggsave(
  filename = file.path(figure_dir, "EDH_intensity_hist.pdf"),
  plot = p_hist,
  width = 85,
  height = 72,
  units = "mm",
  device = cairo_pdf
)

p_bin_box <- ggplot(dt_nz, aes(x = factor(bin), y = val)) +
  geom_boxplot(outlier.size = 0.15, linewidth = 0.25) +
  scale_x_discrete(breaks = as.character(seq(10, 100, 10))) +
  labs(x = "Relative thermal progress (%)", y = "EDH intensity per hour in each GDD bin (°C)") +
  theme_aux(base_size = 9)

ggsave(
  filename = file.path(figure_dir, "EDH_intensity_bin_box.pdf"),
  plot = p_bin_box,
  width = 85,
  height = 72,
  units = "mm",
  device = cairo_pdf
)
