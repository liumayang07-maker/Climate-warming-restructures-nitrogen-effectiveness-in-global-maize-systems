## =========================================================
## Unweighted hierarchical mixed model (nlme::lme)
## Modules A/B/C: compute + plot immediately (no function workflow)
## =========================================================
## -------------------------
## Path & data
## -------------------------
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

root_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
input_dir <- file.path(root_dir, "input_data")
out_dir <- file.path(root_dir, "output_data", "lnrr")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dat <- read.csv(file.path(input_dir, "figure2_model_input.csv"), stringsAsFactors = FALSE)

## -------------------------
## 1) Column mapping
## -------------------------
col_lnRR    <- "lnRR"
col_env     <- "Env_Type"        # Normal/PFR/GFR
col_ref     <- "reference_ID"
col_sy      <- "site_year_ID"
col_EDH_sum <- "EDH35_total"     # EDH_sum column
col_env_ID  <- "env_ID"
col_Nclust <- "N_rate_cluster"
col_Nrate <- "N_rate"

env_levels <- c("Normal","PFR","GFR")
trim_envs  <- c("PFR","GFR")     # only trim these
trial_id   <- col_env_ID             # dedup unit (change to "env_ID" if you have it)

## -------------------------
## 2) Packages
## -------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(lme4)
})

## =========================================================
## 3) Preprocess: build d0 and d_env (keep finite lnRR, IDs)
## =========================================================
d0 <- dat %>%
  mutate(
    lnRR = .data[[col_lnRR]],
    env_group = factor(.data[[col_env]], levels = env_levels),
    reference_ID = as.factor(.data[[col_ref]]),
    site_year_ID = as.factor(.data[[col_sy]]),
    env_ID       = as.factor(.data[[col_env_ID]])
  ) %>%
  filter(
    is.finite(lnRR),
    !is.na(reference_ID),
    !is.na(site_year_ID)
  )

## only rows with env + finite EDH_sum for EDH-QC/trim
d_env <- d0 %>%
  filter(!is.na(env_group)) %>%
  mutate(EDH_sum = .data[[col_EDH_sum]]) %>%
  filter(is.finite(EDH_sum))

cat("\n[INFO] d0 rows =", nrow(d0),
    "| d_env rows (env + finite EDH) =", nrow(d_env), "\n")

## =========================================================
## 4) Trial-level dedup for EDH distribution (avoid overweighting multi-arm trials)
## =========================================================
## consistency check: EDH_sum should be (almost) constant within trial_id
edd_check <- d_env %>%
  group_by(.data[[trial_id]], env_group) %>%
  summarise(
    n_rows = n(),
    n_edd  = n_distinct(EDH_sum),
    edd_min = min(EDH_sum, na.rm = TRUE),
    edd_max = max(EDH_sum, na.rm = TRUE),
    .groups = "drop"
  )

n_inconsistent <- sum(edd_check$n_edd > 1, na.rm = TRUE)
if (n_inconsistent > 0) {
  warning(sprintf(
    "[WARN] %s trial(s) have >1 distinct EDH_sum within the same %s. Using median EDH_sum per trial.",
    n_inconsistent, trial_id
  ))
}

## trial-level dataset: 1 row per trial (per env_group)
trial_df <- d_env %>%
  group_by(.data[[trial_id]], env_group) %>%
  summarise(
    reference_ID = dplyr::first(reference_ID),
    EDH_sum = median(EDH_sum, na.rm = TRUE),
    .groups = "drop"
  )

cat("[INFO] trial_df rows (dedup trials) =", nrow(trial_df), "\n")

## =========================================================
## 5) QC plots: EDH_sum distribution (trial-level; show all envs)
## =========================================================
p_box_edd <- ggplot(trial_df, aes(x = env_group, y = EDH_sum, fill = env_group)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.65) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Env group", y = col_EDH_sum,
       title = "EDH_sum distribution by Env group (trial-level, deduplicated)")
p_box_edd
ggsave(file.path(out_dir, "QC_EDHsum_boxplot_byEnv_trialLevel.pdf"),
       p_box_edd, width = 120, height = 85, units = "mm", dpi = 300)

## =========================================================
## 6) Trim EDH tails (5% & 95%) within PFR/GFR only (trial-level)
##    - Normal group is kept intact
## =========================================================
qs_trim <- trial_df %>%
  filter(as.character(env_group) %in% trim_envs) %>%
  group_by(env_group) %>%
  summarise(
    q05 = quantile(EDH_sum, 0.05, na.rm = TRUE),
    q95 = quantile(EDH_sum, 0.95, na.rm = TRUE),
    n_trials_before = n(),
    .groups = "drop"
  )

## determine keep_trial at trial level
trial_keep <- trial_df %>%
  left_join(qs_trim, by = "env_group") %>%
  mutate(
    keep_trial = dplyr::case_when(
      as.character(env_group) %in% trim_envs ~ (EDH_sum >= q05 & EDH_sum <= q95),
      TRUE ~ TRUE
    )
  )

## kept trial ids for trimmed envs
kept_trials_trimmed <- trial_keep %>%
  filter(as.character(env_group) %in% trim_envs, keep_trial) %>%
  pull(.data[[trial_id]]) %>%
  unique()

## apply to row-level data: keep all Normal; keep only kept trials in PFR/GFR
d_trim <- d_env %>%
  mutate(
    keep_row = dplyr::case_when(
      as.character(env_group) %in% trim_envs ~ (.data[[trial_id]] %in% kept_trials_trimmed),
      TRUE ~ TRUE
    )
  ) %>%
  filter(keep_row) %>%
  select(-keep_row)

cat("\n[QC] EDH trim applied to row-level data.\n",
    "Rows before:", nrow(d_env), "| after:", nrow(d_trim), "\n")

## =========================================================
## 7) Export diagnostics (thresholds + counts)
## =========================================================
## thresholds with n_trials_after
qs_trim_out <- qs_trim %>%
  left_join(
    trial_keep %>%
      filter(as.character(env_group) %in% trim_envs, keep_trial) %>%
      count(env_group, name = "n_trials_after"),
    by = "env_group"
  ) %>%
  mutate(n_trials_after = ifelse(is.na(n_trials_after), 0L, n_trials_after))

count_trials_before <- trial_df %>%
  group_by(env_group) %>%
  summarise(n_trials = n(), .groups = "drop")

count_trials_after <- trial_keep %>%
  filter(keep_trial) %>%
  group_by(env_group) %>%
  summarise(n_trials = n(), .groups = "drop")

count_rows_before <- d_env %>%
  group_by(env_group) %>%
  summarise(k = n(),
            n_ref = n_distinct(reference_ID),
            n_sy  = n_distinct(site_year_ID),
            .groups = "drop")

count_rows_after <- d_trim %>%
  group_by(env_group) %>%
  summarise(k = n(),
            n_ref = n_distinct(reference_ID),
            n_sy  = n_distinct(site_year_ID),
            .groups = "drop")

write.csv(qs_trim_out, file.path(out_dir, "QC_EDHsum_trimThresholds_PFR_GFR_trialLevel.csv"), row.names = FALSE)
write.csv(count_trials_before, file.path(out_dir, "QC_trialCounts_byEnv_beforeTrim.csv"), row.names = FALSE)
write.csv(count_trials_after,  file.path(out_dir, "QC_trialCounts_byEnv_afterTrim.csv"),  row.names = FALSE)
write.csv(count_rows_before,   file.path(out_dir, "QC_rowCounts_byEnv_beforeTrim.csv"),   row.names = FALSE)
write.csv(count_rows_after,    file.path(out_dir, "QC_rowCounts_byEnv_afterTrim.csv"),    row.names = FALSE)

## =========================================================
## 8) Downstream-ready dataset
## =========================================================
d0 <- d_trim

d_env <- d0 %>% filter(!is.na(env_group))

## -------------------------
## Counts: n labels (k, n_ref, n_sy)
## -------------------------
# Overall
n_overall_k   <- nrow(d0)
n_overall_ref <- n_distinct(d0$reference_ID)
n_overall_sy  <- n_distinct(d0$site_year_ID)

# Env
count_env <- d_env %>%
  group_by(env_group) %>%
  summarise(
    k     = n(),
    n_ref = n_distinct(reference_ID),
    n_sy  = n_distinct(site_year_ID),
    .groups = "drop"
  ) %>% mutate(group = as.character(env_group))

## =========================================================
## Part 1) Effect estimation (unweighted hierarchical lme)
## =========================================================

## -------------------------
## Overall mean model
## -------------------------
mO <- nlme::lme(
  fixed  = lnRR ~ 1,
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d0,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

betaO <- as.numeric(fixef(mO)[1])
seO   <- as.numeric(sqrt(vcov(mO)[1,1]))
lwrO  <- betaO - 1.96 * seO
uprO  <- betaO + 1.96 * seO

resO <- tibble(
  group = "Overall",
  est   = betaO,
  se    = seO,
  lwr   = lwrO,
  upr   = uprO,
  k     = n_overall_k,
  n_ref = n_overall_ref,
  n_sy  = n_overall_sy
) %>%
  mutate(
    pct     = (exp(est) - 1) * 100,
    pct_lwr = (exp(lwr) - 1) * 100,
    pct_upr = (exp(upr) - 1) * 100
  )

write.csv(resO, file.path(out_dir, "O_overall_unweighted_lme.csv"), row.names = FALSE)

## -------------------------
## Env mean model: 0 + env_group
## -------------------------
mA <- nlme::lme(
  fixed  = lnRR ~ 0 + env_group,
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d_env,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

betaA <- fixef(mA)
VA    <- vcov(mA)
seA   <- sqrt(diag(VA))

resA <- tibble(
  term = names(betaA),
  est  = as.numeric(betaA),
  se   = as.numeric(seA),
  lwr  = est - 1.96 * se,
  upr  = est + 1.96 * se
) %>%
  mutate(
    group = gsub("^env_group", "", term),
    pct     = (exp(est) - 1) * 100,
    pct_lwr = (exp(lwr) - 1) * 100,
    pct_upr = (exp(upr) - 1) * 100
  ) %>%
  select(group, est, se, lwr, upr, pct, pct_lwr, pct_upr) %>%
  left_join(count_env %>% select(group, k, n_ref, n_sy), by = "group")

write.csv(resA, file.path(out_dir, "A_env_unweighted_lme.csv"), row.names = FALSE)

## Combine Overall + Env
res_all <- bind_rows(
  resO %>% select(group, est, se, lwr, upr, pct, pct_lwr, pct_upr, k, n_ref, n_sy),
  resA %>% select(group, est, se, lwr, upr, pct, pct_lwr, pct_upr, k, n_ref, n_sy)
) %>%
  mutate(
    group = factor(group, levels = c("Overall", env_levels)),
    # n label in plot: choose one
    n_lab = paste0("n=", k)
    # if you want show study counts:
    # n_lab = paste0("k=", k, "\nref=", n_ref, "\nsy=", n_sy)
  )

write.csv(res_all, file.path(out_dir, "OA_overall_plus_env_unweighted_lme.csv"), row.names = FALSE)

## =========================================================
## Part 2) Inference: overall env significance + pairwise P
## =========================================================

## -------------------------
## Overall env significance (approx F test)
## Use model with intercept for inference (env_group as fixed effect)
## -------------------------
m1 <- nlme::lme(
  fixed  = lnRR ~ env_group,     # Normal as reference
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d_env,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

# (a) Approximate F test (Wald-type)
ftab <- anova(m1)
write.csv(as.data.frame(ftab), file.path(out_dir, "ANOVA_env_Ftest_REML.csv"), row.names = FALSE)

# Extract p-value for env_group row
p_env_overall <- as.numeric(ftab["env_group", "p-value"])
# (if rownames differ, fallback)
if (is.na(p_env_overall)) {
  p_env_overall <- as.numeric(ftab[grep("env_group", rownames(ftab))[1], "p-value"])
}

# (b) Likelihood ratio test (ML) as additional confirmation
m0 <- nlme::lme(
  fixed  = lnRR ~ 1,
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d_env,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)
lrt_tab <- anova(update(m0, method="ML"), update(m1, method="ML"))
capture.output(lrt_tab, file = file.path(out_dir, "ANOVA_env_LRT_ML.txt"))

## -------------------------
## Pairwise comparisons (Tukey adjusted)
## -------------------------
emm <- emmeans(m1, ~ env_group)

pw_none  <- pairs(emm, adjust = "none")
pw_tukey <- pairs(emm, adjust = "tukey")

pw_tab_none <- as.data.frame(pw_none) %>%
  mutate(contrast = as.character(contrast)) %>%
  select(contrast, estimate, SE, df, t.ratio, p.value)

pw_tab_tukey <- as.data.frame(pw_tukey) %>%
  mutate(contrast = as.character(contrast)) %>%
  select(contrast, estimate, SE, df, t.ratio, p.value)

pw_both <- bind_rows(
  pw_tab_none  %>% mutate(adjust = "none"),
  pw_tab_tukey %>% mutate(adjust = "tukey")
) %>% select(adjust, contrast, estimate, SE, df, t.ratio, p.value)

write.csv(pw_both, file.path(out_dir, "Pairwise_env_none_and_tukey_lnRR.csv"), row.names = FALSE)


## =========================================================
## Annotation text: overall env P + ordered pairwise P
## =========================================================

## ensure contrast order fixed
## desired order exactly as your table
## ---------------------------------------------------------
## Pairwise annotation text (Tukey) -> bracket style "(...)" 
## ---------------------------------------------------------
## ---- order (top -> bottom) ----
group_levels <- c("Overall", "GFR", "PFR", "Normal")

res_plot <- res_all %>%
  mutate(
    group   = factor(group, levels = group_levels),
    col_key = ifelse(as.character(group) == "Overall", "Overall", as.character(group))
  )

## ---- colors ----
env_cols <- c(
  Normal = "#f8b682",
  PFR    = "#116796",
  GFR    = "#009a76"
)

## ---- style tuned for 30 x 64 mm ----
base_size <- 6.4  # Public-release note
axis_fs   <- 6.2  # Public-release note
tick_fs   <- 6.0  # Public-release note
pt_size   <- 1.55  # Public-release note
pt_stroke <- 0.55  # Public-release note
err_lwd <- 0.32  # Public-release note
cap_lwd <- 0.32  # Public-release note
cap_h   <- 0.10  # Public-release note
n_fs    <- 2.2  # Public-release note
n_x_pad <- 1.2  # Public-release note
vline_lwd <- 0.28  # Public-release note
hsep_lwd  <- 0.28  # Public-release note

right_mar_mm  <- 2.0  # Public-release note
left_mar_mm   <- 1.8  # Public-release note
top_mar_mm    <- 1.2  # Public-release note
bottom_mar_mm <- 1.2  # Public-release note

## ---- x scale (0 -> 90) ----
x_breaks <- c(0, 30, 60, 90)

## x limits: extend a bit to the right for n labels
x_min_data <- min(res_plot$pct_lwr, 0, na.rm = TRUE)
x_max_data <- max(res_plot$pct_upr, 90, na.rm = TRUE)

x_lim_low  <- x_min_data
x_lim_high <- x_max_data + 10   # room for "n=" (tune 6~14)

## ---- separator between GFR and Overall ----
y_sep <- 1.5  # Normal=1, PFR=2, GFR=3, Overall=4

## ---- n label position (to the RIGHT of point) ----
res_plot <- res_plot %>%
  mutate(
    y_num = as.numeric(group),
    n_x   = pct + n_x_pad
  )

p <- ggplot(res_plot, aes(y = group)) +
  
  ## reference line at 0
  geom_vline(xintercept = 0, linetype = 2, linewidth = vline_lwd) +
  
  ## separator line before Overall
  geom_hline(yintercept = y_sep, linetype = 2, linewidth = hsep_lwd) +
  
  ## horizontal CI
  geom_segment(
    aes(x = pct_lwr, xend = pct_upr, yend = group, color = col_key),
    linewidth = err_lwd, show.legend = FALSE
  ) +
  
  ## CI end caps
  geom_segment(
    aes(x = pct_lwr, xend = pct_lwr, y = y_num - cap_h, yend = y_num + cap_h, color = col_key),
    linewidth = cap_lwd, show.legend = FALSE
  ) +
  geom_segment(
    aes(x = pct_upr, xend = pct_upr, y = y_num - cap_h, yend = y_num + cap_h, color = col_key),
    linewidth = cap_lwd, show.legend = FALSE
  ) +
  
  ## hollow points
  geom_point(
    aes(x = pct, color = col_key),
    shape = 21, fill = "white",
    size = pt_size, stroke = pt_stroke,
    show.legend = FALSE
  ) +
  
  ## n labels
  geom_text(
    aes(x = 10, label = n_lab),
    hjust = 0, vjust = 0.35, size = n_fs,
    show.legend = FALSE
  ) +
  
  ## colors
  scale_color_manual(values = c(Overall = "black", env_cols)) +
  
  ## x axis normal direction 0->90
  scale_x_continuous(
    breaks = x_breaks,
    limits = c(x_lim_low, x_lim_high),
    expand = c(0, 0)
  ) +
  
  ## y labels on LEFT (default) + force order
  scale_y_discrete(limits = group_levels) +
  
  theme_bw(base_size = base_size) +
  theme(
    panel.grid = element_blank(),
    
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = tick_fs),
    axis.text.x  = element_text(size = tick_fs),
    axis.title.x = element_text(size = axis_fs, margin = margin(t = 0.8, unit = "mm")),
    
    axis.ticks   = element_line(linewidth = 0.25, color = "black"),
    axis.ticks.length = unit(1.2, "mm"),
    
    plot.margin  = margin(top_mar_mm, right_mar_mm, bottom_mar_mm, left_mar_mm, unit = "mm")
  ) +
  labs(x = "Yield change over control (%)")

print(p)

ggsave(
  file.path(out_dir, "Forest_overall_plus_env_30x64mm.pdf"),
  p, width = 45, height = 65, units = "mm",
  dpi = 300, useDingbats = FALSE
)

cat("\nDone. Outputs saved to: ", normalizePath(out_dir), "\n")

#############################################
# Public-release note
#############################################
## -------------------------

## =========================================================
## 1) Data preprocessing
## =========================================================
d1 <- d0 %>%
  mutate(
    lnRR = .data[[col_lnRR]],
    env_group = factor(.data[[col_env]], levels = env_levels),
    reference_ID = as.factor(.data[[col_ref]]),
    site_year_ID = as.factor(.data[[col_sy]]),
    N_rate_cluster = as.character(.data[[col_Nclust]])
  ) %>%
  filter(
    is.finite(lnRR),
    !is.na(reference_ID),
    !is.na(site_year_ID),
    !is.na(env_group),
    !is.na(N_rate_cluster)
  )

## N_rate_cluster order
Ncl_levels <- c("0_100", "100_200", "200_300", "300_400", "400_600")
Ncl_levels2 <- unique(c(Ncl_levels, sort(unique(d1$N_rate_cluster))))
d1 <- d1 %>% mutate(N_rate_cluster = factor(N_rate_cluster, levels = Ncl_levels2))

## =========================================================
## 2) Counts for n labels (per cell)
## =========================================================
count_cell <- d1 %>%
  group_by(N_rate_cluster, env_group) %>%
  summarise(
    k     = n(),
    n_ref = n_distinct(reference_ID),
    n_sy  = n_distinct(site_year_ID),
    .groups = "drop"
  )

print(count_cell)
write.csv(count_cell, file.path(out_dir, "Cell_counts_Ncluster_by_env.csv"), row.names = FALSE)

## =========================================================
## 3) Cell means estimation (unweighted lme)
##    lnRR ~ 0 + env_group:N_rate_cluster
## =========================================================
cat("\n[Model] Cell means: lnRR ~ 0 + env_group:N_rate_cluster (unweighted lme)\n")

mB <- nlme::lme(
  fixed  = lnRR ~ 0 + env_group:N_rate_cluster,
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d1,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

betaB <- fixef(mB)
VB    <- vcov(mB)
seB   <- sqrt(diag(VB))

resB <- tibble(
  term = names(betaB),
  est  = as.numeric(betaB),
  se   = as.numeric(seB),
  lwr  = est - 1.96 * se,
  upr  = est + 1.96 * se
) %>%
  mutate(
    term2 = gsub("^env_group", "", term),
    env   = sub(":N_rate_cluster.*$", "", term2),
    Ncl   = sub("^.*:N_rate_cluster", "", term2),
    pct     = (exp(est) - 1) * 100,
    pct_lwr = (exp(lwr) - 1) * 100,
    pct_upr = (exp(upr) - 1) * 100
  ) %>%
  mutate(
    env = factor(env, levels = env_levels),
    Ncl = factor(Ncl, levels = levels(d1$N_rate_cluster))
  ) %>%
  select(Ncl, env, est, se, lwr, upr, pct, pct_lwr, pct_upr) %>%
  left_join(
    count_cell %>% transmute(Ncl = N_rate_cluster, env = env_group, k, n_ref, n_sy),
    by = c("Ncl","env")
  ) %>%
  mutate(
    n_lab = paste0("n=", k)
  )

write.csv(resB, file.path(out_dir, "Env_effects_within_Ncluster_unweighted_lme.csv"), row.names = FALSE)

## =========================================================
## 4) Inference (interaction model) + pairwise within each cluster
##    lnRR ~ env_group * N_rate_cluster
## =========================================================
cat("\n[Inference] Interaction model for env_group within each N_rate_cluster\n")

m_int <- nlme::lme(
  fixed  = lnRR ~ env_group * N_rate_cluster,
  random = ~ 1 | reference_ID/site_year_ID,
  data   = d1,
  method = "REML",
  na.action = na.omit,
  control = lmeControl(opt = "optim")
)

ftab_int <- anova(m_int)
print(ftab_int)
write.csv(as.data.frame(ftab_int), file.path(out_dir, "ANOVA_env_byNcluster_Ftest_REML.csv"), row.names = FALSE)

emm_byN <- emmeans(m_int, ~ env_group | N_rate_cluster)

pw_none  <- pairs(emm_byN, adjust = "none")
pw_tukey <- pairs(emm_byN, adjust = "tukey")

pw_tab_none <- as.data.frame(pw_none) %>%
  mutate(
    contrast = as.character(contrast),
    N_rate_cluster = as.character(N_rate_cluster)
  ) %>%
  select(N_rate_cluster, contrast, estimate, SE, df, t.ratio, p.value)

pw_tab_tukey <- as.data.frame(pw_tukey) %>%
  mutate(
    contrast = as.character(contrast),
    N_rate_cluster = as.character(N_rate_cluster)
  ) %>%
  select(N_rate_cluster, contrast, estimate, SE, df, t.ratio, p.value)

write.csv(pw_tab_none,  file.path(out_dir, "Pairwise_env_within_Ncluster_none_lnRR.csv"),  row.names = FALSE)
write.csv(pw_tab_tukey, file.path(out_dir, "Pairwise_env_within_Ncluster_Tukey_lnRR.csv"), row.names = FALSE)

## =========================================================
## 5) Final plot (no facet)
##    X = Ncl, Y = pct, fill = env
##    brackets + P (within cluster) from Tukey table
##    n labels fixed at y=0 line
## =========================================================
cat("\n[Plot] Final plot: N_rate_cluster on X, env as fill, with pairwise P brackets\n")

## resB already has Ncl/env factors; keep them
dodge_width <- 0.72
k_env <- length(env_levels)
offsets <- setNames(((seq_len(k_env) - (k_env + 1)/2) * (dodge_width / k_env)), env_levels)

## bracket base y for each cluster
y_pad_ratio   <- 0.04  # Public-release note
bracket_gap   <- 3.0  # Public-release note
base_offset   <- 1.1  # Public-release note

y_pad <- y_pad_ratio * diff(range(c(resB$pct_lwr, resB$pct_upr), na.rm = TRUE))
y_top_byN <- resB %>%
  group_by(Ncl) %>%
  summarise(y_base = max(pct_upr, 0, na.rm = TRUE) + base_offset*y_pad, .groups = "drop")

## build bracket dataframe (use Tukey; switch to pw_tab_none if needed)
resB <- resB %>%
  mutate(Ncl = trimws(as.character(Ncl))) %>%
  mutate(Ncl = factor(Ncl, levels = unique(Ncl)))

pw_tab_tukey <- pw_tab_tukey %>%
  mutate(N_rate_cluster = trimws(as.character(N_rate_cluster)))


pw_plot <- pw_tab_tukey %>%
  mutate(
    Ncl = factor(N_rate_cluster, levels = levels(resB$Ncl)),
    contrast = as.character(contrast)
  ) %>%
  tidyr::separate(contrast, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    g1 = factor(g1, levels = env_levels),
    g2 = factor(g2, levels = env_levels)
  ) %>%
  left_join(y_top_byN, by = "Ncl") %>%
  group_by(Ncl) %>%
  mutate(
    x_center = as.numeric(Ncl),
    x1 = x_center + offsets[as.character(g1)],
    x2 = x_center + offsets[as.character(g2)],
    x_min = pmin(x1, x2),
    x_max = pmax(x1, x2),
    rank_inN = row_number(),
    y = y_base + (rank_inN - 1) * (bracket_gap * y_pad),
    p_lab = case_when(
      is.na(p.value) ~ "P=NA",
      p.value < 0.001 ~ "<0.001",
      TRUE ~ paste0(sprintf("%.3f", p.value))
    )
  ) %>% ungroup()

## y-limits with bracket space
br_tick_ratio   <- 0.80  # Public-release note
p_text_ratio    <- 1.0  # Public-release note
bracket_headroom <- 1.50  # Public-release note

n_y             <- 5  # Public-release note
n_vjust         <- -0.45  # Public-release note
n_angle         <- 45  # Public-release note

## y-limits with bracket space
y_min <- min(resB$pct_lwr, 0, na.rm = TRUE)
y_max <- max(pw_plot$y + bracket_headroom * y_pad, 130, na.rm = TRUE)

## Manual colors for env
## -------------------------
env_levels <- c("Normal","PFR","GFR")
env_cols <- c(
  Normal = "#f8b682",
  PFR    = "#116796",
  GFR    = "#009a76"
)

p_final <- ggplot(resB, aes(x = Ncl, y = pct, fill = env)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.5) +
  
  ## CI
  geom_errorbar(
    aes(ymin = pct_lwr, ymax = pct_upr, color = env),
    width = 0.2, linewidth = 0.45,
    position = position_dodge(width = dodge_width),
    show.legend = FALSE
  ) +
  
  ## points (hollow)
  geom_point(
    aes(color = env),
    shape = 21, size = 2.0, stroke = 0.55, fill = "white",
    position = position_dodge(width = dodge_width),
    show.legend = FALSE
  ) +
  
  ## n labels (fixed at y = n_y)
  geom_text(
    aes(y = n_y, label = n_lab),
    size = 2.2, vjust = n_vjust, angle = n_angle,
    position = position_dodge(width = dodge_width),
    show.legend = FALSE
  ) +
  
  ## brackets: horizontal line
  geom_segment(
    data = pw_plot,
    aes(x = x_min, xend = x_max, y = y, yend = y),
    inherit.aes = FALSE, linewidth = 0.2
  ) +
  ## brackets: left vertical tick (length controlled by br_tick_ratio)
  geom_segment(
    data = pw_plot,
    aes(x = x_min, xend = x_min, y = y, yend = y - br_tick_ratio * y_pad),
    inherit.aes = FALSE, linewidth = 0.2
  ) +
  ## brackets: right vertical tick (length controlled by br_tick_ratio)
  geom_segment(
    data = pw_plot,
    aes(x = x_max, xend = x_max, y = y, yend = y - br_tick_ratio * y_pad),
    inherit.aes = FALSE, linewidth = 0.2
  ) +
  ## P value text (height controlled by p_text_ratio)
  geom_text(
    data = pw_plot,
    aes(x = (x_min + x_max)/2, y = y + p_text_ratio * y_pad, label = p_lab),
    inherit.aes = FALSE, size = 2.2
  ) +
  
  scale_fill_manual(values = env_cols, breaks = env_levels) +
  scale_color_manual(values = env_cols, breaks = env_levels) +
  
  coord_cartesian(ylim = c(y_min, y_max)) +
  
  ## y-axis ticks (fixed)
  scale_y_continuous(breaks = seq(0, 120, by = 30)) +
  
  theme_bw(base_size = 8.5) +
  theme(
    panel.grid = element_blank(),
    plot.margin = margin(2.5, 2.5, 2.5, 2.5, unit = "mm"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(3.2, "mm"),
    legend.box.spacing = unit(0.8, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    axis.ticks.length = unit(2.0, "mm"),
    axis.ticks = element_line(color = "black", linewidth = 0.50),
    axis.title.x = element_text(size = 8, colour = "black"),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.x  = element_text(size = 7, colour = "black"),
    axis.text.y  = element_text(size = 7, colour = "black"),
    plot.title = element_text(size = 8.5, face = "plain", hjust = 0.5)
  ) +
  labs(
    x = "N rate cluster (kg N/ha)",
    y = "Yield change over control (%)"
  )

print(p_final)

ggsave(
  file.path(out_dir, "Forest_env_by_Ncluster_noFacet_withP_nAtZero.pdf"),
  p_final, width = 110, height = 70, units = "mm", dpi = 300
)


cat("\nDone. Outputs saved to: ", normalizePath(out_dir), "\n")


## =========================================================
## Dose–response by Env (unweighted hierarchical nlme::lme)
## Workflow: read+filter -> compare models (ML) -> refit best (REML)
## -> predict curves (x = original N_rate) -> plot 3 env curves + 95% CI
## (NO N=100/200/300 point evaluation)
## =========================================================


## -------------------------
## Packages
## -------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(nlme)
  library(ggplot2)
  library(splines)
})

## =========================================================
## 1) Data preprocessing (build a common analysis dataset)
## =========================================================
d <- d0 %>%
  transmute(
    lnRR = .data[[col_lnRR]],
    env_group = factor(.data[[col_env]], levels = env_levels),
    N_rate = as.numeric(.data[[col_Nrate]]),
    reference_ID = as.factor(.data[[col_ref]]),
    site_year_ID = as.factor(.data[[col_sy]])
  ) %>%
  filter(
    is.finite(lnRR),
    !is.na(env_group),
    is.finite(N_rate),
    !is.na(reference_ID),
    !is.na(site_year_ID)
  )

## Center N for numerical stability in spline/poly models
N0 <- median(d$N_rate, na.rm = TRUE)
d <- d %>% mutate(Nc = N_rate - N0)

## Use exactly the same dataset for all candidate models
dAIC <- d %>% select(lnRR, env_group, N_rate, Nc, reference_ID, site_year_ID)

cat("\nKept n =", nrow(dAIC),
    "| ref =", n_distinct(dAIC$reference_ID),
    "| site-year =", n_distinct(dAIC$site_year_ID),
    "| N0 (median) =", round(N0, 2), "\n")

write.csv(dAIC, file.path(out_dir, "dAIC_used_for_models.csv"), row.names = FALSE)

## =========================================================
## 2) Model comparison (ML): linear vs quadratic vs spline(df=4)
## =========================================================
cat("\n[Model comparison: ML] linear vs quadratic vs spline(df=4)\n")

m1_ML <- lme(
  lnRR ~ env_group * Nc,
  random = ~ 1 | reference_ID/site_year_ID,
  data = dAIC, method = "ML",
  control = lmeControl(opt = "optim")
)

m2_ML <- lme(
  lnRR ~ env_group * poly(Nc, 2, raw = TRUE),
  random = ~ 1 | reference_ID/site_year_ID,
  data = dAIC, method = "ML",
  control = lmeControl(opt = "optim")
)

m3_ML <- lme(
  lnRR ~ env_group * ns(Nc, df = 4),
  random = ~ 1 | reference_ID/site_year_ID,
  data = dAIC, method = "ML",
  control = lmeControl(opt = "optim")
)

aic_tab <- AIC(m1_ML, m2_ML, m3_ML)
print(aic_tab)
write.csv(as.data.frame(aic_tab), file.path(out_dir, "AIC_model_comparison_ML.csv"))

best_name <- rownames(aic_tab)[which.min(aic_tab$AIC)]
cat("Best model by AIC:", best_name, "\n")

## =========================================================
## 3) Refit best model using REML (for reporting/prediction)
## =========================================================
cat("\n[Refit best model: REML]\n")

m_best_ML <- switch(best_name,
                    "m1_ML" = m1_ML,
                    "m2_ML" = m2_ML,
                    "m3_ML" = m3_ML
)
m_best <- update(m_best_ML, method = "REML")

sink(file.path(out_dir, "Best_model_summary_REML.txt"))
cat("Best model:", best_name, "\n\n")
summary(m_best)
cat("\n--- Wald-type ANOVA (nlme::anova) ---\n")
print(anova(m_best))
sink()
ftab_best <- anova(m_best)
write.csv(ftab_best, file.path(out_dir, "best_model_ANOVA.csv"))

## =========================================================
## 4) Predict curves on original N_rate axis (fixed-effect 95% CI)
## =========================================================
cat("\n[Predict curves]\n")

N_seq <- seq(min(dAIC$N_rate, na.rm = TRUE),
             max(dAIC$N_rate, na.rm = TRUE),
             length.out = 250)

newdat <- expand.grid(
  N_rate = N_seq,
  env_group = factor(env_levels, levels = env_levels)
) %>%
  mutate(Nc = N_rate - N0)

## Build RHS-only formula for model.matrix (avoid 'lnRR not found')
rhs_form <- if (best_name == "m1_ML") {
  ~ env_group * Nc
} else if (best_name == "m2_ML") {
  ~ env_group * poly(Nc, 2, raw = TRUE)
} else {
  ~ env_group * ns(Nc, df = 4)
}

X    <- model.matrix(rhs_form, data = newdat)
beta <- fixef(m_best)
V    <- vcov(m_best)

eta <- as.numeric(X %*% beta)
se  <- sqrt(diag(X %*% V %*% t(X)))

pred <- newdat %>%
  mutate(
    lnRR_hat = eta,
    lnRR_lwr = eta - 1.96 * se,
    lnRR_upr = eta + 1.96 * se,
    pct      = (exp(lnRR_hat) - 1) * 100,
    pct_lwr  = (exp(lnRR_lwr) - 1) * 100,
    pct_upr  = (exp(lnRR_upr) - 1) * 100
  )

write.csv(pred, file.path(out_dir, "Predicted_curves_env_pct.csv"), row.names = FALSE)

## =========================================================
## 5) Plot curves (3 env) with 95% CI ribbon
## =========================================================
cat("\n[Plot curves]\n")
## ---- compact styling params (easy to tweak) ----
base_fs      <- 9.5  # Public-release note
axis_fs      <- 9.0
tick_fs      <- 8.6
legend_fs    <- 8.6
legend_keymm <- 3.6  # Public-release note
line_w       <- 0.55
ribbon_alpha <- 0.25
hline_w      <- 0.45
marg_mm      <- 2.0  # Public-release note

## Public-release note
## Public-release note
fmt_p1 <- function(p){
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}
ftab_best <- anova(m_best)
p_env <- as.numeric(ftab_best["env_group","p-value"])
p_int <- as.numeric(ftab_best[grep("^env_group:", rownames(ftab_best))[1], "p-value"])
anno_text <- paste0("Env P=", fmt_p1(p_env), "\nEnv×N P=", fmt_p1(p_int))

## annotation position inside panel (top-left)
x_anno <- min(pred$N_rate, na.rm = TRUE)
y_anno <- max(pred$pct_upr, na.rm = TRUE)

env_levels <- c("Normal","PFR","GFR")
env_cols <- c(
  Normal = "#f8b682",
  PFR    = "#116796",
  GFR    = "#009a76"
)

p_curve <- ggplot(pred, aes(x = N_rate, y = pct, color = env_group, fill = env_group)) +
  geom_ribbon(aes(ymin = pct_lwr, ymax = pct_upr),
              alpha = ribbon_alpha, colour = NA) +
  geom_line(linewidth = line_w) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = hline_w) +
  
  ## Public-release note
  annotate(
    "label",
    x = x_anno, y = y_anno,
    label = anno_text,
    hjust = 0, vjust = 1,
    size = 2.6,  # Public-release note
    label.size = 0.2,  # Public-release note
    label.padding = unit(1.0, "mm")
  ) +
  
  scale_color_manual(values = env_cols, breaks = env_levels) +
  scale_fill_manual(values = env_cols, breaks = env_levels) +
  
  theme_bw(base_size = base_fs) +
  theme(
    panel.grid = element_blank(),
    
    ## Public-release note
    plot.margin = margin(marg_mm, marg_mm, marg_mm, marg_mm, unit = "mm"),
    
    ## Public-release note
    axis.title.x = element_text(size = axis_fs, margin = margin(t = 1.2, unit = "mm")),
    axis.title.y = element_text(size = axis_fs, margin = margin(r = 1.2, unit = "mm")),
    axis.text    = element_text(size = tick_fs, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks.length = unit(2.0, "mm"),
    
    ## Public-release note
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = legend_fs),
    legend.key.size = unit(legend_keymm, "mm"),
    legend.spacing.x = unit(1.2, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  ) +
  labs(
    x = "N rate (kg/ha)",
    y = "Yield change over control (%)"
  )

print(p_curve)

ggsave(
  file.path(out_dir, "DoseResponse_env_curves_curveOnly_40x30mm.pdf"),
  p_curve, width = 90, height = 75, units = "mm", dpi = 300
)

cat("\nDone. Outputs saved to:", normalizePath(out_dir), "\n")
