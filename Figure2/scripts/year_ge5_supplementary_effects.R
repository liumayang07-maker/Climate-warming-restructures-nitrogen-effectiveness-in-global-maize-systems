## =========================================================
## Supplementary Figure 2 analysis:
## effects within long-term experiments (Exp_year_duration >= 5)
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
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

root_dir <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
input_file <- file.path(root_dir, "meta_data_v2_filled_Data_year_duration.csv")
out_dir <- file.path(root_dir, "Figure2", "output_data", "year_ge5")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dat <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

env_levels <- c("Normal", "PFR", "GFR")
ncl_levels <- c("0_100", "100_200", "200_300", "300_400", "400_600")

required_cols <- c(
  "reference_ID", "site_year_ID", "Env_Type", "N_rate_cluster",
  "Exp_year_duration", "Data_year_duration", "lnRR", "AEN"
)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

fit_model <- function(fixed_formula, data) {
  attempts <- list(
    list(type = "lme_reference_ID_site_year_ID", expr = quote(lme(
      fixed = fixed_formula,
      random = ~ 1 | reference_ID/site_year_ID,
      data = data,
      method = "REML",
      na.action = na.omit,
      control = lmeControl(opt = "optim")
    ))),
    list(type = "lme_reference_ID", expr = quote(lme(
      fixed = fixed_formula,
      random = ~ 1 | reference_ID,
      data = data,
      method = "REML",
      na.action = na.omit,
      control = lmeControl(opt = "optim")
    ))),
    list(type = "lm_no_random_effect", expr = quote(lm(
      formula = fixed_formula,
      data = data,
      na.action = na.omit
    )))
  )

  errors <- character()
  for (a in attempts) {
    fit <- try(eval(a$expr), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      return(list(fit = fit, model_type = a$type, errors = errors))
    }
    errors <- c(errors, paste0(a$type, ": ", conditionMessage(attr(fit, "condition"))))
  }
  stop("All model attempts failed:\n", paste(errors, collapse = "\n"))
}

d_base <- dat %>%
  mutate(
    Data_year_duration = as.integer(Data_year_duration),
    Exp_year_duration = as.integer(Exp_year_duration),
    env_group = factor(Env_Type, levels = env_levels),
    reference_ID = factor(reference_ID),
    site_year_ID = factor(site_year_ID),
    N_rate_cluster = as.character(N_rate_cluster)
  ) %>%
  filter(
    Exp_year_duration >= 5,
    !is.na(env_group),
    !is.na(reference_ID),
    !is.na(site_year_ID)
  )

write.csv(
  d_base %>%
    count(env_group, name = "k") %>%
    mutate(n_ref = as.integer(sapply(env_group, function(x) n_distinct(d_base$reference_ID[d_base$env_group == x]))),
           n_sy = as.integer(sapply(env_group, function(x) n_distinct(d_base$site_year_ID[d_base$env_group == x])))),
  file.path(out_dir, "long_term_counts_by_env.csv"),
  row.names = FALSE
)

write.csv(
  d_base %>%
    filter(!is.na(N_rate_cluster), N_rate_cluster != "") %>%
    mutate(N_rate_cluster = factor(N_rate_cluster, levels = unique(c(ncl_levels, sort(unique(N_rate_cluster)))))) %>%
    count(N_rate_cluster, env_group, name = "k") %>%
    group_by(N_rate_cluster, env_group) %>%
    summarise(k = sum(k), .groups = "drop") %>%
    left_join(
      d_base %>%
        filter(!is.na(N_rate_cluster), N_rate_cluster != "") %>%
        group_by(N_rate_cluster, env_group) %>%
        summarise(
          n_ref = n_distinct(reference_ID),
          n_sy = n_distinct(site_year_ID),
          .groups = "drop"
        ),
      by = c("N_rate_cluster", "env_group")
    ),
  file.path(out_dir, "long_term_counts_by_Ncluster_env.csv"),
  row.names = FALSE
)

fit_metric <- function(metric, is_lnrr) {
  metric_dir <- file.path(out_dir, metric)
  if (!dir.exists(metric_dir)) dir.create(metric_dir, recursive = TRUE)

  d <- d_base %>%
    mutate(response = as.numeric(.data[[metric]])) %>%
    filter(is.finite(response))

  count_env <- d %>%
    group_by(env_group) %>%
    summarise(
      k = n(),
      n_ref = n_distinct(reference_ID),
      n_sy = n_distinct(site_year_ID),
      raw_mean = mean(response, na.rm = TRUE),
      raw_sd = sd(response, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(count_env, file.path(metric_dir, "counts_and_raw_means_by_env.csv"), row.names = FALSE)

  m_env_obj <- fit_model(response ~ 0 + env_group, d)
  m_env <- m_env_obj$fit

  env_emm <- as.data.frame(emmeans(m_env, ~ env_group, data = d)) %>%
    rename(est = emmean, se = SE, lwr = lower.CL, upr = upper.CL) %>%
    mutate(
      metric = metric,
      effect = if (is_lnrr) (exp(est) - 1) * 100 else est,
      effect_lwr = if (is_lnrr) (exp(lwr) - 1) * 100 else lwr,
      effect_upr = if (is_lnrr) (exp(upr) - 1) * 100 else upr,
      effect_unit = if (is_lnrr) "relative yield response (%)" else "AEN",
      model_type = m_env_obj$model_type
    ) %>%
    left_join(count_env %>% select(env_group, k, n_ref, n_sy), by = "env_group")
  write.csv(env_emm, file.path(metric_dir, "env_effects_long_term_lme.csv"), row.names = FALSE)

  env_pairs <- as.data.frame(pairs(emmeans(m_env, ~ env_group, data = d), adjust = "tukey")) %>%
    mutate(metric = metric)
  write.csv(env_pairs, file.path(metric_dir, "pairwise_env_long_term_Tukey.csv"), row.names = FALSE)

  d_n <- d %>%
    filter(!is.na(N_rate_cluster), N_rate_cluster != "") %>%
    mutate(N_rate_cluster = factor(N_rate_cluster, levels = unique(c(ncl_levels, sort(unique(N_rate_cluster))))))

  count_cell <- d_n %>%
    group_by(N_rate_cluster, env_group) %>%
    summarise(
      k = n(),
      n_ref = n_distinct(reference_ID),
      n_sy = n_distinct(site_year_ID),
      raw_mean = mean(response, na.rm = TRUE),
      raw_sd = sd(response, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(count_cell, file.path(metric_dir, "counts_and_raw_means_by_Ncluster_env.csv"), row.names = FALSE)

  m_cell_obj <- fit_model(response ~ 0 + env_group:N_rate_cluster, d_n)
  m_cell <- m_cell_obj$fit

  cell_emm <- as.data.frame(emmeans(m_cell, ~ env_group | N_rate_cluster, data = d_n)) %>%
    rename(est = emmean, se = SE, lwr = lower.CL, upr = upper.CL) %>%
    mutate(
      metric = metric,
      effect = if (is_lnrr) (exp(est) - 1) * 100 else est,
      effect_lwr = if (is_lnrr) (exp(lwr) - 1) * 100 else lwr,
      effect_upr = if (is_lnrr) (exp(upr) - 1) * 100 else upr,
      effect_unit = if (is_lnrr) "relative yield response (%)" else "AEN",
      model_type = m_cell_obj$model_type
    ) %>%
    left_join(count_cell %>% select(N_rate_cluster, env_group, k, n_ref, n_sy), by = c("N_rate_cluster", "env_group")) %>%
    arrange(N_rate_cluster, env_group)
  write.csv(cell_emm, file.path(metric_dir, "env_effects_within_Ncluster_long_term_lme.csv"), row.names = FALSE)

  m_int_obj <- fit_model(response ~ env_group * N_rate_cluster, d_n)
  m_int <- m_int_obj$fit

  write.csv(as.data.frame(anova(m_int)), file.path(metric_dir, "ANOVA_env_by_Ncluster_long_term_REML.csv"))

  pairwise_n <- as.data.frame(pairs(emmeans(m_int, ~ env_group | N_rate_cluster, data = d_n), adjust = "tukey")) %>%
    mutate(metric = metric) %>%
    arrange(N_rate_cluster, contrast)
  write.csv(pairwise_n, file.path(metric_dir, "pairwise_env_within_Ncluster_long_term_Tukey.csv"), row.names = FALSE)

  writeLines(
    c(
      paste0("metric: ", metric),
      paste0("env model: ", m_env_obj$model_type),
      paste0("cell model: ", m_cell_obj$model_type),
      paste0("interaction model: ", m_int_obj$model_type),
      "fallback errors:",
      paste(c(m_env_obj$errors, m_cell_obj$errors, m_int_obj$errors), collapse = "\n")
    ),
    file.path(metric_dir, "model_fit_notes.txt")
  )

  invisible(list(env = env_emm, cell = cell_emm))
}

fit_metric("lnRR", TRUE)
fit_metric("AEN", FALSE)

cat("Done. Outputs saved to: ", normalizePath(out_dir, winslash = "/"), "\n")
