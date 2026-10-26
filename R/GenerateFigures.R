# GenerateFigures.R
# Figures of the site-level Bayesian regressions fitted in BayesianSiteLevel.Rmd.
#
#   Figure 2  Figures/Bayesian_post_slopes.png (and .svg)
#   Figure 3  Figures/Bayesian_TL_SVDC_LD_vs_Area_partial_origscale.png
#   Fig. S    Figures/Bayesian_RD_Mod_ST_vs_Area_partial.png
#   Fig. S    Figures/Bayesian_all_vs_latitude_partial.png
#   Fig. S    Figures/Bayesian_all_vs_S_partial.png
#   Fig. S    Figures/Bayesian_all_vs_impact_mean_partial.png
#   Fig. S    Figures/Bayesian_pp_checks.png
#
# Reads: Data/brms_site_level_fits.rds, Data/site_level_n7.rds,
#        Data/site_replicates_n7.rds (see section 2), Data/support_table_n7.rds,
#        Data/loo_site_results.rds (all written by BayesianSiteLevel.Rmd).
#
# Requires: tidyverse, brms, patchwork, bayesplot, posterior, ggtext

library(tidyverse)
library(brms)
library(patchwork)
library(bayesplot)
library(posterior)
library(ggtext)

dir.create("Figures", showWarnings = FALSE)

# ---- 1. Fitted models and site-level data ----
fits_site  <- readRDS("Data/brms_site_level_fits.rds")
site_level <- readRDS("Data/site_level_n7.rds")

responses <- c("C_logit", "Mod_logit", "Rank_logit", "LD_log", "TL_z", "SVD_z", "ME_log")

nice_names <- c(
  C_logit    = "Connectance (logit)",
  Mod_logit  = "Modularity (logit)",
  Rank_logit = "Rank Deficiency (logit)",
  LD_log     = "Link Density (log)",
  TL_z       = "Trophic Level (z)",
  SVD_z      = "SVD Complexity (z)",
  ME_log     = "Dynamic Stability (log)"
)

predictor_names <- c(
  S_s           = "Number of trophic species (standardized)",
  log_area_s    = "Log(area) (standardized)",
  latitude_s    = "Latitude (standardized)",
  impact_mean_s = "Human impact index (standardized)"
)

# ---- 2. Replicates used to fit the models (for the 95% bars) ----
# The same 200 metaweb-derived replicates per food web that BayesianSiteLevel.Rmd
# collapsed into site_level_n7.rds. If that file does not exist (models fitted
# before it was saved), the replicates are rebuilt; they are then a different
# random subsample from the one used to fit the models.
if (file.exists("Data/site_replicates_n7.rds")) {
  d <- readRDS("Data/site_replicates_n7.rds")
} else {
  warning("Data/site_replicates_n7.rds not found: rebuilding the replicates. ",
          "Re-run BayesianSiteLevel.Rmd so that figures and models use the same replicates.")
  source("R/site_replicates.R")
  d <- build_site_replicates()
}

site_quantiles <- d %>%
  group_by(site) %>%
  summarise(across(all_of(responses),
                   list(q025 = ~quantile(.x, 0.025, na.rm = TRUE),
                        q975 = ~quantile(.x, 0.975, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")

site_plot_data <- site_level %>% left_join(site_quantiles, by = "site")

# ---- 3. Posterior predictive checks ----
pp_plots <- lapply(responses, function(r) {
  brms::pp_check(fits_site[[r]], ndraws = 100) + ggtitle(nice_names[[r]])
})
fig_pp <- wrap_plots(pp_plots, ncol = 3)
ggsave("Figures/Bayesian_pp_checks.png", fig_pp, width = 12, height = 9, dpi = 300)

# ---- 4. Back-transformation to the original scale of each metric ----
#   C_logit, Mod_logit, Rank_logit  -> plogis()  (proportion 0-1)
#   LD_log, ME_log                  -> exp()     (original units)
#   TL_z, SVD_z                     -> z*sd+mean (original units)
# For ME_log the back-transform gives -Re(lambda_max) (positive, larger = more stable).
TL_center  <- attr(scale(d$TLmean), "scaled:center")
TL_scale   <- attr(scale(d$TLmean), "scaled:scale")
SVD_center <- attr(scale(d$SVDComplexity), "scaled:center")
SVD_scale  <- attr(scale(d$SVDComplexity), "scaled:scale")

inv_transform <- list(
  C_logit    = plogis,
  Mod_logit  = plogis,
  Rank_logit = plogis,
  LD_log     = exp,
  ME_log     = exp,
  TL_z       = function(z) z * TL_scale + TL_center,
  SVD_z      = function(z) z * SVD_scale + SVD_center
)

raw_nice_names <- c(
  C_logit    = "Connectance",
  Mod_logit  = "Modularity",
  Rank_logit = "Rank Deficiency",
  LD_log     = "Link Density",
  TL_z       = "Mean Trophic Level",
  SVD_z      = "SVD Complexity",
  ME_log     = "Dynamic Stability (-Re(λmax))"
)

# ---- 5. Partial residuals in the original scale (Figure 3) ----
# Partial residual of each food web for predictor `pred`: observed value (mean of
# the replicates, transformed scale) minus the contribution of the other three
# predictors, using the posterior mean of their slopes. Points are then on the
# same scale as the marginal-effect curve (other predictors = 0). Bars (2.5 and
# 97.5% quantiles of the replicates) are shifted by the same amount. Everything is
# back-transformed with inv_transform. For area the x axis is labelled in km2
# (log scale), while the position is the standardized log(area) used by the model.
partial_plot_orig <- function(resp, pred) {
  fit <- fits_site[[resp]]
  all_preds <- c("S_s", "log_area_s", "latitude_s", "impact_mean_s")
  other_preds <- setdiff(all_preds, pred)

  newdata <- tibble(!!pred := seq(min(site_level[[pred]]), max(site_level[[pred]]), length.out = 100))
  for (op in other_preds) newdata[[op]] <- 0
  newdata[[paste0(resp, "_sd")]] <- 0  # mean function, without the measurement-error term
  f <- inv_transform[[resp]]
  pred_fit <- fitted(fit, newdata = newdata, re_formula = NA) %>% as_tibble() %>%
    mutate(across(c(Estimate, Q2.5, Q97.5), f))
  newdata <- bind_cols(newdata, pred_fit)

  b <- fixef(fit)[, "Estimate"]
  adj <- as.vector(as.matrix(site_plot_data[, other_preds]) %*% b[other_preds])
  pr <- site_plot_data %>%
    transmute(site, x = .data[[pred]],
              y   = f(.data[[paste0(resp, "_mean")]] - adj),
              ylo = f(.data[[paste0(resp, "_q025")]] - adj),
              yhi = f(.data[[paste0(resp, "_q975")]] - adj))

  p <- ggplot() +
    geom_ribbon(data = newdata, aes(x = .data[[pred]], ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
    geom_line(data = newdata, aes(x = .data[[pred]], y = Estimate), linewidth = 0.8) +
    geom_errorbar(data = pr, aes(x = x, ymin = ylo, ymax = yhi, color = site), width = 0, alpha = 0.6) +
    geom_point(data = pr, aes(x = x, y = y, color = site), size = 2.2) +
    labs(x = predictor_names[[pred]], y = raw_nice_names[[resp]], color = "Food web") +
    theme_bw(base_size = 11)

  if (pred == "log_area_s") {
    # The standardization was done over the replicates, not over the 7 sites,
    # so the exact relation is taken from the data.
    std <- coef(lm(log_area_s ~ log_area, data = site_level))
    km2 <- 10^(1:6)
    p <- p + scale_x_continuous(name = expression(Area~(km^2)),
                                breaks = std[1] + std[2] * log(km2),
                                labels = parse(text = paste0("10^", 1:6)))
  }
  p
}

fig3_partial <- (partial_plot_orig("TL_z", "log_area_s") + partial_plot_orig("SVD_z", "log_area_s") +
                 partial_plot_orig("C_logit", "log_area_s") + partial_plot_orig("LD_log", "log_area_s")) +
  plot_layout(ncol = 2, guides = "collect")
ggsave("Figures/Bayesian_TL_SVDC_LD_vs_Area_partial_origscale.png", fig3_partial, width = 9, height = 7, dpi = 300)

# ---- 6. Supplementary figures of partial residuals for the four predictors ----
# Same calculation as partial_plot_orig(), with the x axis labelled in the
# original units of each predictor.
raw_pred <- c(S_s = "S", log_area_s = "log_area", latitude_s = "latitude", impact_mean_s = "impact_mean")
raw_pred_names <- list(
  S_s           = "Number of trophic species",
  log_area_s    = expression(Area~(km^2)),
  latitude_s    = "Latitude (°S)",
  impact_mean_s = "Human impact index"
)

x_scale_orig <- function(pred) {
  std <- coef(lm(as.formula(paste(pred, "~", raw_pred[[pred]])), data = site_level))
  if (pred == "log_area_s") {
    raw_breaks <- log(10^(1:6)); labs_x <- parse(text = paste0("10^", 1:6))
  } else {
    raw_breaks <- pretty(site_level[[raw_pred[[pred]]]], n = 5); labs_x <- format(raw_breaks)
  }
  scale_x_continuous(name = raw_pred_names[[pred]], breaks = std[1] + std[2] * raw_breaks, labels = labs_x)
}

partial_plot_supp <- function(resp, pred) {
  suppressMessages(partial_plot_orig(resp, pred) + x_scale_orig(pred))
}

supp_resp <- c("C_logit", "LD_log", "TL_z", "SVD_z", "Rank_logit", "Mod_logit", "ME_log")

fig_area_supp <- wrap_plots(lapply(c("Rank_logit", "Mod_logit", "ME_log"), partial_plot_supp, pred = "log_area_s"), ncol = 3) +
  plot_layout(guides = "collect")
ggsave("Figures/Bayesian_RD_Mod_ST_vs_Area_partial.png", fig_area_supp, width = 11, height = 3.8, dpi = 300)

for (pp in c("latitude_s", "S_s", "impact_mean_s")) {
  fig_pp <- wrap_plots(lapply(supp_resp, partial_plot_supp, pred = pp), ncol = 3) +
    plot_layout(guides = "collect")
  ggsave(paste0("Figures/Bayesian_all_vs_", sub("_s$", "", pp), "_partial.png"), fig_pp, width = 11, height = 10, dpi = 300)
}

# ---- 7. Figure 2: posterior slopes, one panel per predictor ----
# Coloured bands mark effects with posterior probability of a consistent sign
# >= 0.90 that are robust to leave-one-site-out (fragile in <= 3 of 7 refits).
support_table    <- readRDS("Data/support_table_n7.rds")
loo_site_results <- readRDS("Data/loo_site_results.rds")

site_nice_names <- c(
  S_s = "Trophic Species", log_area_s = "Area (log)",
  latitude_s = "Latitude", impact_mean_s = "Human Impact",
  C_logit = "Connectance", Mod_logit = "Modularity",
  Rank_logit = "Rank Deficiency", LD_log = "Link Density",
  TL_z = "Trophic Level", SVD_z = "SVD Complexity", ME_log = "Dyn. Stability"
)
rev_nice_names <- setNames(names(site_nice_names), site_nice_names)

# Same calculation as fragile_effects in section 5 of BayesianSiteLevel.Rmd
fragile_effects <- support_table %>%
  filter(prob_max >= 0.90) %>%
  select(response, predictor, direction_full = direction) %>%
  left_join(loo_site_results, by = c("response", "predictor")) %>%
  mutate(fragile = prob_max < 0.90 | direction != direction_full) %>%
  group_by(response, predictor, direction_full) %>%
  summarise(n_fragile_of_7 = sum(fragile), .groups = "drop")

robust_table <- support_table %>%
  left_join(fragile_effects %>% select(response, predictor, n_fragile_of_7),
            by = c("response", "predictor")) %>%
  mutate(
    robust   = !is.na(n_fragile_of_7) & prob_max >= 0.90 & n_fragile_of_7 <= 3,
    resp_raw = rev_nice_names[response],
    pred_raw = rev_nice_names[predictor],
    param    = paste0(resp_raw, "__", pred_raw)
  )

# Posterior draws of the slopes of the 7 fits, columns "<resp>__<predictor>"
site_preds <- c("S_s", "log_area_s", "latitude_s", "impact_mean_s")
draws_wide <- purrr::imap(fits_site, function(fit, resp) {
  dd <- as_draws_df(fit) %>% dplyr::select(all_of(paste0("b_", site_preds)))
  names(dd) <- paste0(resp, "__", site_preds)
  dd
}) %>% dplyr::bind_cols()

make_band <- function(rows, fill, alpha) {
  if (nrow(rows) == 0) return(NULL)
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = rows$ypos - 0.5, ymax = rows$ypos + 0.5,
           fill = fill, alpha = alpha)
}

metric_order <- c("C_logit", "LD_log", "TL_z", "SVD_z", "Rank_logit", "Mod_logit", "ME_log")
pred_order   <- c("log_area_s", "latitude_s", "S_s", "impact_mean_s")
pred_titles  <- c(log_area_s = "Area (log)", latitude_s = "Latitude",
                  S_s = "Number of trophic species", impact_mean_s = "Human impact")
x_lim <- range(sapply(draws_wide, quantile, probs = c(0.005, 0.995)))

# prob = 0.8: the shaded central interval excludes 0 exactly when the posterior
# probability of a consistent sign exceeds 0.90. prob_outer = 0.99 shows almost
# the whole posterior.
post_panel <- function(pred, show_y = TRUE) {
  pars <- paste0(metric_order, "__", pred)
  lims <- rev(pars)
  rt <- robust_table %>% filter(param %in% pars)
  lab <- setNames(site_nice_names[metric_order], pars)
  bg <- tibble(param = lims, ypos = seq_along(lims)) %>%
    left_join(rt %>% select(param, robust, direction), by = "param") %>%
    mutate(band = dplyr::case_when(robust %in% TRUE & direction == "positive" ~ "positive",
                                   robust %in% TRUE & direction == "negative" ~ "negative",
                                   TRUE ~ "neutral"))
  p <- mcmc_areas(draws_wide, pars = pars, prob = 0.8, prob_outer = 0.99)
  bands <- list(
    make_band(dplyr::filter(bg, band == "positive"), "#2166ac", 0.6),
    make_band(dplyr::filter(bg, band == "negative"), "#d6604d", 0.6),
    make_band(dplyr::filter(bg, band == "neutral"),  "#f5f5f5", 0.6)
  ) %>% purrr::compact()
  p$layers <- c(bands, p$layers)
  p + scale_y_discrete(limits = lims, labels = lab) +
    coord_cartesian(xlim = x_lim) +
    geom_vline(xintercept = 0, linewidth = 0.8) +
    labs(x = pred_titles[[pred]], y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          axis.text.y = if (show_y) element_text() else element_blank(),
          axis.ticks.y = if (show_y) element_line() else element_blank())
}

# Metric names only in the left column (same order in all panels)
fig_post_slopes_4 <- wrap_plots(mapply(post_panel, pred_order, c(TRUE, FALSE, TRUE, FALSE),
                                       SIMPLIFY = FALSE), ncol = 2)
ggsave("Figures/Bayesian_post_slopes.png", fig_post_slopes_4, width = 7, height = 6, dpi = 300, bg = "white")
ggsave("Figures/Bayesian_post_slopes.svg", fig_post_slopes_4, width = 7, height = 6, bg = "white")
