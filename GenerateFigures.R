# GenerateFigures.R
# Regenera las figuras Bayesian_*.png a partir de los modelos site-level
# ya ajustados en BayesianSiteLevel.Rmd (Data/brms_site_level_fits.rds).
# No corrido por mi (sin acceso a R/cmdstan en este entorno) -- probable que
# necesites ajustar algun detalle (theme, tamanios, nombres de columnas si
# cambiaste algo en BayesianSiteLevel.Rmd).
#
# Requiere: tidyverse, brms, patchwork (install.packages("patchwork") si falta)

library(tidyverse)
library(brms)
library(patchwork)

dir.create("Figures", showWarnings = FALSE)

# ---- 1. Cargar fits y datos site-level ----
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

# ---- 2. Reconstruir `d` (replicas individuales) para las barras de 95% ----
# Identico al chunk load-and-transform de BayesianSiteLevel.Rmd.
network_info <- readRDS("Data/network_info.rds")
stability_results_meig <- readRDS("Data/stability_results_meig.rds")
simMetaWebMetrics <- readRDS("Data/simulations_metaWeb_metrics.rds") %>%
  mutate(SVDComplexity = 1 - Entropy, RankDeficiency = 1 - relRank) %>%
  filter(!(Metaweb %in% c("WeddellSea", "BeagleChannel", "NorthernScotia", "SouthernScotia")))

stability_id <- stability_results_meig %>%
  select(site, MEing_stable) %>%
  group_by(site) %>%
  mutate(Sim = row_number()) %>%
  ungroup()

sim_final <- simMetaWebMetrics %>%
  select(-S) %>%
  rename(site = Metaweb) %>%
  inner_join(stability_id, by = c("site", "Sim")) %>%
  group_by(site) %>%
  slice_sample(n = 200) %>%
  ungroup()

env <- network_info %>%
  select(site, latitude, area_km2, depth_m, impact_mean, S) %>%
  mutate(log_area = log(area_km2))

df_all <- sim_final %>%
  left_join(env, by = "site") %>%
  drop_na()

d <- df_all %>%
  mutate(
    C_logit    = qlogis(pmin(pmax(C, 1e-6), 1 - 1e-6)),
    Mod_logit  = qlogis(pmin(pmax(Modularity, 1e-6), 1 - 1e-6)),
    Rank_logit = qlogis(pmin(pmax(RankDeficiency, 1e-6), 1 - 1e-6)),
    LD_log     = log(LD),
    TL_z       = as.numeric(scale(TLmean)),
    SVD_z      = as.numeric(scale(SVDComplexity)),
    MEpos      = -MEing_stable,
    ME_log     = log(MEpos)
  )

site_quantiles <- d %>%
  group_by(site) %>%
  summarise(across(all_of(responses),
                    list(q025 = ~quantile(.x, 0.025, na.rm = TRUE),
                         q975 = ~quantile(.x, 0.975, na.rm = TRUE)),
                    .names = "{.col}_{.fn}"),
            .groups = "drop")

site_plot_data <- site_level %>% left_join(site_quantiles, by = "site")

# ---- 3. Funcion generica: efecto marginal de un predictor sobre una metrica ----
# Grafica en la escala transformada del modelo (misma escala en la que se
# ajusto se(); evita un paso extra de back-transformacion por metrica).
marginal_plot <- function(resp, pred) {
  fit <- fits_site[[resp]]
  all_preds <- c("S_s", "log_area_s", "latitude_s", "impact_mean_s")
  other_preds <- setdiff(all_preds, pred)

  newdata <- tibble(!!pred := seq(min(site_level[[pred]]), max(site_level[[pred]]), length.out = 100))
  for (op in other_preds) newdata[[op]] <- 0

  mean_var <- paste0(resp, "_mean")
  sd_var   <- paste0(resp, "_sd")
  newdata[[sd_var]] <- 0  # predecir la funcion media, sin sumar el termino de error de medicion

  pred_fit <- fitted(fit, newdata = newdata, re_formula = NA) %>% as_tibble()
  newdata <- bind_cols(newdata, pred_fit)

  q025_col <- paste0(resp, "_q025")
  q975_col <- paste0(resp, "_q975")

  ggplot() +
    geom_ribbon(data = newdata, aes(x = .data[[pred]], ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
    geom_line(data = newdata, aes(x = .data[[pred]], y = Estimate), linewidth = 0.8) +
    geom_errorbar(data = site_plot_data,
                  aes(x = .data[[pred]], ymin = .data[[q025_col]], ymax = .data[[q975_col]], color = site),
                  width = 0, alpha = 0.6) +
    geom_point(data = site_plot_data,
               aes(x = .data[[pred]], y = .data[[mean_var]], color = site), size = 2.2) +
    labs(x = predictor_names[[pred]], y = nice_names[[resp]], color = "Food web") +
    theme_bw(base_size = 11)
}

# ---- 4. Figura 2 del texto principal: efectos robustos de Area ----
fig2 <- (marginal_plot("TL_z", "log_area_s") + marginal_plot("SVD_z", "log_area_s") +
         marginal_plot("C_logit", "log_area_s") + marginal_plot("LD_log", "log_area_s")) +
  plot_layout(ncol = 2, guides = "collect")
ggsave("Figures/Bayesian_TL_SVDC_LD_vs_Area_mv.png", fig2, width = 9, height = 7, dpi = 300)

# ---- 5. SuppMat: Area sobre C, Rank Deficiency, Modularity, Stability ----
fig_area_supp <- (marginal_plot("C_logit", "log_area_s") + marginal_plot("Rank_logit", "log_area_s") +
                   marginal_plot("Mod_logit", "log_area_s") + marginal_plot("ME_log", "log_area_s")) +
  plot_layout(ncol = 2, guides = "collect")
ggsave("Figures/Bayesian_C_RD_Mod_ST_vs_Area_mv.png", fig_area_supp, width = 9, height = 7, dpi = 300)

# ---- 6. SuppMat: Latitud sobre C, Rank Deficiency, Modularity, Stability ----
fig_lat_supp <- (marginal_plot("C_logit", "latitude_s") + marginal_plot("Rank_logit", "latitude_s") +
                  marginal_plot("Mod_logit", "latitude_s") + marginal_plot("ME_log", "latitude_s")) +
  plot_layout(ncol = 2, guides = "collect")
ggsave("Figures/Bayesian_C_RD_Mod_ST_vs_Latitude_mv.png", fig_lat_supp, width = 9, height = 7, dpi = 300)

# ---- 7. SuppMat: numero de especies (S) sobre las 7 metricas ----
fig_S <- wrap_plots(lapply(responses, marginal_plot, pred = "S_s"), ncol = 3) +
  plot_layout(guides = "collect")
ggsave("Figures/Bayesian_all_vs_S_mv.png", fig_S, width = 12, height = 9, dpi = 300)

# ---- 8. SuppMat: impacto humano sobre las 7 metricas ----
fig_impact <- wrap_plots(lapply(responses, marginal_plot, pred = "impact_mean_s"), ncol = 3) +
  plot_layout(guides = "collect")
ggsave("Figures/Bayesian_all_vs_impact_mv.png", fig_impact, width = 12, height = 9, dpi = 300)

# ---- 9. Posterior predictive checks (reemplaza al pp_checks del modelo viejo) ----
pp_plots <- lapply(responses, function(r) {
  brms::pp_check(fits_site[[r]], ndraws = 100) + ggtitle(nice_names[[r]])
})
fig_pp <- wrap_plots(pp_plots, ncol = 3)
ggsave("Figures/Bayesian_pp_checks.png", fig_pp, width = 12, height = 9, dpi = 300)


# ---- 10. Figura S6 (SuppMat): correlaciones entre covariables ambientales ----
# Esta figura NO se regenera aca. Se genera en metawebassemblysimulations.Rmd,
# chunk "CheckCorrelation_createTables" (eval=FALSE, hay que correrlo a mano
# en RStudio): ese chunk ya hace el join que agrega temp_C via la columna
# `name`, guarda network_info.rds actualizado, y escribe
# Figures/correlation_network_info.png con el estilo de ggpairs original
# (cor en el panel superior, puntos abajo, densidad en la diagonal). Se le
# agrego S al select() ahi para cumplir el pedido del Revisor 1. No lo
# duplico aca para no tener dos versiones de la misma logica divergiendo.

cat("Listo. Figuras escritas en Figures/. Revisar visualmente antes de compilar.\n")

# ---- 11. Distribucion cruda por red de cada metrica (sin modelo) ----
# Complementa la tabla de metricas: dado que solo el area dio efectos
# robustos, esto muestra la variacion cruda (replicas de metaweb) red por
# red para las 7 metricas, en su escala original (no logit/log/z), para
# poder discutir cualitativamente si alguna red se separa visiblemente del
# resto. Es puramente descriptivo, sin ajuste de modelo ni test estadistico
# -- ojo con no convertir esto en una interpretacion post-hoc de "por que
# esta red es distinta" sin dejar claro que es exploratorio.
raw_vars <- c("C", "Modularity", "RankDeficiency", "LD", "TLmean", "SVDComplexity", "MEing_stable")

raw_nice_names <- c(
  C              = "Connectance",
  Modularity     = "Modularity",
  RankDeficiency = "Rank Deficiency",
  LD             = "Link Density",
  TLmean         = "Mean Trophic Level",
  SVDComplexity  = "SVD Complexity",
  MEing_stable   = "Dynamic Stability (Re(λmax))"
)

site_order <- site_level %>% arrange(latitude) %>% pull(site)
d_dist <- d %>% mutate(site = factor(site, levels = site_order))

distribution_plot <- function(var) {
  ggplot(d_dist, aes(x = site, y = .data[[var]], fill = site)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.7) +
    labs(x = NULL, y = raw_nice_names[[var]]) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

fig_dist <- wrap_plots(lapply(raw_vars, distribution_plot), ncol = 3)
ggsave("Figures/metric_distributions_by_site.png", fig_dist, width = 12, height = 9, dpi = 300, bg = "white")

# ---- 12. Efectos marginales en escala original (para comparar con las transformadas) ----
# Mismas figuras que los bloques 4 a 8, pero con la curva ajustada
# destransformada de vuelta a la escala natural de cada metrica:
#   C_logit, Mod_logit, Rank_logit  -> plogis()  (proporcion 0-1)
#   LD_log, ME_log                  -> exp()     (unidades originales)
#   TL_z, SVD_z                     -> z*sd+mean (unidades originales)
# Para ME_log el back-transform da MEpos = -Re(lambda_max) (positivo, mayor
# = mas estable), mismo signo que se usa en el resto del script.
#
# Decision a tener en cuenta: los puntos/barras observados por sitio NO son
# un back-transform de site_plot_data (que ya esta en escala transformada);
# se recalculan directo desde `d` en escala original. Promediar en escala
# logit/log y despues destransformar no da lo mismo que promediar directo
# en escala original (desigualdad de Jensen), y para "comparar visualmente"
# la escala natural es la que un lector va a interpretar, asi que use el
# promedio calculado ahi. Si preferis que el punto observado sea el
# back-transform exacto del promedio ya usado en las figuras actuales,
# avisame y lo cambio.

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

raw_col <- c(
  C_logit    = "C",
  Mod_logit  = "Modularity",
  Rank_logit = "RankDeficiency",
  LD_log     = "LD",
  TL_z       = "TLmean",
  SVD_z      = "SVDComplexity",
  ME_log     = "MEpos"
)

raw_nice_names2 <- c(
  C_logit    = "Connectance",
  Mod_logit  = "Modularity",
  Rank_logit = "Rank Deficiency",
  LD_log     = "Link Density",
  TL_z       = "Mean Trophic Level",
  SVD_z      = "SVD Complexity",
  ME_log     = "Dynamic Stability (-Re(λmax))"
)

site_raw_summary <- d %>%
  group_by(site) %>%
  summarise(across(all_of(unname(raw_col)),
                    list(mean = ~mean(.x, na.rm = TRUE),
                         q025 = ~quantile(.x, 0.025, na.rm = TRUE),
                         q975 = ~quantile(.x, 0.975, na.rm = TRUE)),
                    .names = "{.col}_{.fn}"),
            .groups = "drop")

site_plot_data_raw <- site_level %>% left_join(site_raw_summary, by = "site")

marginal_plot_orig <- function(resp, pred) {
  fit <- fits_site[[resp]]
  all_preds <- c("S_s", "log_area_s", "latitude_s", "impact_mean_s")
  other_preds <- setdiff(all_preds, pred)

  newdata <- tibble(!!pred := seq(min(site_level[[pred]]), max(site_level[[pred]]), length.out = 100))
  for (op in other_preds) newdata[[op]] <- 0

  sd_var <- paste0(resp, "_sd")
  newdata[[sd_var]] <- 0

  pred_fit <- fitted(fit, newdata = newdata, re_formula = NA) %>% as_tibble()
  f <- inv_transform[[resp]]
  pred_fit <- pred_fit %>% mutate(across(c(Estimate, Q2.5, Q97.5), f))
  newdata <- bind_cols(newdata, pred_fit)

  rc <- raw_col[[resp]]
  mean_col <- paste0(rc, "_mean")
  q025_col <- paste0(rc, "_q025")
  q975_col <- paste0(rc, "_q975")

  ggplot() +
    geom_ribbon(data = newdata, aes(x = .data[[pred]], ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
    geom_line(data = newdata, aes(x = .data[[pred]], y = Estimate), linewidth = 0.8) +
    geom_errorbar(data = site_plot_data_raw,
                  aes(x = .data[[pred]], ymin = .data[[q025_col]], ymax = .data[[q975_col]], color = site),
                  width = 0, alpha = 0.6) +
    geom_point(data = site_plot_data_raw,
               aes(x = .data[[pred]], y = .data[[mean_col]], color = site), size = 2.2) +
    labs(x = predictor_names[[pred]], y = raw_nice_names2[[resp]], color = "Food web") +
    theme_bw(base_size = 11)
}

fig2_orig <- (marginal_plot_orig("TL_z", "log_area_s") + marginal_plot_orig("SVD_z", "log_area_s") +
              marginal_plot_orig("C_logit", "log_area_s") + marginal_plot_orig("LD_log", "log_area_s")) +
  plot_layout(ncol = 2, guides = "collect")
ggsave("Figures/Bayesian_TL_SVDC_LD_vs_Area_mv_origscale.png", fig2_orig, width = 9, height = 7, dpi = 300)

fig_area_supp_orig <- (marginal_plot_orig("C_logit", "log_area_s") + marginal_plot_orig("Rank_logit", "log_area_s") +
                        marginal_plot_orig("Mod_logit", "log_area_s") + marginal_plot_orig("ME_log", "log_area_s")) +
  plot_layout(ncol = 2, guides = "collect")
ggsave("Figures/Bayesian_C_RD_Mod_ST_vs_Area_mv_origscale.png", fig_area_supp_orig, width = 9, height = 7, dpi = 300)

fig_lat_supp_orig <- (marginal_plot_orig("C_logit", "latitude_s") + marginal_plot_orig("Rank_logit", "latitude_s") +
                       marginal_plot_orig("Mod_logit", "latitude_s") + marginal_plot_orig("ME_log", "latitude_s")) +
  plot_layout(ncol = 2, guides = "collect")
ggsave("Figures/Bayesian_C_RD_Mod_ST_vs_Latitude_mv_origscale.png", fig_lat_supp_orig, width = 9, height = 7, dpi = 300)

fig_S_orig <- wrap_plots(lapply(responses, marginal_plot_orig, pred = "S_s"), ncol = 3) +
  plot_layout(guides = "collect")
ggsave("Figures/Bayesian_all_vs_S_mv_origscale.png", fig_S_orig, width = 12, height = 9, dpi = 300)

fig_impact_orig <- wrap_plots(lapply(responses, marginal_plot_orig, pred = "impact_mean_s"), ncol = 3) +
  plot_layout(guides = "collect")
ggsave("Figures/Bayesian_all_vs_impact_mv_origscale.png", fig_impact_orig, width = 12, height = 9, dpi = 300)

cat("Listo. Figuras en escala original: Figures/*_origscale.png\n")

# ---- 13. Bayesian_post_slopes.png con el criterio nuevo (soporte + LOSO) ----
# Recrea la figura de bandas de fondo (Figures/Bayesian_post_slopes.png,
# celeste/rosa/gris) que antes generaba plot_mcmc_areas_ordered() en
# R/network_fun.r a partir de fit_mv (modelo multivariado viejo). Esa
# funcion asume UN fit conjunto con parametros b_<resp>_<pred>; con 7
# fits_site separados no aplica directo, asi que este bloque reimplementa
# la misma logica (mismos colores, mismo estilo) combinando los 7 modelos,
# y cambia el criterio de "banda de color" de solo prob_max >= 0.90 (como
# estaba) a prob_max >= 0.90 Y robusto a LOSO (<=3 de 7 exclusiones
# fragiles) -- el mismo criterio que ya se usa en el resto del SuppMat.
# NO LO PROBE (sin R aca): antes de usarla para el manuscrito, compara
# contra la figura vieja para chequear que el orden, los colores y las
# etiquetas tengan sentido.
library(bayesplot)
library(posterior)
if (!requireNamespace("ggtext", quietly = TRUE)) install.packages("ggtext")
library(ggtext)

if (!file.exists("Data/support_table_n7.rds") || !file.exists("Data/loo_site_results.rds")) {
  stop("Faltan Data/support_table_n7.rds o Data/loo_site_results.rds -- corre primero las secciones 4 y 5 de BayesianSiteLevel.Rmd.")
}

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

# Mismo calculo que la seccion 5 de BayesianSiteLevel.Rmd (fragile_effects),
# repetido aca para no depender de correr ese Rmd entero, solo de sus RDS.
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

# Draws combinados de los 7 fits, columnas unicas "<resp>__<predictor>"
site_preds <- c("S_s", "log_area_s", "latitude_s", "impact_mean_s")
draws_wide <- purrr::imap(fits_site, function(fit, resp) {
  dd <- as_draws_df(fit) %>% dplyr::select(all_of(paste0("b_", site_preds)))
  names(dd) <- paste0(resp, "__", site_preds)
  dd
}) %>% dplyr::bind_cols()

ord_site <- tibble(param = names(draws_wide)) %>%
  mutate(
    mean      = sapply(draws_wide, mean),
    resp_raw  = sub("__.*$", "", param),
    pred_raw  = sub("^.*__", "", param),
    response  = site_nice_names[resp_raw],
    predictor = site_nice_names[pred_raw]
  ) %>%
  left_join(robust_table %>% select(param, robust), by = "param") %>%
  mutate(
    robust = ifelse(is.na(robust), FALSE, robust),
    # Negrita (markdown, via ggtext) para los efectos que cumplen el
    # criterio nuevo (prob_max >= 0.90 y robusto a LOSO), reemplaza el
    # retoque manual que se le hacia al SVG en la version vieja.
    label  = ifelse(robust,
                     paste0("**", predictor, " → ", response, "**"),
                     paste(predictor, "→", response))
  ) %>%
  arrange(desc(mean))

limits_vec_site <- rev(ord_site$param)

p_post_slopes_site <- mcmc_areas(draws_wide, pars = ord_site$param, prob_outer = 0.9)

pos_df <- tibble(param = limits_vec_site, ypos = seq_along(limits_vec_site))

bg_df <- pos_df %>%
  left_join(robust_table %>% select(param, robust, direction), by = "param") %>%
  mutate(
    band = dplyr::case_when(
      !is.na(robust) & robust & direction == "positive" ~ "positive",
      !is.na(robust) & robust & direction == "negative" ~ "negative",
      TRUE ~ "neutral"
    )
  )

make_band <- function(rows, fill, alpha) {
  if (nrow(rows) == 0) return(NULL)
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = rows$ypos - 0.5, ymax = rows$ypos + 0.5,
           fill = fill, alpha = alpha)
}

bg_layers_site <- list(
  make_band(dplyr::filter(bg_df, band == "positive"), "#2166ac", 0.6),
  make_band(dplyr::filter(bg_df, band == "negative"), "#d6604d", 0.6),
  make_band(dplyr::filter(bg_df, band == "neutral"),  "#f5f5f5", 0.6)
) %>% purrr::compact()

p_post_slopes_site$layers <- c(bg_layers_site, p_post_slopes_site$layers)

fig_post_slopes <- p_post_slopes_site +
  scale_y_discrete(limits = limits_vec_site, labels = setNames(ord_site$label, ord_site$param)) +
  theme_bw(base_size = 13) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = ggtext::element_markdown()) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 1) +
  labs(x = "Posterior distribution", y = NULL)

ggsave("Figures/Bayesian_post_slopes.png", fig_post_slopes, width = 5, height = 10, dpi = 300, bg = "white")
ggsave("Figures/Bayesian_post_slopes.svg", fig_post_slopes, width = 5, height = 10, bg = "white")

# ---- 14. Figura de Modularidad por red (NO USADA en SuppMat.md) ----
# Boxplot por sitio, aislado del panel de 7 metricas del bloque 11. Lo que
# Leonardo pidio documentar en la seccion de Modularity era en realidad la
# comparacion Infomap vs spinglass (CompareModularityAlgorithms.R), no esto.
# Dejo el bloque por si sirve para otra cosa, pero no esta referenciado en
# SuppMat.md.
fig_modularity_dist <- distribution_plot("Modularity")
ggsave("Figures/Modularity_distribution_by_site.png", fig_modularity_dist, width = 6, height = 5, dpi = 300, bg = "white")
