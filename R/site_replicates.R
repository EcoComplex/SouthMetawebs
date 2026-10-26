# build_site_replicates()
# Builds the replicate-level data used by the site-level Bayesian regressions:
# for each of the seven food webs, `n` metaweb-derived replicates with the six
# structural metrics, paired with one dynamic-stability replicate, the four
# predictors, and the transformed responses.
#
# Used by BayesianSiteLevel.Rmd (which saves the result to
# Data/site_replicates_n7.rds) and, as a fallback, by R/GenerateFigures.R.

library(tidyverse)

analysis_sites <- c("BurdwoodBank", "GulfSanJorge", "PotterCove",
                    "Std BeagleChannel", "Std Weddell Sea",
                    "Std NorthernScotia", "Std SouthernScotia")

build_site_replicates <- function(n = 200, seed = 2025) {
  network_info <- readRDS("Data/network_info.rds")
  stability_results_meig <- readRDS("Data/stability_results_meig.rds")
  simMetaWebMetrics <- readRDS("Data/simulations_metaWeb_metrics.rds") %>%
    mutate(SVDComplexity = 1 - Entropy, RankDeficiency = 1 - relRank) %>%
    filter(Metaweb %in% analysis_sites)

  stability_id <- stability_results_meig %>%
    select(site, MEing_stable) %>%
    group_by(site) %>%
    mutate(Sim = row_number()) %>%
    ungroup()

  set.seed(seed)
  sim_final <- simMetaWebMetrics %>%
    select(-S) %>%
    rename(site = Metaweb) %>%
    inner_join(stability_id, by = c("site", "Sim")) %>%
    group_by(site) %>%
    slice_sample(n = n) %>%
    ungroup()

  env <- network_info %>%
    select(site, latitude, area_km2, depth_m, impact_mean, S) %>%
    mutate(log_area = log(area_km2))

  d <- sim_final %>%
    left_join(env, by = "site") %>%
    drop_na() %>%
    mutate(
      C_logit    = qlogis(pmin(pmax(C, 1e-6), 1 - 1e-6)),
      Mod_logit  = qlogis(pmin(pmax(Modularity, 1e-6), 1 - 1e-6)),
      Rank_logit = qlogis(pmin(pmax(RankDeficiency, 1e-6), 1 - 1e-6)),
      LD_log     = log(LD),
      TL_z       = as.numeric(scale(TLmean)),
      SVD_z      = as.numeric(scale(SVDComplexity)),
      MEpos      = -MEing_stable,
      ME_log     = log(MEpos),
      S_s           = as.numeric(scale(S)),
      log_area_s    = as.numeric(scale(log_area)),
      latitude_s    = as.numeric(scale(latitude)),
      impact_mean_s = as.numeric(scale(impact_mean))
    )

  stopifnot(n_distinct(d$site) == 7)
  d
}
