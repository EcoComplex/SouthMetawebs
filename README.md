# Ecosystem area and the structure of marine food webs along a latitudinal gradient in the Southwest Atlantic and Antarctica

[![DOI](https://zenodo.org/badge/1009246443.svg)](https://doi.org/10.5281/zenodo.21045337)

Code and data for the manuscript.

## Abstract

Understanding how food web structure varies across broad spatial gradients is essential for identifying general principles of ecosystem organization and resilience under global change. We analyzed seven highly resolved marine food webs spanning a latitudinal gradient from 48°S to 78°S in the Southwest Atlantic and Antarctic region to evaluate how ecosystem area, latitude, human impact, and number of trophic species are associated with network structure and dynamic stability. We quantified connectance, link density, mean trophic level, modularity, SVD complexity, rank deficiency, and dynamic stability, and related each metric to the four predictors using site-level Bayesian regressions that incorporate the uncertainty derived from a metaweb assembly model as a measurement-error term. Because the dataset comprises only seven independent food webs, we assessed the robustness of each relationship by refitting the models with each food web excluded in turn.

Ecosystem area showed the most consistent associations: larger ecosystems had higher mean trophic level, link density, connectance, and SVD complexity, and these associations held in all seven exclusions. The number of trophic species was positively associated with rank deficiency and modularity and negatively with connectance, but these associations lost support when a single food web, Burdwood Bank, was excluded. Latitude was positively associated with link density and SVD complexity, although these associations depended on three of the seven food webs, and its associations with connectance and trophic level did not survive exclusion. Human impact showed no clear association with any metric, possibly because the cumulative impact index integrates stressors whose effects on food web structure differ in magnitude and direction.

Among the predictors considered, ecosystem area showed the most consistent associations with marine food web structure in the Southwest Atlantic and Antarctic. Patterns associated with species richness and latitude should be treated as exploratory and require evaluation with additional independent food webs.

## Workflow

Scripts are run from the repository root (also the `.R` files in `R/`), in this order. Each step reads the `.rds` files written by the previous ones in `Data/`. Chunks that run long simulations, or that overwrite `Data/network_info.rds`, are marked `eval = FALSE` and are run by hand; the other chunks only read saved results.

| Step | Script | What it does | Writes |
|---|---|---|---|
| 0 | `R/compare_networks_refined.R` | Compares the original and standardized food webs of Beagle Channel, Weddell Sea and the Scotia Sea (diagnostic, not needed by the next steps). | `Data/standardization_comparison/` |
| 1 | `LoadNetworksDescription.Rmd` | Loads the seven food webs and computes their empirical metrics. | `network_info.rds` |
| 2 | `HalpernHumanImpact.Rmd` | Mean human impact index of each area from `Data/ssp245_current.tif`. | adds `impact_mean` to `network_info.rds` |
| 3 | `Stabilitysimulations.Rmd` | Self-damping at QSS = 0.05 and maximum eigenvalues of random community matrices (dynamic stability). | `stability_results.rds`, `stability_results_meig.rds`, adds `selfDamping` to `network_info.rds` |
| 4 | `metawebassemblysimulations.Rmd` | Metaweb assembly model (`meweasmo`): simulations, ABC fit, 1,000 simulated food webs per site and their metrics; Figures S3 to S6 and Table S1. | `simulations_metaWebAssembly_all_networks.rds`, `fittedMetaWeb_all_networks.rds`, `gof_all_networks.rds`, `simulations_metaWeb_metrics.rds`, `simulations_metaWeb_metrics_spinglass.rds` |
| 5 | `BayesianSiteLevel.Rmd` | Site-level Bayesian regressions (`brms`, one model per metric with the four predictors and a measurement-error term), posterior support, convergence, leave-one-site-out refits and variance inflation factors. | `site_replicates_n7.rds`, `site_level_n7.rds`, `brms_site_level_fits.rds`, `support_table_n7.rds`, `convergence_n7.rds`, `loo_site_results.rds` |
| 6 | `R/GenerateFigures.R` | Figures 2 and 3, partial-residual figures for the four predictors and posterior predictive checks. | `Figures/Bayesian_*.png` |
| 7 | `R/CompareModularityAlgorithms.R` | Modularity with Infomap vs spinglass (supplementary figure). | `Figures/modularity_infomap_vs_spinglass_comparison.png` |
| 8 | `R/RandomNetworksSVD.R` | SVD entropy of the observed food webs vs random networks (supplementary figure). | `Figures/random_vs_empirical_SVDE.png` |

Shared code:

- `R/network_fun.r`: functions for the metaweb assembly simulations, the stability simulations and the density plots of Figures S3 to S5.
- `R/site_replicates.R`: builds the replicate-level data of step 5 (200 metaweb-derived replicates per food web, seed 2025).

The manuscript is `Manuscript.Rmd` and the supplementary material `SuppMat.md`; both are rendered with `rmarkdown::render()`. Figure 1 (`Figures/Fig1_mapa.pdf`) and the schematic of SVD complexity and rank deficiency (`Figures/svd_complexity_figure_v4.svg`) are not generated by these scripts.

### Requirements

R packages: `tidyverse`, `igraph`, `readxl`, `brms` (with `cmdstanr`), `posterior`, `bayesplot`, `patchwork`, `ggtext`, `GGally`, `xtable`, `lhs`, `abc`, `future`, `future.apply`, `furrr`, `sf`, `terra`, `exactextractr`, and from GitHub `lsaravia/multiweb` and `lsaravia/meweasmo`.

## Data

- `Data/*_links_original_*.xlsx`: original link lists of the seven food webs. San Jorge Gulf, Burdwood Bank and Potter Cove are used as they are.
- `Data/*_links_standarized.csv`, `Data/*_links_standardized.csv`: standardized versions of Beagle Channel, Weddell Sea and the Scotia Sea food webs, used in the analysis.
- `Data/Links_data_2024_06_24_Time_04_31_40.xlsx`: name updates used by `R/compare_networks_refined.R`.
- `Data/ssp245_current.tif`: cumulative human impact, current layer of the SSP2-4.5 scenario (Halpern et al. 2025, <https://knb.ecoinformatics.org/view/doi%3A10.5063%2FF18K77KZ>).

In the result files the standardized food webs keep the site identifiers `Std BeagleChannel`, `Std Weddell Sea`, `Std NorthernScotia` and `Std SouthernScotia`.

## Main result files in `Data/`

- `network_info.rds`: the seven food webs (`igraph`), their empirical metrics and the predictors.
- `simulations_metaWeb_metrics.rds`: metrics of the 1,000 simulated food webs per site (modularity with Infomap).
- `stability_results_meig.rds`: maximum eigenvalues of the stable community matrices.
- `site_replicates_n7.rds`: the 200 replicates per food web used in the regressions, with the transformed responses.
- `site_level_n7.rds`: site-level means and standard deviations of each transformed metric, with the predictors.
- `brms_site_level_fits.rds`: the seven fitted `brms` models.
- `support_table_n7.rds`: posterior probability of a consistent sign for each metric and predictor.
- `loo_site_results.rds`: posterior support in each leave-one-site-out refit.
- `convergence_n7.rds`: R-hat, effective sample sizes and divergent transitions of the seven models.
