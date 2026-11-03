# Ecosystem area predicts marine food web complexity and trophic structure across the Southwest Atlantic and Antarctica

[![DOI](https://zenodo.org/badge/1009246443.svg)](https://doi.org/10.5281/zenodo.21045337)

Code and data for the manuscript.

## Abstract

**Aim:** To assess how ecosystem area, latitude, human impact and species richness are associated with the structure and dynamic stability of marine food webs.

**Location:** Southwest Atlantic Ocean and Antarctica, from San Jorge Gulf to the Weddell Sea (46°S to 76°S).

**Taxon:** Marine organisms across all trophic levels, from primary producers to top predators (fishes, seabirds and marine mammals).

**Methods:** We quantified connectance, link density, mean trophic level, SVD complexity, rank deficiency, modularity and dynamic stability in seven highly resolved food webs. We estimated their variability by simulating local food webs using a metaweb assembly model fitted via approximate Bayesian computation. Each metric was related to ecosystem area, latitude, human impact and species richness using a Bayesian regression, including the simulated variability as a measurement-error term. The robustness of each association was assessed by refitting the models with each food web excluded in turn.

**Results:** Ecosystem area was positively associated with mean trophic level, link density, connectance and SVD complexity, and these associations persisted across all seven food web exclusions. Species richness was positively associated with rank deficiency and modularity and negatively associated with connectance, while latitude was positively associated with link density and SVD complexity. All these associations lost support when Burdwood Bank was excluded, and some also when one or two other food webs were excluded. Human impact showed no clear association with any metric, and no predictor was associated with dynamic stability.

**Main conclusions:** Ecosystem area was the most consistent correlate of marine food web structure along the Southwest Atlantic-Antarctic gradient, whereas the patterns associated with species richness and latitude remain exploratory. In macroecological syntheses based on few, non-replicated food webs, sensitivity analyses that exclude individual systems in turn are needed to distinguish consistent relationships from those that depend on particular systems and their ecological context.

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
| 9 | `R/FoodChainLength.R` | Compares mean trophic level with food chain length (maximum and short-weighted trophic level) in the seven observed food webs (diagnostic). | `food_chain_length.csv` |

Shared code:

- `R/network_fun.r`: functions for the metaweb assembly simulations, the stability simulations and the density plots of Figures S3 to S5.
- `R/site_replicates.R`: builds the replicate-level data of step 5 (200 metaweb-derived replicates per food web, seed 2025).

The manuscript is `Manuscript.Rmd` and the supplementary material `SuppMat.md`; both are rendered with `rmarkdown::render()`. `R/make_submission.R` builds the anonymized `Main.Rmd` from `Manuscript.Rmd`, copies the title and author block into `Title.Rmd`, and renders both for submission. Figure 1 (`Figures/Fig1_mapa.pdf`) and the schematic of SVD complexity and rank deficiency (`Figures/svd_complexity_figure_v4.svg`) are not generated by these scripts.

### Requirements

R packages: `tidyverse`, `igraph`, `readxl`, `brms` (with `cmdstanr`), `posterior`, `bayesplot`, `patchwork`, `ggtext`, `GGally`, `NetIndices`, `xtable`, `lhs`, `abc`, `future`, `future.apply`, `furrr`, `sf`, `terra`, `exactextractr`, and from GitHub `lsaravia/multiweb` and `lsaravia/meweasmo`.

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
