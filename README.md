# Ecosystem area and the structure of marine food webs along a latitudinal gradient in the Southwest Atlantic and Antarctica

[![DOI](https://zenodo.org/badge/1009246443.svg)](https://doi.org/10.5281/zenodo.21045337)

Code and data for the manuscript.

## Abstract

Understanding how food web structure varies across broad spatial gradients is essential for identifying general principles of ecosystem organization and resilience under global change. We analyzed seven highly resolved marine food webs spanning a latitudinal gradient from 48°S to 78°S in the Southwest Atlantic and Antarctic region to evaluate how ecosystem area, latitude, human impact, and number of trophic species are associated with network structure and dynamic stability. We quantified connectance, link density, mean trophic level, modularity, SVD complexity, rank deficiency, and dynamic stability, and related each metric to the four predictors using site-level Bayesian regressions that incorporate the uncertainty derived from a metaweb assembly model as a measurement-error term. Because the dataset comprises only seven independent food webs, we assessed the robustness of each relationship by refitting the models with each food web excluded in turn.

Ecosystem area showed the most consistent associations: larger ecosystems had higher mean trophic level, link density, connectance, and SVD complexity, and these associations held in all seven exclusions. The number of trophic species was positively associated with rank deficiency and modularity and negatively with connectance, but these associations lost support when a single food web, Burdwood Bank, was excluded. Latitude was positively associated with link density and SVD complexity, although these associations depended on three of the seven food webs, and its associations with connectance and trophic level did not survive exclusion. Human impact showed no clear association with any metric, possibly because the cumulative impact index integrates stressors whose effects on food web structure differ in magnitude and direction.

Among the predictors considered, ecosystem area showed the most consistent associations with marine food web structure in the Southwest Atlantic and Antarctic. Patterns associated with species richness and latitude should be treated as exploratory and require evaluation with additional independent food webs.

## Workflow

The analyses run in this order; each step reads the `.rds` files written by the previous ones in `Data/`.

1. `R/compare_networks_refined.R`: compares the original and standardized food webs (added top predators in Beagle Channel and the Scotia Sea; phytoplankton collapsed into functional groups in the Weddell Sea).
2. `LoadNetworksDescription.Rmd`: loads the seven food webs and describes them (trophic levels and modules).
3. `HalpernHumanImpact.Rmd`: computes the mean human impact index for each area.
4. `metawebassemblysimulations.Rmd`: fits the metaweb assembly model by approximate Bayesian computation (package `meweasmo`), simulates local networks and computes their metrics (Figures S3 to S5 and S7).
5. `Stabilitysimulations.Rmd`: dynamic stability (maximum eigenvalue of randomly parameterized community matrices).
6. `BayesianSiteLevel.Rmd`: site-level Bayesian regressions (`brms`, one model per metric with the four predictors and a measurement-error term), posterior support, convergence diagnostics, leave-one-site-out refits and variance inflation factors.
7. `GenerateFigures.R`: Figures 2 and 3 and the supplementary figures of partial effects and posterior predictive checks.
8. `CompareModularityAlgorithms.R`: Infomap versus spinglass modularity (supplementary figure).

`R/network_fun.r` holds shared functions and `R/generate_figure_SVD.R` the schematic of SVD complexity and rank deficiency. The manuscript is `Manuscript.Rmd` and the supplementary material `SuppMat.md`; both are rendered with `rmarkdown::render()`.

## Main result files in `Data/`

- `site_level_n7.rds`: site-level means and standard deviations of each transformed metric, with the predictors.
- `brms_site_level_fits.rds`: the seven fitted `brms` models.
- `support_table_n7.rds`: posterior probability of a consistent sign for each metric and predictor.
- `loo_site_results.rds`: posterior support in each leave-one-site-out refit.
- `convergence_n7.rds`: R-hat, effective sample sizes and divergent transitions of the seven models.
