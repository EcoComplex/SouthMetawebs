---
title: Supplementary Materials Latitudinal gradients in the structure and stability of marine food webs in the Southwest Atlantic
bibliography: SouthMarineFoodWebs.bib
csl: "functional-ecology.csl"
output:
  pdf_document:
    latex_engine: lualatex
    keep_tex: true
    includes:
      in_header: Appendices.sty
---
## Metric Formulas

The food web is represented by an **adjacency matrix** $\mathbf{A} = [a_{ij}]$, where each element is defined as:

$$
a_{ij} =
\begin{cases}
1, & \text{if species } i \text{ consumes species } j \\
0, & \text{otherwise}.
\end{cases}
$$

Here, rows represent **predators**, and columns represent **prey**. This directed binary matrix forms the basis for computing all structural metrics.


### Mean Trophic Level (MTL)

The trophic level of species *i* is defined recursively as:

$$
TL_i = 1 + \frac{1}{k_i} \sum_{j} a_{ij} \, TL_j,
\quad \text{where} \quad k_i = \sum_{j} a_{ij}.
$$

Then, the mean trophic level for the food web is:

\begin{equation}
MTL = \frac{1}{S} \sum_{i=1}^{S} TL_i.
\end{equation}

being $S$ the total number of trophic species in the food web.

### Connectance (C)

Connectance quantifies the proportion of realized trophic interactions relative to all possible ones:

\begin{equation}
C = \frac{L}{S^2}, \quad L = \sum_{i,j} a_{ij}.
\end{equation}


### Link Density (LD)

Link density measures the average number of trophic links per species:

\begin{equation}
LD = \frac{L}{S}.
\end{equation}


###  Modularity (Q)

Modularity describes the extent to which the network is organized into modules with dense intra-module links and sparse inter-module links.
We detect modules using the **Infomap algorithm**, which minimizes the description length of a random walker’s trajectory to reveal community structure (Rosvall and Bergstrom, 2008).

Formally, modularity is reported as the fraction of links within modules relative to the whole network, given by the Infomap solution.

A widely used definition for modularity in ecological networks is:

\begin{equation}
Q = \sum_{s=1}^{M} \left[ \frac{l_s}{L} - \left( \frac{d_s}{2L} \right)^2 \right]
\end{equation}


where:

* $M$ = total number of modules.
* $l_s$ = number of links within module *s*.
* $d_s$ = sum of degrees of nodes in module *s*.
* $L$ = total number of links in the network.

(This is the standard Newman-Girvan modularity.)

### References

* Infomap: Rosvall M., Bergstrom C.T. (2008). Maps of random walks on complex networks reveal community structure. *PNAS* 105(4): 1118–1123.

\newpage


![Distribution of metaweb-simulated network metrics (SVD complexity and modularity) with the empirical value indicated by the red dashed line.](Figures/empirical_vs_simulated_metrics_SVDC_Mod.png) 

![Distribution of metaweb-simulated network metrics (Connectance and Rank Deficiency) with the empirical value indicated by the red dashed line.](Figures/empirical_vs_simulated_metrics_C_RD.png) 

![Distribution of metaweb-simulated network metrics (Link Density and Trophic Level) with the empirical value indicated by the red dashed line.](Figures/empirical_vs_simulated_metrics_LD_TL.png)

![Pairwise relationships among environmental and network covariates used in the statistical analyses. Panels show bivariate scatterplots (lower triangle), smoothed density distributions (diagonal), and Pearson correlation coefficients (upper triangle) for the number of trophic species (S), latitude, log-transformed area, mean human impact index, and mean depth. This exploratory analysis was used to assess potential collinearity among predictors prior to model fitting](Figures/correlation_network_info.png)

# **Supplementary Methods S1 — Multivariate Bayesian Multilevel Model**

Let $i = 1,\dots,N$ index observations (metaweb-derived replicates) and $j = 1,\dots,J$ index sampling sites.
For each observation we measured a vector of $K = 7$ network metrics:

$$
\mathbf{y}_{ij} =
\big(
C_{ij}, \text{SVD}_{ij}, M_{ij}, TL_{ij}, LD_{ij}, RD_{ij}, ME_{ij}
\big),
$$

representing respectively: connectance, SVD complexity, modularity, mean trophic level, link density, rank deficiency, and the dynamic stability metric (Maximum eigenvalue).

Following the principle of using monotone transformations that preserve order and scale structure, we applied the following transformations to each metric

1. Proportions → logit-transform

$$
C^{\ast}_{ij} = \operatorname{logit}(C_{ij}), \qquad
M^{\ast}_{ij} = \operatorname{logit}(M_{ij}), \qquad
RD^{\ast}_{ij} = \operatorname{logit}(RD_{ij}).
$$

2. Positive-valued variables → log-transform

$$
LD^{\ast}_{ij} = \log(LD_{ij}), \qquad
ME^{\ast}_{ij} = \log(-ME_{ij}),
$$

where the stability metric is the *negative* maximum eigenvalue, so $-ME_{ij} > 0$.

3. Approximately symmetric metrics → standardization

$$
\text{SVD}^{\ast}_{ij} = \frac{\text{SVD}_{ij} - \mu_{\text{SVD}}}{\sigma_{\text{SVD}}},
\qquad
TL^{\ast}_{ij} = \frac{TL_{ij} - \mu_{TL}}{\sigma_{TL}} .
$$

Let $\mathbf{y}^{\ast}_{ij}$ denote the transformed response vector.

Continuous predictors were centered and scaled. Representing number of trophic species = $S$ , natural logarithm of the Area, Latitude, Human Impact Index. Depth was initially included as a covariate; however, it showed a strong correlation with log-transformed area (Figure S4). This collinearity was further confirmed after fitting the model with both predictors, by examining the posterior correlation between their regression coefficients. To avoid multicollinearity and improve interpretability of parameter estimates, depth was therefore excluded from the final model and log-transformed area was retained as the sole spatial size–related predictor.

For each observation:

$$
\mathbf{x}_{ij} =
\big(
S_{ij},
\log A_{ij},
\text{Lat}_{ij},
\text{HumImp}_{ij},
\big),
$$

all standardised to mean 0 and unit variance.

For each response $k \in {1,\dots,K}$ the likelihood is:

$$
y^{\ast}_{ijk} \sim \mathcal{N}(\mu_{ijk}, \sigma_k),
$$

$$
\mu_{ijk} = \alpha_{jk} + \mathbf{x}_{ij}\boldsymbol{\beta}_k ,
$$

where $\boldsymbol{\beta}_k$ = vector of regression coefficients for response $k$ and  $\alpha_{jk}$ = varying intercept for site $j$ and metric $k$

Each site has a $K$-dimensional vector of intercepts (varying effects):

$$ 
\boldsymbol{\alpha}_{j} =
(\alpha_{j1}, \dots, \alpha_{jK})
\sim
\mathcal{MVN}
\big(
\boldsymbol{\alpha}_{0},
\Sigma_{\alpha}
\big),
$$

allowing correlated among-site shifts across metrics.

Residuals across the seven responses share a multivariate normal structure:

$$
\boldsymbol{\varepsilon}_{ij}
= \mathbf{y}^{\ast}_{ij} - \boldsymbol{\mu}_{ij}
\sim
\mathcal{MVN}\left(\mathbf{0}, \Sigma_{\varepsilon}\right).
$$

Thus $\Sigma_{\varepsilon}$ is a full $7\times 7$ covariance matrix capturing residual co-variation among metrics after accounting for covariates and site effects.


We use weakly informative priors for regularised estimation:

$$
\beta_{kl} \sim \mathcal{N}(0, 1), \qquad
\alpha_{0k} \sim \mathcal{N}(0, 1).
$$

$$
\sigma_k \sim \text{Exponential}(1).
$$

$$
\Sigma_{\alpha} \sim \text{LKJcorr}(2), \qquad
\Sigma_{\varepsilon} \sim \text{LKJcorr}(2).
$$

Standardization ensures regression coefficients lie on comparable scales and priors correspond to plausible effect sizes.

Because a fully joint multivariate model is computationally demanding, we used **200 simulated replicates per site**, which balances precision and computational feasibility. The hierarchical structure still correctly models within-site dependence and avoids pseudo-replication.

The model was fitted in **brms** [@Burkner2017] using the **CmdStan** backend, with  4 chains and 4000 iterations (2000 warmup).

\newpage


![Posterior predictive checks for the seven response variables. Each panel shows density overlays of observed (black) and posterior predicted (blue) values for each metric.](Figures/Bayesian_pp_checks.png)


\newpage

![](Figures/Bayesian_C_RD_Mod_ST_vs_Area_mv.png)

