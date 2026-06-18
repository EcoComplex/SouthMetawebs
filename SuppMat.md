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

## Food Web Standardization

### Beagle Channel

The standardized version increases the resolution of the top trophic level by adding 21 species of seabirds and marine mammals absent from the original network, along with their diets. These include penguins (*Spheniscus magellanicus*, *Pygoscelis papua*), albatrosses (*Thalassarche melanophris*), cormorants (*Phalacrocorax atriceps*, *P. magellanicus*), gulls (*Larus dominicanus*, *L. scoresbii*), skuas (*Stercorarius chilensis*), terns (*Sterna hirundinacea*), petrels (*Macronectes giganteus*), cetaceans (*Cephalorhynchus commersonii*, *Grampus griseus*, *Lagenorhynchus australis*, *L. cruciger*, *L. obscurus*, *Lissodelphis peronii*, *Megaptera novaeangliae*, *Phocoena spinipinnis*, *Pseudorca crassidens*), otariids (*Otariidae*), and the marine otter *Lontra provocax*.

**Summary of structural changes:**

```{=latex}
\begin{center}
\begin{tabular}{lcccc}
\toprule
Version & S & L & Connectance & Link Density\\
\midrule
Original & 145 & 1,115 & 0.0530 & 7.69\\
Standardized & 166 & 1,304 & 0.0473 & 7.86\\
\bottomrule
\end{tabular}
\end{center}
```

The increase from 1,115 to 1,304 links (+185) reflects exclusively the addition of trophic interactions involving the 21 new top-level taxa.

### Weddell Sea

The standardized version reduces the resolution of basal producers by collapsing 62 individually resolved phytoplankton species — primarily diatoms (*Chaetoceros*, *Fragilariopsis*, *Thalassiosira*, *Nitzschia*, *Proboscia*, *Pseudo-Nitzschia*, *Porosira*, *Rhizosolenia*, *Trichotoxon*, and others) — into two functional groups: *Bacillariophyceae* and *Phytoplankton_other*. This reduction reflects the higher taxonomic resolution of the original Weddell Sea network relative to the other networks in the dataset, and brings basal resolution in line with the remaining study sites.

**Summary of structural changes:**

```{=latex}
\begin{center}
\begin{tabular}{lcccc}
\toprule
Version & S & L & Connectance & Link Density\\
\midrule
Original & 490 & 16,041 & 0.0668 & 32.74\\
Standardized & 430 & 11,284 & 0.0610 & 26.24\\
\bottomrule
\end{tabular}
\end{center}
```

The reduction from 490 to 430 species (−60) and from 16,041 to 11,284 links (−4,757 net) is entirely attributable to the collapse of 62 phytoplankton taxa into 2 functional groups. The 112 links added in the standardized version correspond to interactions reassigned to the two new functional groups that were not directly recoverable from the original link list.

### Northern Scotia Sea

The standardized version increases the resolution of the top trophic level by adding 35 species of seabirds and marine mammals absent from the original network, along with their diets.

**Summary of structural changes:**

```{=latex}
\begin{center}
\begin{tabular}{lcccc}
\toprule
Version & S & L & Connectance & Link Density\\
\midrule
Original & 218 & 10,008 & 0.2106 & 45.91\\
Standardized & 253 & 10,561 & 0.1650 & 41.74\\
\bottomrule
\end{tabular}
\end{center}
```

The increase from 218 to 253 species (+35) and 10,008 to 10,561 links (+553) reflects the addition of the 35 new top-level taxa and their interactions.

### Southern Scotia Sea

Identical to Northern Scotia Sea: the same 35 species of seabirds and marine mammals were added with their respective diets (see Northern Scotia Sea section above for the full species list).

**Summary of structural changes:**

```{=latex}
\begin{center}
\begin{tabular}{lcccc}
\toprule
Version & S & L & Connectance & Link Density\\
\midrule
Original & 192 & 7,241 & 0.1964 & 37.71\\
Standardized & 227 & 7,730 & 0.1500 & 34.05\\
\bottomrule
\end{tabular}
\end{center}
```

The increase from 192 to 227 species (+35) and 7,241 to 7,730 links (+489) is entirely attributable to the added top trophic level taxa.

\newpage


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

### Connectance (C)

Connectance quantifies the proportion of realized trophic interactions relative to all possible ones:

\begin{equation}
C = \frac{L}{S^2}, \quad L = \sum_{i,j} a_{ij}.
\end{equation}

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

### Link Density (LD)

Link density measures the average number of trophic links per species:

\begin{equation}
LD = \frac{L}{S}.
\end{equation}


###  Modularity (Q)

Modularity describes the extent to which the network is organized into modules with dense intra-module links and sparse inter-module links.
We detect modules using the **Infomap algorithm**, which minimizes the description length of a random walker’s trajectory to reveal community structure [@Rosvall2008].

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

### Estimating Network Complexity

We characterized structural complexity using two complementary metrics derived from the adjacency matrix. The first, SVD complexity, is based on the 
distribution of singular values obtained by Singular Value Decomposition (SVD) of the adjacency matrix. Ecological networks with strong trophic hierarchy, 
modularity, or body-size constraints tend to concentrate structural information in a few dominant dimensions, resulting in an uneven distribution of singular values. We captured this organized heterogeneity as $E = 1 - J$, where $J$ is the normalized Shannon entropy of the singular value spectrum-- Pielou's evenness [@Pielou1975]. Higher values of $E$ indicate that network structure is dominated by fewer independent dimensions, reflecting stronger ecological constraints, whereas random networks maximize $J$ and thus minimize $E$ (Figure S5).

The second metric, **rank deficiency** ($D$), quantifies the proportion of linearly dependent rows and columns in the adjacency matrix, relative to the maximum possible rank (number of trophic species). A fully ranked matrix ($D = 0$) implies that every species has a unique interaction profile, whereas high deficiency ($D \to 1$) indicates substantial redundancy in trophic strategies. Together, SVD complexity and rank deficiency provide complementary views of the external and internal dimensionality of food web structure [@Strydom2021]. For a graphical illustration of how these metrics capture different aspects of network organization, see Figure S1.

#### SVD Complexity

Singular Value Decomposition (SVD) factorizes the adjacency matrix $\mathbf{A}$ as

$$
\mathbf{A} = \mathbf{U} \, \boldsymbol{\Sigma} \, \mathbf{V}^T,
$$

where $\mathbf{U}$ and $\mathbf{V}$ are orthogonal matrices of left- and 
right-singular vectors, and $\boldsymbol{\Sigma}$ is a diagonal matrix 
containing the non-negative singular values $\sigma_i$ in descending order. 
Each singular value represents the weight of one independent structural 
dimension of the network.

Normalized singular values are defined as:

$$
s_i = \frac{\sigma_i}{\sum_{j=1}^{k} \sigma_j},
$$

where $k$ is the number of non-zero singular values. The SVD entropy $J$ is 
then computed following @Shannon1948:

$$
J = -\frac{1}{\ln(k)} \sum_{i=1}^{k} s_i \ln(s_i),
$$

where the normalization by $\ln(k)$ ensures a maximum of 1, equivalent to 
Pielou's evenness index [@Pielou1975]. We define SVD complexity as $E = 1 - J$.

A random network with fixed connectance distributes structural information 
evenly across many dimensions (high $J$, low $E$), whereas ecological networks 
with trophic hierarchy or modular organization concentrate most information in 
a few dimensions (low $J$, high $E$; see Figure S4).

#### Rank Deficiency

The rank $r$ of the adjacency matrix is the number of linearly independent 
rows (or columns). For a square matrix of dimension $M$ (the number of trophic 
species), rank deficiency is defined as:

$$
D = 1 - \frac{r}{M},
$$

where $D = 0$ indicates a full-rank matrix (all species have unique interaction 
profiles) and $D \to 1$ indicates high redundancy in trophic strategies. 
Dividing by $M$ controls for differences in species richness across networks, 
enabling cross-system comparisons [@Strydom2021].

#### Dynamic stability

We used the eigenvalue with the maximum real part of the community matrix Jacobian, $\lambda_{\max}$ [@Allesina2015a], for randomly parameterized systems, preserving the predator–prey (sign) structure and conditioning on stability [@Barabas2017]. This corresponds to the rightmost eigenvalue in the complex plane and determines local asymptotic stability. The system is stable when $\lambda_{\max} < 0$, and more negative values indicate faster return to equilibrium following perturbations.  

The entries of the Jacobian were sampled from uniform distributions: positive effects of prey on predators were bounded above by 1, and negative effects of predators on prey were bounded below by $-10$ (reflecting a 10% ecological efficiency). Diagonal entries representing self-regulation were sampled from a uniform distribution between $-\textit{SelfReg}$ and $0$. The value of $\textit{SelfReg}$ was calibrated via simulation so that 5% of randomly parameterized systems were stable (i.e., $\lambda_{\max} < 0$).  

To estimate $\textit{SelfReg}$, we performed a grid search over values from $-7$ to $-60$ (step size = 1). For each value, we generated 10,000 Jacobians and computed the proportion of stable systems. We then interpolated these results to identify the value of $\textit{SelfReg}$ yielding a stability probability of 0.05.  

Using this calibrated value, we generated 10,000 additional realizations and retained only stable systems to obtain the distribution of $\lambda_{\max}$ for each network. Because $\lambda_{\max} < 0$ for all retained systems, we multiplied values by $-1$ and modeled $\log(-\lambda_{\max})$, placing the stability metric on a positive, approximately log-normal scale.


\newpage

![Figure explaining SVD complexity and Rank Deficiency](Figures/svd_complexity_figure_v4.png)

\newpage


![Distribution of metaweb-simulated network metrics (SVD complexity and modularity) with the observed value indicated by the red dashed line.](Figures/empirical_vs_simulated_metrics_SVDC_Mod.png) 

\newpage

![Distribution of metaweb-simulated network metrics (Connectance and Rank Deficiency) with the observed value indicated by the red dashed line.](Figures/empirical_vs_simulated_metrics_C_RD.png) 

\newpage

![Distribution of metaweb-simulated network metrics (Link Density and Trophic Level) with the observed value indicated by the red dashed line.](Figures/empirical_vs_simulated_metrics_LD_TL.png)

\newpage

\begin{table}[ht]
\centering
\scriptsize
\begin{tabular}{lccccccc}
  \toprule
site & C & LD & TLmean & SVDComplexity & RankDeficiency & Modularity & Max. Eigenvalue\\ 
  \midrule
BurdwoodBank & \shortstack{0.014 \\ {\tiny(0.012--0.018)}} & \shortstack{4.463 \\ {\tiny(3.937--4.810)}} & \shortstack{2.470 \\ {\tiny(2.377--2.551)}} & \shortstack{0.079 \\ {\tiny(0.075--0.082)}} & \shortstack{0.635 \\ {\tiny(0.572--0.660)}} & \shortstack{0.334 \\ {\tiny(0.301--0.371)}} & \shortstack{-0.083 \\ {\tiny(-0.233---0.006)}} \\ 
  GulfSanJorge & \shortstack{0.038 \\ {\tiny(0.035--0.042)}} & \shortstack{5.590 \\ {\tiny(4.750--6.121)}} & \shortstack{2.906 \\ {\tiny(2.653--3.025)}} & \shortstack{0.080 \\ {\tiny(0.075--0.084)}} & \shortstack{0.447 \\ {\tiny(0.419--0.481)}} & \shortstack{0.137 \\ {\tiny(0.074--0.212)}} & \shortstack{-0.153 \\ {\tiny(-0.454---0.006)}} \\ 
  PotterCove & \shortstack{0.055 \\ {\tiny(0.050--0.060)}} & \shortstack{5.461 \\ {\tiny(4.639--5.935)}} & \shortstack{2.233 \\ {\tiny(2.094--2.530)}} & \shortstack{0.091 \\ {\tiny(0.085--0.096)}} & \shortstack{0.493 \\ {\tiny(0.451--0.533)}} & \shortstack{0.164 \\ {\tiny(0.054--0.243)}} & \shortstack{-0.172 \\ {\tiny(-0.543---0.006)}} \\ 
  Std BeagleChannel & \shortstack{0.048 \\ {\tiny(0.045--0.052)}} & \shortstack{7.256 \\ {\tiny(6.433--7.855)}} & \shortstack{2.607 \\ {\tiny(2.480--2.798)}} & \shortstack{0.092 \\ {\tiny(0.088--0.096)}} & \shortstack{0.411 \\ {\tiny(0.386--0.439)}} & \shortstack{0.119 \\ {\tiny(0.058--0.180)}} & \shortstack{-0.180 \\ {\tiny(-0.514---0.005)}} \\ 
  Std NorthernScotia & \shortstack{0.165 \\ {\tiny(0.157--0.174)}} & \shortstack{38.458 \\ {\tiny(34.339--41.574)}} & \shortstack{3.549 \\ {\tiny(3.431--3.683)}} & \shortstack{0.172 \\ {\tiny(0.167--0.175)}} & \shortstack{0.455 \\ {\tiny(0.434--0.476)}} & \shortstack{0.076 \\ {\tiny(0.001--0.113)}} & \shortstack{-0.375 \\ {\tiny(-1.220---0.011)}} \\ 
  Std SouthernScotia & \shortstack{0.151 \\ {\tiny(0.143--0.158)}} & \shortstack{31.339 \\ {\tiny(27.805--33.969)}} & \shortstack{3.513 \\ {\tiny(3.374--3.649)}} & \shortstack{0.167 \\ {\tiny(0.162--0.171)}} & \shortstack{0.424 \\ {\tiny(0.400--0.447)}} & \shortstack{0.071 \\ {\tiny(0.049--0.120)}} & \shortstack{-0.344 \\ {\tiny(-1.062---0.012)}} \\ 
  Std Weddell Sea & \shortstack{0.063 \\ {\tiny(0.059--0.075)}} & \shortstack{24.084 \\ {\tiny(21.360--26.175)}} & \shortstack{3.443 \\ {\tiny(3.170--3.757)}} & \shortstack{0.155 \\ {\tiny(0.148--0.159)}} & \shortstack{0.707 \\ {\tiny(0.669--0.724)}} & \shortstack{0.219 \\ {\tiny(0.171--0.245)}} & \shortstack{-0.255 \\ {\tiny(-0.765---0.010)}} \\ 
   \bottomrule
\end{tabular}
\caption{Summary of network metrics across 1000 simulations for each metaweb. For Maximal eigenvalues, the variability is generated by randomizations of interaction strength. Values are means with 95\% intervals in parentheses.} 
\label{tab:metaweb_metrics}
\end{table}


\begin{table}[ht]
\centering
\begin{tabular}{llrrrr}
  \toprule
  Response & Predictor & Prob pos & Prob neg & Prob max & Direction \\
  \midrule
  Trophic Level   & Area (log) & 0.976 & 0.024 & 0.976 & positive \\
  SVD Complexity  & Area (log) & 0.971 & 0.029 & 0.971 & positive \\
  Link Density    & Area (log) & 0.960 & 0.040 & 0.960 & positive \\
  Connectance     & Area (log) & 0.959 & 0.041 & 0.959 & positive \\
  Modularity      & Area (log) & 0.041 & 0.959 & 0.959 & negative \\
  Modularity      & Species    & 0.949 & 0.051 & 0.949 & positive \\
  Stability       & Area (log) & 0.944 & 0.056 & 0.944 & positive \\
  SVD Complexity  & Latitude   & 0.942 & 0.058 & 0.942 & positive \\
  Stability       & Latitude   & 0.922 & 0.078 & 0.922 & positive \\
  Rank Deficiency & Species    & 0.922 & 0.078 & 0.922 & positive \\
  Connectance     & Latitude   & 0.922 & 0.079 & 0.922 & positive \\
  Connectance     & Species    & 0.083 & 0.917 & 0.917 & negative \\
  Link Density    & Latitude   & 0.915 & 0.085 & 0.915 & positive \\
  Stability       & Species    & 0.086 & 0.914 & 0.914 & negative \\
  Trophic Level   & Latitude   & 0.909 & 0.091 & 0.909 & positive \\
  \bottomrule
\end{tabular}
\caption{Posterior probabilities of directional effects for predictor--response
pairs retained in the multivariate Bayesian model (posterior probability of
consistent sign $\geq 0.90$). Columns show the probability that each slope is
positive (Prob pos) or negative (Prob neg), the maximum of the two (Prob max),
and the inferred direction. Pairs are ordered by decreasing Prob max.}
\label{tab:posterior_slopes}
\end{table}

\newpage

![SVD entropy of the observed food webs analyzed in this study (filled triangles) compared with directed random networks (filled circles) constructed over the same range of trophic species numbers (S). Lower entropy values indicate greater structural organization relative to random expectations. ](Figures/random_vs_empirical_SVDE.png)

\newpage

![Pairwise relationships among environmental and network covariates used in the statistical analyses. Panels show bivariate scatterplots (lower triangle), smoothed density distributions (diagonal), and Pearson correlation coefficients (upper triangle) for the number of trophic species (S), latitude, log-transformed area, mean human impact index, and mean depth. This exploratory analysis was used to assess potential collinearity among predictors prior to model fitting](Figures/correlation_network_info.png)

\newpage


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

Continuous predictors were centered and scaled. These included the number of trophic species ($S$), the natural logarithm of area, latitude, and the human impact index. Depth and temperature were initially included as covariates; however, depth was strongly correlated with log-transformed area, and temperature was strongly correlated with latitude (Figure S6). To avoid multicollinearity and improve parameter interpretability, both variables were excluded from the final model. This collinearity was further supported by strong posterior correlations among regression coefficients when the full model was fitted.

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

Because a fully joint multivariate model is computationally demanding, we used **200 simulated replicates per site**, which balances precision and computational feasibility. To assess the robustness of this subsampling strategy, we refitted the model using progressively larger subsets of 50, 100, 150, and 200 replicates per food web. Parameter estimates stabilized after 100 replicates, indicating that the results were insensitive to further increases in the number of simulations.


The model was fitted in **brms** [@Burkner2017] using the **CmdStan** backend, with  4 chains and 4000 iterations (2000 warmup).

\newpage


![Posterior predictive checks for the seven response variables. Each panel shows density overlays of observed (black) and posterior predicted (blue) values for each metric.](Figures/Bayesian_pp_checks.png)

\newpage



\newpage

![Marginal effects of log-transformed area on food-web structure and dynamic stability from the multivariate Bayesian model. Solid black lines show posterior mean predictions, with shaded regions representing 95% credible intervals. Colored points indicate observed values for each site, and vertical bars denote the 95% quantiles summarizing variability generated by the metaweb assembly model. Panels show the effects of area on (top left) connectance (C), (top right) rank deficiency, (bottom left) modularity, and (bottom right) log dynamic stability (maximum eigenvalue). Overall, larger areas are associated with increased connectance and stability, alongside reduced modularity and rank deficiency, consistent with more interconnected and dynamically resilient food webs.](Figures/Bayesian_C_RD_Mod_ST_vs_Area_mv.png)


![Marginal effects of environmental and network size predictors on food-web structural metrics from the multivariate Bayesian model. Solid black lines show posterior mean predictions, with shaded regions indicating 95% credible intervals. Colored points represent observed values for each site, and vertical bars denote the 95% quantiles summarizing variability generated by the metaweb assembly model. Panels show the effects of latitude on (top left) connectance (C), (top right) rank deficiency, (bottom left) modularity, and (bottom right) log dynamic stability (maximum eigenvalue).](Figures/Bayesian_C_RD_Mod_ST_vs_Latitude_mv.png)


![Marginal effects of environmental and network size predictors on food-web structural metrics from the multivariate Bayesian model. Solid black lines show posterior mean predictions, with shaded regions indicating 95% credible intervals. Colored points represent observed values for each site, and vertical bars denote the 95% quantiles summarizing variability generated by the metaweb assembly model. Panels illustrate the effects of the number of trophic species on the full set of food-web metrics.](Figures/Bayesian_all_vs_S_mv.png)

![Marginal effects of environmental and network size predictors on food-web structural metrics from the multivariate Bayesian model. Solid black lines show posterior mean predictions, with shaded regions indicating 95% credible intervals. Colored points represent observed values for each site, and vertical bars denote the 95% quantiles summarizing variability generated by the metaweb assembly model. Panels illustrate the effects of the human impact index on the full set of food-web metrics.](Figures/Bayesian_all_vs_impact_mv.png)

![Residual correlations of network properties after accounting for fixed effects and site‐level shifts.](Figures/Bayesian_ResCor.png)

### References
