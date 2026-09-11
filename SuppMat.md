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

The standardized version reduces the resolution of basal producers by collapsing 62 individually resolved phytoplankton species, primarily diatoms (*Chaetoceros*, *Fragilariopsis*, *Thalassiosira*, *Nitzschia*, *Proboscia*, *Pseudo-Nitzschia*, *Porosira*, *Rhizosolenia*, *Trichotoxon*, and others), into two functional groups: *Bacillariophyceae* and *Phytoplankton_other*. This reduction reflects the higher taxonomic resolution of the original Weddell Sea network relative to the other networks in the dataset, and brings basal resolution in line with the remaining study sites.

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

<!-- FIGURA: ya generada, correr CompareModularityAlgorithms.R si hace falta reproducirla. Requiere
Data/simulations_metaWeb_metrics.rds (Infomap) y Data/simulations_metaWeb_metrics_spinglass.rds
(ya calculado). El resumen numerico (Data/modularity_infomap_vs_spinglass_summary.csv) no esta en
la carpeta Data -- si lo segui necesitando, correr el script de nuevo lo regenera junto con la
figura. -->

![Comparison of modularity (Newman-Girvan Q) distributions obtained with the Infomap algorithm (used throughout this study) and with the spinglass algorithm, computed on the same metaweb-derived replicates for each of the seven food webs. For most sites the two algorithms produce closely overlapping distributions. For Northern Scotia, spinglass returns a disproportionate fraction of replicates with modularity at or near zero that Infomap does not, indicating that community detection is comparatively less stable for this network regardless of the algorithm used.](Figures/modularity_infomap_vs_spinglass_comparison.png)

### Estimating Network Complexity

We characterized structural complexity using two complementary metrics derived from the adjacency matrix. The first, SVD complexity, is based on the 
distribution of singular values obtained by Singular Value Decomposition (SVD) of the adjacency matrix. Ecological networks with strong trophic hierarchy, 
modularity, or body-size constraints tend to concentrate structural information in a few dominant dimensions, resulting in an uneven distribution of singular values. We captured this organized heterogeneity as $E = 1 - J$, where $J$ is the normalized Shannon entropy of the singular value spectrum-- Pielou's evenness [@Pielou1975]. Higher values of $E$ indicate that network structure is dominated by fewer independent dimensions, reflecting stronger ecological constraints, whereas random networks maximize $J$ and thus minimize $E$ (Figure S6).

The second metric, **rank deficiency** ($D$), quantifies the proportion of linearly dependent rows and columns in the adjacency matrix, relative to the maximum possible rank (number of trophic species). A fully ranked matrix ($D = 0$) implies that every species has a unique interaction profile, whereas high deficiency ($D \to 1$) indicates substantial redundancy in trophic strategies. Together, SVD complexity and rank deficiency provide complementary views of the external and internal dimensionality of food web structure [@Strydom2021]. For a graphical illustration of how these metrics capture different aspects of network organization, see Figure S2.

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
a few dimensions (low $J$, high $E$; see Figure S2).

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

<!-- TABLA S2 ACTUALIZADA: reemplaza la vieja tabla del modelo multivariado (tab:posterior_slopes).
Incluye ahora la columna de leave-one-site-out (fallas de 7 exclusiones) y el criterio de robustez
"falla en <=3 de 7 exclusiones". Para las 5 relaciones de Área los conteos exactos de fallas (2 o 3
de 7) surgen directamente del texto del Manuscript.Rmd. Para las 7 filas de Latitud/Especies, el
texto original solo reporta la cota "fallan en al menos 5 de 7"; no tengo el conteo exacto por fila
(6 vs 7, por ejemplo) porque no quedó registrada la tabla de fragilidad completa que pegaste en el
chat. Dejé ">=5" en esas filas, consistente con lo ya escrito en el manuscrito -- si tenés a mano el
conteo exacto por efecto, decímelo y lo actualizo. -->

\begin{table}[ht]
\centering
\scriptsize
\begin{tabular}{llrrrrcc}
  \toprule
  Response & Predictor & Prob pos & Prob neg & Prob max & Direction & LOSO failures (of 7) & Robust ($\leq$3 failures) \\
  \midrule
  Trophic Level    & Area (log) & 0.974 & 0.026 & 0.974 & positive & 2        & Yes \\
  SVD Complexity   & Area (log) & 0.965 & 0.035 & 0.965 & positive & 2        & Yes \\
  Connectance      & Area (log) & 0.956 & 0.044 & 0.956 & positive & 3        & Yes \\
  Link Density     & Area (log) & 0.949 & 0.051 & 0.949 & positive & 2        & Yes \\
  Modularity       & Area (log) & 0.069 & 0.931 & 0.931 & negative & 7        & No \\
  SVD Complexity   & Latitude   & 0.939 & 0.061 & 0.939 & positive & $\geq$5  & No \\
  Modularity       & Species    & 0.938 & 0.062 & 0.938 & positive & $\geq$5  & No \\
  Connectance      & Species    & 0.075 & 0.925 & 0.925 & negative & $\geq$5  & No \\
  Connectance      & Latitude   & 0.923 & 0.077 & 0.923 & positive & $\geq$5  & No \\
  Rank Deficiency  & Species    & 0.914 & 0.086 & 0.914 & positive & $\geq$5  & No \\
  Link Density     & Latitude   & 0.908 & 0.092 & 0.908 & positive & $\geq$5  & No \\
  Trophic Level    & Latitude   & 0.900 & 0.100 & 0.900 & positive & $\geq$5  & No \\
  \bottomrule
\end{tabular}
\caption{Posterior probabilities of directional effects for predictor--response
pairs reaching the nominal 0.90 support threshold (posterior probability of
consistent sign) in the site-level Bayesian regressions ($n=7$ food webs per
model), together with the number of leave-one-site-out refits (of seven) in
which the effect failed to reach this threshold with a consistent sign. An
effect is classified as robust when it failed to hold in no more than three of
the seven exclusions. Pairs are ordered by decreasing robustness and then by
decreasing Prob max. Effects of dynamic stability on area, latitude and
species richness, and all effects of the human impact index, did not reach the
0.90 support threshold under the full dataset and are not shown here (see
Results).}
\label{tab:posterior_slopes}
\end{table}

\begin{table}[ht]
\centering
\begin{tabular}{lr}
  \toprule
  Predictor & VIF \\
  \midrule
  Number of trophic species (S) & 3.10 \\
  Area (log)                    & 1.91 \\
  Latitude                      & 1.81 \\
  Human impact (mean)           & 1.68 \\
  \bottomrule
\end{tabular}
\caption{Variance inflation factors (VIF) for the four predictors used in the
site-level Bayesian regressions ($n=7$ food webs), computed as
$1/(1-R^2_l)$ from an ordinary least-squares regression of each standardized
predictor $l$ on the remaining three. All values are well below conventional
thresholds for concerning collinearity (commonly VIF $>5$--10).}
\label{tab:vif}
\end{table}

\newpage

![SVD entropy of the observed food webs analyzed in this study (filled triangles) compared with directed random networks (filled circles) constructed over the same range of trophic species numbers (S). Lower entropy values indicate greater structural organization relative to random expectations. ](Figures/random_vs_empirical_SVDE.png)

\newpage

<!-- FIGURA A REGENERAR: correr a mano en RStudio el chunk "CheckCorrelation_createTables" de
metawebassemblysimulations.Rmd (tiene eval=FALSE, no corre solo al knitear). Ese chunk agrega
temp_C via join por `name`, guarda network_info.rds, y escribe la figura con
select(temp_C, S, latitude, log_area, impact_mean, depth_m). La version actual del PNG todavia
no incluye S. El caption ya esta actualizado para las 6 variables. -->

![Pairwise relationships among environmental and network covariates considered prior to model fitting. Panels show bivariate scatterplots (lower triangle), smoothed density distributions (diagonal), and Pearson correlation coefficients (upper triangle) for mean sea surface temperature, the number of trophic species (S), latitude, log-transformed area, mean human impact index, and mean depth. Temperature and depth were strongly correlated with latitude and log-transformed area, respectively, and were excluded from the statistical model on that basis (see Supplementary Methods S1).](Figures/correlation_network_info.png)

\newpage


# **Supplementary Methods S1: Site-Level Bayesian Measurement-Error Models**

Let $j = 1,\dots,7$ index the seven food webs (sites) analyzed in this study. For each site $j$ and each of $K = 7$ network metrics, we summarized the distribution of metaweb-derived replicate values (see Methods) into a site-level point estimate and its associated uncertainty, after applying the transformations described below.

$$
\mathbf{y}_{j} =
\big(
C_{j}, \text{SVD}_{j}, M_{j}, TL_{j}, LD_{j}, RD_{j}, ME_{j}
\big),
$$

representing respectively: connectance, SVD complexity, modularity, mean trophic level, link density, rank deficiency, and the dynamic stability metric (maximum eigenvalue).

Following the principle of using monotone transformations that preserve order and scale structure, we applied the following transformations to each metric, computed across the retained metaweb-derived replicates of each site:

1. Proportions → logit-transform

$$
C^{\ast} = \operatorname{logit}(C), \qquad
M^{\ast} = \operatorname{logit}(M), \qquad
RD^{\ast} = \operatorname{logit}(RD).
$$

2. Positive-valued variables → log-transform

$$
LD^{\ast} = \log(LD), \qquad
ME^{\ast} = \log(-ME),
$$

where the stability metric is the *negative* maximum eigenvalue, so $-ME > 0$.

3. Approximately symmetric metrics → standardization

$$
\text{SVD}^{\ast} = \frac{\text{SVD} - \mu_{\text{SVD}}}{\sigma_{\text{SVD}}},
\qquad
TL^{\ast} = \frac{TL - \mu_{TL}}{\sigma_{TL}} .
$$

For each site $j$ and each transformed metric $k$, we then computed the replicate mean $\bar y_{jk}$ and replicate standard deviation $s_{jk}$ across the retained metaweb-derived replicates, floored at a small positive value ($10^{-6}$) to avoid degenerate zero-uncertainty sites. These site-level pairs $(\bar y_{jk}, s_{jk})$, with $n=7$ per metric, are the unit of analysis for all regression models described below.

Continuous predictors were centered and scaled. These included the number of trophic species ($S$), the natural logarithm of area, latitude, and the human impact index. Depth and temperature were initially considered as covariates; however, depth was strongly correlated with log-transformed area, and temperature was strongly correlated with latitude (Figure S7). To avoid multicollinearity and improve parameter interpretability, both variables were excluded from the final model.

For each site $j$:

$$
\mathbf{x}_{j} =
\big(
S_{j},
\log A_{j},
\text{Lat}_{j},
\text{HumImp}_{j}
\big),
$$

all standardised to mean 0 and unit variance.

### Site-level measurement-error model

Each of the $K=7$ response metrics was modeled **independently**, with one row of data per site ($n=7$). For metric $k$, the site-level mean $\bar y_{jk}$ is modeled as:

$$
\bar y_{jk} \sim \mathcal{N}\!\left(\mu_{jk},\; \sqrt{s_{jk}^{2} + \sigma_k^{2}}\right),
$$

$$
\mu_{jk} = \beta_{0k} + \mathbf{x}_{j}\boldsymbol{\beta}_k ,
$$

where $s_{jk}$ is the (known, data-derived) replicate standard deviation for site $j$ and metric $k$, incorporated directly as a measurement-error term via the `se(sd_y, sigma = TRUE)` specification in `brms` [@Burkner2017], and $\sigma_k$ is an additional residual standard deviation estimated from the data, capturing site-level variability not explained by the predictors or by replicate-derived measurement error. Because there is exactly one observation per site once replicates are summarized, no site-level grouping term (random intercept) is included.

Each of the seven metrics ($k = 1,\dots,7$) was fit as a **separate univariate model**. We modeled the seven responses independently rather than jointly with a shared residual correlation matrix, since a joint multivariate specification would require estimating a full $7\times7$ residual correlation matrix (21 pairwise correlations) from only seven independent food webs, which the data cannot identify. Correlations among metrics are therefore not estimated or reported.

We used weakly informative priors for regularized estimation:

$$
\beta_{kl} \sim \mathcal{N}(0, 1), \qquad
\beta_{0k} \sim \mathcal{N}(0, 1), \qquad
\sigma_k \sim \text{Exponential}(1).
$$

Standardization ensures regression coefficients lie on comparable scales and priors correspond to plausible effect sizes.

Each of the seven models was fitted in **brms** [@Burkner2017] using the **CmdStan** backend, with 4 chains and 4,000 iterations (2,000 warmup), `adapt_delta = 0.95` and `max_treedepth = 15`. Convergence was assessed with R-hat and effective sample size diagnostics for all parameters.

### Posterior support

For each predictor–response pair, we summarized the posterior distribution of the corresponding slope as the posterior probability of a consistent sign:

$$
p_{kl} = \max\!\big(\Pr(\beta_{kl} > 0 \mid \mathbf{y}),\; \Pr(\beta_{kl} < 0 \mid \mathbf{y})\big),
$$

and treated $p_{kl} \geq 0.90$ as the nominal threshold for discussing an effect (Table S2).

### Leave-one-site-out robustness check

Because $n = 7$, any relationship reaching the nominal support threshold under the full dataset could in principle be driven by a small number of influential food webs. For every predictor–response pair with $p_{kl} \geq 0.90$ under the full dataset, we refit the corresponding univariate model seven times, each time excluding one site ($n = 6$ per refit), and recorded whether the sign of the effect and $p_{kl} \geq 0.90$ were preserved. An effect was classified as **robust** if this held in at least four of the seven leave-one-site-out refits (i.e., it failed to hold in at most three exclusions); effects failing to hold in four or more refits are reported in Table S2 but are not interpreted as general patterns in the main text.

### Variance inflation factors

Collinearity among the four standardized predictors ($S$, log-area, latitude, human impact) was assessed at the site level ($n=7$) using variance inflation factors,

$$
VIF_l = \frac{1}{1 - R^2_l},
$$

where $R^2_l$ is the coefficient of determination from an ordinary least-squares regression of predictor $l$ on the remaining three predictors. Resulting values (Table S3) ranged from 1.68 to 3.10, well below conventional thresholds for concerning collinearity.

\newpage

<!-- FIGURA A REGENERAR: correr GenerateFigures.R (pp_check() sobre cada uno de los 7 modelos
univariados de fits_site). El caption ya esta actualizado para el modelo site-level. -->

![Posterior predictive checks for the seven site-level Bayesian regression models. Each panel shows density overlays of observed (black) and posterior predicted (blue) values for one network metric.](Figures/Bayesian_pp_checks.png)

\newpage

<!-- FIGURA A REGENERAR: correr GenerateFigures.R (bloque 12, marginal_plot_orig() para C_logit,
Rank_logit, Mod_logit y ME_log vs log_area_s). -->

![Marginal effects of log-transformed area on connectance, rank deficiency, modularity, and dynamic stability, in each metric's original units. Solid black lines show posterior mean predictions from the site-level Bayesian regressions, back-transformed from the model scale; shaded regions represent 95% credible intervals. Colored points show the observed mean for each food web, computed directly from the metaweb-derived replicates, with vertical bars denoting their 95% quantiles. Connectance is the only one of these four relationships that reached the 0.90 posterior support threshold and remained robust to leave-one-site-out exclusion (Table S2); it is also shown, together with the other area-robust relationships, in the main text (Figure 3). The negative association between area and modularity reached the nominal 0.90 threshold under the full dataset but did not survive exclusion of any single food web and is not interpreted as a genuine effect (see Discussion). Rank deficiency and dynamic stability did not reach the 0.90 support threshold for area in the site-level model.](Figures/Bayesian_C_RD_Mod_ST_vs_Area_mv_origscale.png)

<!-- FIGURA A REGENERAR: correr GenerateFigures.R (bloque 12, marginal_plot_orig() para C_logit,
Rank_logit, Mod_logit y ME_log vs latitude_s). -->

![Marginal effects of latitude on connectance, rank deficiency, modularity, and dynamic stability, in each metric's original units. Solid black lines show posterior mean predictions from the site-level Bayesian regressions, back-transformed from the model scale; shaded regions indicate 95% credible intervals. Colored points show the observed mean for each food web, computed directly from the metaweb-derived replicates, with vertical bars denoting their 95% quantiles. Connectance reached the 0.90 posterior support threshold for latitude under the full dataset but did not survive leave-one-site-out exclusion (Table S2); rank deficiency, modularity, and dynamic stability did not reach this threshold.](Figures/Bayesian_C_RD_Mod_ST_vs_Latitude_mv_origscale.png)

<!-- FIGURA A REGENERAR: correr GenerateFigures.R (bloque 12, marginal_plot_orig() para las 7
metricas vs S_s). -->

![Marginal effects of the number of trophic species (S) on the full set of food-web metrics, in each metric's original units. Solid black lines show posterior mean predictions from the site-level Bayesian regressions, back-transformed from the model scale; shaded regions indicate 95% credible intervals. Colored points show the observed mean for each food web, computed directly from the metaweb-derived replicates, with vertical bars denoting their 95% quantiles. Modularity, rank deficiency, and connectance reached the 0.90 posterior support threshold for the number of species under the full dataset, but none of these effects survived leave-one-site-out exclusion (Table S2) and none is interpreted as a general pattern.](Figures/Bayesian_all_vs_S_mv_origscale.png)

<!-- FIGURA A REGENERAR: correr GenerateFigures.R (bloque 12, marginal_plot_orig() para las 7
metricas vs impact_mean_s). -->

![Marginal effects of the human impact index on the full set of food-web metrics, in each metric's original units. Solid black lines show posterior mean predictions from the site-level Bayesian regressions, back-transformed from the model scale; shaded regions indicate 95% credible intervals. Colored points show the observed mean for each food web, computed directly from the metaweb-derived replicates, with vertical bars denoting their 95% quantiles. No predictor-response combination reached the 0.90 posterior support threshold (maximum posterior probability = 0.819, for SVD complexity).](Figures/Bayesian_all_vs_impact_mv_origscale.png)

<!-- FIGURA ELIMINADA: Bayesian_ResCor.png (matriz de correlacion residual) ya no corresponde: el
modelo site-level ajusta cada metrica de forma independiente y no estima una matriz de correlacion
residual conjunta (ver Supplementary Methods S1). -->

### References
