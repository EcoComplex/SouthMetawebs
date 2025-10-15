---
title: Supplementary Materials Latitudinal gradients in the structure and stability of marine food webs in the Southwest Atlantic
bibliography: SouthMarineFoodWebs.bib
csl: "functional-ecology.csl"
output:
  pdf_document:
    latex_engine: lualatex
    keep_tex: true
    includes:
---


|Metric         |Term           |    edf| p-value| PartialDeviance| TotalDevianceExplained|
|:--------------|:--------------|------:|-------:|---------------:|----------------------:|
|C              |s(S)           | 1.0003|  0.4121|          2.2282|                   89.4|
|C              |s(log_area)    | 1.0003|  0.7780|          0.2638|                   89.4|
|C              |s(latitude)    | 1.0496|  0.3471|          3.4241|                   89.4|
|C              |s(impact_mean) | 1.0006|  0.7286|          0.3990|                   89.4|
|C              |s(depth_m)     | 1.0001|  0.0019|         31.8231|                   89.4|
|C              |s(site)        | 0.9492|  0.0000|         61.8618|                   89.4|
|SVDComplexity  |s(S)           | 1.0002|  0.6288|          0.6628|                   94.4|
|SVDComplexity  |s(log_area)    | 1.0002|  0.3721|          2.2609|                   94.4|
|SVDComplexity  |s(latitude)    | 1.0363|  0.0322|         13.9596|                   94.4|
|SVDComplexity  |s(impact_mean) | 1.0004|  0.6884|          0.4561|                   94.4|
|SVDComplexity  |s(depth_m)     | 1.0000|  0.0757|          8.9453|                   94.4|
|SVDComplexity  |s(site)        | 0.9629|  0.0000|         73.7153|                   94.4|
|Modularity     |s(S)           | 1.0000|  0.0002|          1.9378|                   74.3|
|Modularity     |s(log_area)    | 1.0000|  0.0188|          0.7733|                   74.3|
|Modularity     |s(latitude)    | 1.0006|  0.0340|          0.6310|                   74.3|
|Modularity     |s(impact_mean) | 1.0000|  0.2208|          0.2101|                   74.3|
|Modularity     |s(depth_m)     | 1.0000|  0.8680|          0.0039|                   74.3|
|Modularity     |s(site)        | 0.9985|  0.0000|         96.4439|                   74.3|
|LD             |s(S)           | 1.0000|  0.8878|          0.0267|                   90.6|
|LD             |s(log_area)    | 1.0000|  0.6082|          0.3514|                   90.6|
|LD             |s(latitude)    | 1.0136|  0.3504|          1.2175|                   90.6|
|LD             |s(impact_mean) | 1.0001|  0.7834|          0.1011|                   90.6|
|LD             |s(depth_m)     | 1.0000|  0.1267|          3.1191|                   90.6|
|LD             |s(site)        | 0.9862|  0.0000|         95.1843|                   90.6|
|TLmean         |s(S)           | 1.0000|  0.1551|          0.0832|                   86.3|
|TLmean         |s(log_area)    | 1.0000|  0.0085|          0.2856|                   86.3|
|TLmean         |s(latitude)    | 1.0003|  0.0577|          0.1484|                   86.3|
|TLmean         |s(impact_mean) | 1.0000|  0.9238|          0.0004|                   86.3|
|TLmean         |s(depth_m)     | 1.0000|  0.8646|          0.0012|                   86.3|
|TLmean         |s(site)        | 0.9996|  0.0000|         99.4812|                   86.3|
|RankDeficiency |s(S)           | 1.0001|  0.0000|         16.8583|                   83.2|
|RankDeficiency |s(log_area)    | 1.0001|  0.0716|          1.6758|                   83.2|
|RankDeficiency |s(latitude)    | 1.0057|  0.4396|          0.3148|                   83.2|
|RankDeficiency |s(impact_mean) | 1.0002|  0.0023|          4.7819|                   83.2|
|RankDeficiency |s(depth_m)     | 1.0000|  0.7967|          0.0343|                   83.2|
|RankDeficiency |s(site)        | 0.9933|  0.0000|         76.3349|                   83.2|
|MEing_stable   |s(S)           | 1.0000|  0.1344|         12.3301|                   19.1|
|MEing_stable   |s(log_area)    | 1.0000|  0.2778|          6.4802|                   19.1|
|MEing_stable   |s(latitude)    | 1.0446|  0.0451|         23.1881|                   19.1|
|MEing_stable   |s(impact_mean) | 1.0002|  0.5918|          1.5806|                   19.1|
|MEing_stable   |s(depth_m)     | 1.0000|  0.4690|          2.8839|                   19.1|
|MEing_stable   |s(site)        | 0.9418|  0.0000|         53.5371|                   19.1|

Table S1. Results from generalized additive mixed models (GAMMs) testing the effects of number of species latitude, log-transformed area, mean human impact and depth on food web structural metrics. Reported are estimated degrees of freedom (edf), approximate p-values, approximate partial deviance contributions of each smooth term, and the total deviance explained by the model.

\newpage

![Relationship between dynamical stability, measured as the maximum eigenvalue of the community matrix, and latitude, log-transformed network area, and human impact across marine food webs. Negative values of the maximum eigenvalue indicate greater stability. Points show site-level means with 95% confidence intervals, colored by site, while dashed lines and shaded envelopes represent the fitted smooths (±95% confidence bands) from generalized additive mixed models (GAMMs). Model results indicate significant effects of latitude ($p = 0.0015$) and log area ($p = 0.0005$) on stability, while human impact showed no detectable effect ($p = 0.44$).](Figures/metric_Q_TL_MEing_vs_S_ci_gam.png)


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


### Mean Trophic Level (MTL)

The trophic level of species *i* is defined recursively as:

$$
TL_i = 1 + \frac{1}{k_i} \sum_{j} a_{ij} \, TL_j,
\quad \text{where} \quad k_i = \sum_{j} a_{ij}.
$$

Then, the mean trophic level for the food web is:

$$
MTL = \frac{1}{S} \sum_{i=1}^{S} TL_i.
$$


### Connectance (C)

Connectance quantifies the proportion of realized trophic interactions relative to all possible ones:

$$
C = \frac{L}{S^2}, \quad L = \sum_{i,j} a_{ij}.
$$


### Link Density (LD)

Link density measures the average number of trophic links per species:

$$
LD = \frac{L}{S}.
$$



###  Modularity (Q)

Modularity describes the extent to which the network is organized into modules with dense intra-module links and sparse inter-module links.
We detect modules using the **Infomap algorithm**, which minimizes the description length of a random walker’s trajectory to reveal community structure (Rosvall and Bergstrom, 2008).

Formally, modularity is reported as the fraction of links within modules relative to the whole network, given by the Infomap solution.

A widely used definition for modularity in ecological networks is:

$$
Q = \sum_{s=1}^{M} \left[ \frac{l_s}{L} - \left( \frac{d_s}{2L} \right)^2 \right]
$$

where:

* $M$ = total number of modules.
* $l_s$ = number of links within module *s*.
* $d_s$ = sum of degrees of nodes in module *s*.
* $L$ = total number of links in the network.

(This is the standard Newman-Girvan modularity.)

### References

* Infomap: Rosvall M., Bergstrom C.T. (2008). Maps of random walks on complex networks reveal community structure. *PNAS* 105(4): 1118–1123.


