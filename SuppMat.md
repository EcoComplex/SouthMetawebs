---
title: Supplementary Materials Gradientes latitudinales en la estructura y estabilidad de redes tróficas marinas
  del Atlántico Sudoccidental
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
|C              |s(S)           | 1.0000|  0.1321|          0.0440|                   89.4|
|C              |s(log_area)    | 1.0002|  0.0000|          0.7472|                   89.4|
|C              |s(latitude)    | 2.9996|  0.0000|         99.2037|                   89.4|
|C              |s(impact_mean) | 1.0001|  0.6048|          0.0052|                   89.4|
|SVDComplexity  |s(S)           | 1.0002|  0.8332|          0.1146|                   94.9|
|SVDComplexity  |s(log_area)    | 2.0623|  0.0778|         16.9324|                   94.9|
|SVDComplexity  |s(latitude)    | 1.9368|  0.0000|         71.8385|                   94.9|
|SVDComplexity  |s(impact_mean) | 1.0006|  0.0378|         11.1144|                   94.9|
|Modularity     |s(S)           | 1.0001|  0.0000|         52.4972|                   74.5|
|Modularity     |s(log_area)    | 1.9384|  0.0000|         20.2582|                   74.5|
|Modularity     |s(latitude)    | 2.0542|  0.0001|         17.6223|                   74.5|
|Modularity     |s(impact_mean) | 1.0051|  0.0002|          9.6223|                   74.5|
|LD             |s(S)           | 1.0023|  0.1170|         18.0881|                   90.6|
|LD             |s(log_area)    | 2.0435|  0.3747|         12.7398|                   90.6|
|LD             |s(latitude)    | 1.9499|  0.0632|         41.5093|                   90.6|
|LD             |s(impact_mean) | 1.0042|  0.0518|         27.6628|                   90.6|
|TLmean         |s(S)           | 1.0068|  0.0367|         24.4954|                   86.6|
|TLmean         |s(log_area)    | 1.0000|  0.0202|         30.2805|                   86.6|
|TLmean         |s(latitude)    | 1.8901|  0.0538|         34.0190|                   86.6|
|TLmean         |s(impact_mean) | 2.1028|  0.4346|         11.2051|                   86.6|
|RankDeficiency |s(S)           | 1.0002|  0.0000|         68.1208|                   84.3|
|RankDeficiency |s(log_area)    | 1.9171|  0.0000|         25.5057|                   84.3|
|RankDeficiency |s(latitude)    | 2.0776|  0.0503|          0.9691|                   84.3|
|RankDeficiency |s(impact_mean) | 1.0025|  0.0000|          5.4044|                   84.3|
|MEing_stable   |s(S)           | 1.0122|  0.0000|         56.2627|                   19.2|
|MEing_stable   |s(log_area)    | 1.2684|  0.0129|         16.5487|                   19.2|
|MEing_stable   |s(latitude)    | 1.0001|  0.0000|         26.6094|                   19.2|
|MEing_stable   |s(impact_mean) | 2.6914|  0.5414|          0.5792|                   19.2|

Table S1. Results from generalized additive models (GAMs) testing the effects of latitude, log-transformed area, and mean human impact on food web structural metrics. Reported are estimated degrees of freedom (edf), approximate p-values, approximate partial deviance contributions of each smooth term, and the total deviance explained by the model.

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


**1️⃣ Mean Trophic Level (MTL)**
The trophic level of species *i* is defined recursively as:

$$
TL_i = 1 + \frac{1}{k_i} \sum_{j} a_{ij} \, TL_j,
\quad \text{where} \quad k_i = \sum_{j} a_{ij}.
$$

Then, the mean trophic level for the food web is:

$$
MTL = \frac{1}{S} \sum_{i=1}^{S} TL_i.
$$

---

**2️⃣ Connectance (C)**

Connectance quantifies the proportion of realized trophic interactions relative to all possible ones:

$$
C = \frac{L}{S^2}, \quad L = \sum_{i,j} a_{ij}.
$$

---

**3️⃣ Link Density (LD)**
Link density measures the average number of trophic links per species:

$$
LD = \frac{L}{S}.
$$


---

**4️⃣ Modularity (Q)**

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

**References:**

* Infomap: Rosvall M., Bergstrom C.T. (2008). Maps of random walks on complex networks reveal community structure. *PNAS* 105(4): 1118–1123.


