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


|metric         |term        | estimate| p_value|
|:--------------|:-----------|--------:|-------:|
|C             |latitude    |  -0.0004|       0|
|C             |log_area    |   0.0195|       0|
|C             |impact_mean |   0.1493|       0|
|SVDComplexity |latitude    |   0.0019|       0|
|SVDComplexity |log_area    |   0.0032|       0|
|SVDComplexity |impact_mean |   0.0063|       0|
|Modularity    |latitude    |   0.0013|       0|
|Modularity    |log_area    |   0.0059|       0|
|Modularity    |impact_mean |  -0.1134|       0|
|LD             |latitude    |   0.3913|   0e+00|
|LD             |log_area    |   1.1366|   0e+00|
|LD             |impact_mean |   0.7398|   0e+00|
|TLmean         |latitude    |   0.0133|   0e+00|
|TLmean         |log_area    |   0.0947|   0e+00|
|TLmean         |impact_mean |   0.0652|   3e-04|
|RankDeficiency |latitude    |   0.0083|   0e+00|
|RankDeficiency |log_area    |   0.0122|   0e+00|
|RankDeficiency |impact_mean |  -0.1394|   0e+00|
|MEing_stable   |latitude    |  -0.0018|  0.0000|
|MEing_stable   |log_area    |  -0.0009|  0.2227|
|MEing_stable   |impact_mean |  -0.0288|  0.0460|

Table S1. Results from linear models testing the effects of latitude, log-transformed area, and mean human impact on food web structural metrics. Significant p-values (<0.05) are highlighted in bold.

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


