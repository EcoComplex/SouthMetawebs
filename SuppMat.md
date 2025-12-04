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

|Metric         |Term           |    edf| p-value| PartialDeviance| TotalDevianceExplained|
|:--------------|:--------------|------:|-------:|---------------:|----------------------:|
|C              |s(S)           | 1.0003|  0.4104|          2.2166|                   89.3|
|C              |s(log_area)    | 1.0003|  0.8398|          0.1341|                   89.3|
|C              |s(latitude)    | 1.0480|  0.3401|          3.4578|                   89.3|
|C              |s(impact_mean) | 1.0006|  0.7787|          0.2588|                   89.3|
|C              |s(depth_m)     | 1.0001|  0.0021|         30.8584|                   89.3|
|C              |s(site)        | 0.9507|  0.0000|         63.0743|                   89.3|
|SVDComplexity  |s(S)           | 1.0002|  0.6333|          0.6153|                   93.8|
|SVDComplexity  |s(log_area)    | 1.0002|  0.4130|          1.8131|                   93.8|
|SVDComplexity  |s(latitude)    | 1.0337|  0.0335|         13.0496|                   93.8|
|SVDComplexity  |s(impact_mean) | 1.0004|  0.7195|          0.3486|                   93.8|
|SVDComplexity  |s(depth_m)     | 1.0000|  0.0789|          8.3490|                   93.8|
|SVDComplexity  |s(site)        | 0.9655|  0.0000|         75.8244|                   93.8|
|Modularity     |s(S)           | 1.0000|  0.0002|          1.9771|                   73.8|
|Modularity     |s(log_area)    | 1.0000|  0.0221|          0.7333|                   73.8|
|Modularity     |s(latitude)    | 1.0006|  0.0306|          0.6549|                   73.8|
|Modularity     |s(impact_mean) | 1.0000|  0.2107|          0.2193|                   73.8|
|Modularity     |s(depth_m)     | 1.0000|  0.8865|          0.0029|                   73.8|
|Modularity     |s(site)        | 0.9986|  0.0000|         96.4126|                   73.8|
|LD             |s(S)           | 1.0000|  0.8852|          0.0264|                   90.7|
|LD             |s(log_area)    | 1.0000|  0.6393|          0.2778|                   90.7|
|LD             |s(latitude)    | 1.0128|  0.3455|          1.1715|                   90.7|
|LD             |s(impact_mean) | 1.0001|  0.8092|          0.0738|                   90.7|
|LD             |s(depth_m)     | 1.0000|  0.1325|          2.8606|                   90.7|
|LD             |s(site)        | 0.9869|  0.0000|         95.5899|                   90.7|
|TLmean         |s(S)           | 1.0000|  0.1558|          0.0797|                   86.2|
|TLmean         |s(log_area)    | 1.0000|  0.0098|          0.2641|                   86.2|
|TLmean         |s(latitude)    | 1.0003|  0.0556|          0.1452|                   86.2|
|TLmean         |s(impact_mean) | 1.0000|  0.9132|          0.0005|                   86.2|
|TLmean         |s(depth_m)     | 1.0000|  0.8530|          0.0014|                   86.2|
|TLmean         |s(site)        | 0.9996|  0.0000|         99.5092|                   86.2|
|RankDeficiency |s(S)           | 1.0001|  0.0000|         16.9191|                   83.5|
|RankDeficiency |s(log_area)    | 1.0001|  0.0783|          1.5955|                   83.5|
|RankDeficiency |s(latitude)    | 1.0056|  0.4253|          0.3338|                   83.5|
|RankDeficiency |s(impact_mean) | 1.0002|  0.0023|          4.7862|                   83.5|
|RankDeficiency |s(depth_m)     | 1.0000|  0.7989|          0.0334|                   83.5|
|RankDeficiency |s(site)        | 0.9933|  0.0000|         76.3319|                   83.5|
|MEing_stable   |s(S)           | 1.0000|  0.1318|         12.3887|                   19.2|
|MEing_stable   |s(log_area)    | 1.0001|  0.2921|          6.0545|                   19.2|
|MEing_stable   |s(latitude)    | 1.0436|  0.0434|         23.3134|                   19.2|
|MEing_stable   |s(impact_mean) | 1.0002|  0.6029|          1.4748|                   19.2|
|MEing_stable   |s(depth_m)     | 1.0000|  0.4711|          2.8321|                   19.2|
|MEing_stable   |s(site)        | 0.9428|  0.0000|         53.9365|                   19.2|

Table: Results from generalized additive mixed models (GAMMs) testing the effects of number of species latitude, log-transformed area, mean human impact and depth on food web structural metrics. Reported are estimated degrees of freedom (edf), approximate p-values, approximate partial deviance contributions of each smooth term, and the total deviance explained by the model.

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


