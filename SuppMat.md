# Supplementary Materials

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


