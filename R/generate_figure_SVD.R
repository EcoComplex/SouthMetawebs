require(tidyverse)
require(igraph)
require(multiweb)

g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 4,3 -+ 5,1 -+ 5 -+7, 2 -+ 6, 3 -+ 6 -+7, simplify = FALSE)
#g <-   graph_from_literal( 1 -+ 4 -+ 7,2 -+ 4,3 -+ 4,3 -+ 5,1 -+ 5 -+7, 2 -+ 6, 3 -+ 6 -+7, simplify = FALSE)

plot_troph_level_ggplot(g,modules = TRUE)
(a <- as_adjacency_matrix(g)[order(rownames(as_adjacency_matrix(g))), order(colnames(as_adjacency_matrix(g)))])
calc_topological_indices(g)
calc_svd_entropy(g, threshold_factor = 1e-15)
g1 <- generate_er_basal(g,nsim=1)[[1]]
as_adjacency_matrix(g1)

calc_svd_entropy(g1)
calc_topological_indices(g1)

n_basal <- function(g) {
  sum(degree(g, mode = "in") == 0)
}

rank_deficiency_from_g <- function(g) {
  res <- calc_svd_entropy(g)
  tibble(
    rank_def = 1 - (res$Rank / vcount(g)),
    complexity = 1 - res$Entropy
    )
}
set.seed(123)
gn  <- generate_niche(20,0.1,nsim=1) 
plot_troph_level_ggplot(gn,modules = TRUE)

as_adjacency_matrix(gn)

calc_svd_entropy(gn, threshold_factor = 1e-15)
n_basal(gn)
rank_deficiency_from_g(gn)
gn
g1 <- generate_er_basal(gn,nsim=1)[[1]]
as_adjacency_matrix(g1)
calc_svd_entropy(g1, threshold_factor = 1e-15)
n_basal(g1)
rank_deficiency_from_g(g1)


#
# Dependencia de la deficiencia de rango/complejidad con el número de especies basales
#
gn <- generate_niche(100, 0.05, nsim = 1000)

df <- do.call(rbind, lapply(gn, function(g) {
  res <- calc_svd_entropy(g)
  
  tibble(
    basal   = n_basal(g),
    rank    = res$Rank,
    entropy = res$Entropy,
    rank_def = 1 - (res$Rank / vcount(g)),
    complexity = 1 - entropy
  )
}))

library(ggplot2)

ggplot(df, aes(x = basal, y = rank_def)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_minimal() +
  labs(
    x = "Number of basal species",
    y = "Rank deficiency",
    title = "Rank deficiency vs basal species"
  )

ggplot(df, aes(x = basal, y = complexity)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_minimal() +
  labs(
    x = "Number of basal species",
    y = "Rank deficiency",
    title = "Rank deficiency vs basal species"
  )
