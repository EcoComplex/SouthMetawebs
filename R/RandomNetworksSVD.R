# RandomNetworksSVD.R
# Supplementary figure Figures/random_vs_empirical_SVDE.png: SVD entropy of the
# seven observed food webs compared with directed random networks (Erdos-Renyi)
# with numbers of species and connectance in the observed range.
#
# Reads: Data/network_info.rds
# Requires: tidyverse, igraph, multiweb, future, furrr

library(tidyverse)
library(igraph)
library(multiweb)
library(future)
library(furrr)

network_info <- readRDS("Data/network_info.rds")

# Observed range of S and C
S_min <- min(network_info$S)
S_max <- max(network_info$S)
C_min <- min(network_info$C)
C_max <- max(network_info$C)

# SVD entropy and rank deficiency of one random directed network
calc_sim_metrics <- function(S, C) {
  g <- sample_gnm(S, floor(C * S^2), directed = TRUE)
  A <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
  svde <- calc_svd_entropy(A)
  tibble(
    Entropy = svde$Entropy,
    rankDeficiency = 1 - (svde$Rank / vcount(g))
  )
}

plan(multisession, workers = future::availableCores())

set.seed(123)
sim_results <- tibble(
  S = sample(S_min:S_max, 100, replace = TRUE),
  C = runif(100, C_min, C_max)
) %>%
  mutate(metrics = future_map2(S, C, calc_sim_metrics, .progress = TRUE,
                               .options = furrr_options(seed = TRUE))) %>%
  unnest(metrics)

plan(sequential)

sim_results %>%
  ggplot(aes(x = S, y = Entropy, color = S)) +
  geom_point(alpha = 0.6) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(y = "SVD Entropy", color = "Species (S)") +
  geom_point(data = network_info, aes(x = S, y = Entropy), size = 3, shape = 17) +
  guides(color = "none")

ggsave("Figures/random_vs_empirical_SVDE.png", height = 6, width = 10, bg = "white")
