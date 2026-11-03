# FoodChainLength.R
# For the seven observed food webs, computes prey-averaged trophic levels and
# compares mean trophic level (our response) with food chain length
# (maximum TL and related measures). Run from the project root
# (open SouthMetawebs.Rproj and source this file).
#
# Reads:  Data/network_info.rds (igraph objects and TLmean from multiweb)
# Writes: Data/food_chain_length.csv

library(tidyverse)
library(igraph)

network_info <- readRDS("Data/network_info.rds")

# Trophic levels with NetIndices::TrophInd, the same function multiweb uses
# for TLmean, so TLmax and TLmean share one definition (TL = 1 + mean TL of
# prey; Flow[i, j] = 1 if j eats i; self-feeding ignored). TLmax equals
# the TLmax column returned by multiweb::calc_topological_indices.
# Species with no path from a basal species (closed loops of consumers that
# only feed on each other) are reported as a data check; they are kept, as in
# the TLmean used in the paper.
trophic_levels <- function(g, site = "") {
  # Same matrix as multiweb::calc_topological_indices (not simplified)
  A <- as_adjacency_matrix(g, sparse = FALSE)               # A[prey, predator]
  tl <- setNames(NetIndices::TrophInd(A)$TL, V(g)$name)
  if (any(A > 1)) message(site, ": ", sum(A > 1), " duplicated links")
  gs <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  basal <- V(gs)[degree(gs, mode = "in") == 0]
  reach <- unique(unlist(lapply(basal, function(b)
    subcomponent(gs, b, mode = "out")$name)))
  orphan <- setdiff(V(gs)$name, reach)
  if (length(orphan) > 0)
    message(site, ": no path from basal species: ",
            paste(sprintf("%s (TL = %.2f)", orphan, tl[orphan]), collapse = ", "))
  tl
}

# Shortest TL: 1 + shortest path length from any basal species.
# Short-weighted TL (Williams & Martinez 2004): mean of shortest and
# prey-averaged TL; the most accurate binary estimate of flow-based TL, and
# the closest analogue of isotope-based food chain length (trophic position
# of the top predator). Species with no path from a basal species get NA.
short_weighted_TL <- function(g, tl_pa) {
  gs <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  basal <- V(gs)[degree(gs, mode = "in") == 0]
  d <- distances(gs, v = basal, to = V(gs), mode = "out")
  tl_short <- 1 + apply(d, 2, min)
  tl_short[is.infinite(tl_short)] <- NA
  names(tl_short) <- V(gs)$name
  (tl_short + tl_pa[names(tl_short)]) / 2
}

fcl <- network_info %>%
  mutate(
    tl     = map2(igraph, site, trophic_levels),           # prey-averaged
    sw     = map2(igraph, tl, short_weighted_TL),          # short-weighted
    top    = map(igraph, function(g) {
      g <- simplify(g)
      V(g)$name[degree(g, mode = "out") == 0 & degree(g, mode = "in") > 0]
    }),
    TLmean_check = map_dbl(tl, mean),
    TLmax   = map_dbl(tl, max),                            # prey-averaged max
    TLtop   = map2_dbl(tl, top, ~ mean(.x[intersect(.y, names(.x))])),
    swTLmean = map_dbl(sw, mean, na.rm = TRUE),
    swTLmax  = map_dbl(sw, max, na.rm = TRUE),             # food chain length
    swTLtop  = map2_dbl(sw, top, ~ mean(.x[intersect(.y, names(.x))], na.rm = TRUE)),
    propBasal = map_dbl(tl, ~ mean(abs(.x - 1) < 1e-8))
  ) %>%
  transmute(site, name, latitude, area_km2, S, TLmean, TLmean_check,
            TLmax, TLtop, swTLmean, swTLmax, swTLtop, propBasal)

print(fcl %>% arrange(latitude), width = Inf)

# Should be 0: same function as multiweb's TLmean.
cat("\nMax |TLmean - TLmean_check|:",
    round(max(abs(fcl$TLmean - fcl$TLmean_check)), 3), "\n")

# How comparable is mean TL (our response) with food chain length?
cat("\nCorrelation of TLmean with other vertical-structure measures (n = 7)\n")
for (v in c("TLmax", "TLtop", "swTLmean", "swTLmax", "swTLtop")) {
  p <- cor(fcl$TLmean, fcl[[v]], method = "pearson")
  s <- cor(fcl$TLmean, fcl[[v]], method = "spearman")
  cat(sprintf("  %-8s Pearson r = %5.2f   Spearman rho = %5.2f\n", v, p, s))
}

write_csv(fcl, "Data/food_chain_length.csv")
