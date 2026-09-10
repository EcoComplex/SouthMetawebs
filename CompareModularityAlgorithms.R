# CompareModularityAlgorithms.R
# Compara la distribucion de Modularidad (Q) obtenida con Infomap vs cluster_spinglass
# sobre las mismas ~1000 replicas del metaweb assembly model, por sitio.
#
# Motivacion: el IC de modularidad de Northern/Southern Scotia en las Figuras S8/S9
# es desproporcionadamente ancho (roza Q=0 en el limite inferior). Esta comparacion
# evalua si el patron se debe a inestabilidad de Infomap en redes grandes y densas,
# o si se repite con otro algoritmo de deteccion de comunidades (spinglass), en cuyo
# caso el patron es una propiedad real de las redes ensambladas y no un artefacto
# del metodo.
#
# Requiere: Data/simulations_metaWeb_metrics.rds (Infomap, pipeline principal)
#           Data/simulations_metaWeb_metrics_spinglass.rds (ya calculado, las 7 redes)

library(tidyverse)

keep_sites <- c("BurdwoodBank", "GulfSanJorge", "PotterCove",
                "Std BeagleChannel", "Std Weddell Sea",
                "Std NorthernScotia", "Std SouthernScotia")

nice_names <- c(
  "Std NorthernScotia" = "Northern Scotia",
  "Std SouthernScotia" = "Southern Scotia",
  "Std BeagleChannel"  = "Beagle Channel",
  "Std Weddell Sea"    = "Weddell Sea",
  "GulfSanJorge"       = "Gulf San Jorge",
  "PotterCove"         = "Potter Cove",
  "BurdwoodBank"       = "Burdwood Bank"
)

im <- readRDS("Data/simulations_metaWeb_metrics.rds") %>%
  filter(Metaweb %in% keep_sites) %>%
  mutate(Method = "Infomap")

sg <- readRDS("Data/simulations_metaWeb_metrics_spinglass.rds") %>%
  filter(Metaweb %in% keep_sites) %>%
  mutate(Method = "Spinglass")

both <- bind_rows(im, sg) %>%
  mutate(site = factor(nice_names[Metaweb], levels = unname(nice_names[keep_sites])))

# ---- Tabla resumen: media, sd, CV, cola cercana a cero, por metodo y sitio ----
summary_tbl <- both %>%
  group_by(site, Method) %>%
  summarise(
    mean            = mean(Modularity),
    sd              = sd(Modularity),
    CV              = sd / mean,
    frac_below_0.01 = mean(Modularity < 0.01),
    frac_exact_0    = mean(Modularity == 0),
    .groups = "drop"
  ) %>%
  arrange(desc(CV))

print(summary_tbl, n = Inf)
write_csv(summary_tbl, "Data/modularity_infomap_vs_spinglass_summary.csv")

# ---- Figura comparativa: densidad de Q por sitio, Infomap vs Spinglass ----
p <- ggplot(both, aes(x = Modularity, fill = Method)) +
  geom_density(alpha = 0.45, color = NA) +
  facet_wrap(~site, scales = "free", ncol = 4) +
  scale_fill_manual(values = c(Infomap = "#3B7EA1", Spinglass = "#D1615D")) +
  labs(x = "Modularity (Q)", y = "Density",
       title = "Distribucion de modularidad por replica: Infomap vs cluster_spinglass") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("Figures/modularity_infomap_vs_spinglass_comparison.png", p,
       width = 12, height = 7, dpi = 300, bg = "white")

cat("Listo: Figures/modularity_infomap_vs_spinglass_comparison.png y",
    "Data/modularity_infomap_vs_spinglass_summary.csv\n")
