## compare_networks_refined.R
## Compara redes "original" con "standardized" usando los nombres actualizados
## del archivo Links_data_2024_06_24_Time_04_31_40.xlsx.
##
## La actualización de nombres se aplica link-a-link (IDlink_original_network),
## lo que maneja correctamente splits (una especie vieja -> varias nuevas).
##
## Normalización de nombres:
##   - Puntos -> espacios (CSV BeagleChannel usa "Genus.species")
##   - Espacios múltiples colapsados, trim
##
## Outputs:
##   {net}_links_original_updated.csv   — original con nombres actualizados
##   {net}_links_standardized_updated.csv — standardized normalizado
##   network_differences_refined.csv    — diffs link a link por sitio
##   network_metrics_summary_refined.csv — métricas S/L/C/LD por versión y sitio

library(readxl)
library(dplyr)
library(igraph)

# ── Rutas ──────────────────────────────────────────────────────────────────
DATA_DIR   <- "Data"
LINKS_DATA <- file.path(DATA_DIR, "Links_data_2024_06_24_Time_04_31_40.xlsx")

files <- list(
  BeagleChannel = list(
    orig = file.path(DATA_DIR, "BeagleChannel_links_original_2024_06_24_Time_04_30_52.xlsx"),
    std  = file.path(DATA_DIR, "BeagleChannel_links_standarized.csv")
  ),
  WeddellSea = list(
    orig = file.path(DATA_DIR, "WeddellSea_links_original_2024_06_24_Time_04_31_27.xlsx"),
    std  = file.path(DATA_DIR, "WeddellSea_links_standarized.csv")
  ),
  NorthernScotia = list(
    orig = file.path(DATA_DIR, "NorthernScotia_links_original_2024_06_24_Time_04_31_01.xlsx"),
    std  = file.path(DATA_DIR, "NorthernScotia_links_standardized.csv")
  ),
  SouthernScotia = list(
    orig = file.path(DATA_DIR, "SouthernScotia_links_original_2024_06_24_Time_04_31_14.xlsx"),
    std  = file.path(DATA_DIR, "SouthernScotia_links_standardized.csv")
  )
)

# ── Funciones ──────────────────────────────────────────────────────────────

# Normaliza nombres: puntos->espacios, trim, espacios múltiples
normalize_name <- function(x) {
  x <- trimws(gsub("\\.", " ", x))
  trimws(gsub("\\s+", " ", x))
}

read_links <- function(path) {
  df <- if (grepl("\\.csv$", path)) read.csv(path, stringsAsFactors = FALSE) \
        else readxl::read_excel(path)
  df |>
    dplyr::select(resource, consumer) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), normalize_name))
}

basic_metrics <- function(df) {
  g <- igraph::graph_from_data_frame(df, directed = TRUE)
  S <- igraph::vcount(g); L <- igraph::ecount(g)
  data.frame(S = S, L = L,
             Connectance = round(L / S^2, 4),
             LinkDensity = round(L / S, 2))
}

# ── Cargar tabla de actualizaciones de nombres ─────────────────────────────
links_data <- readxl::read_excel(LINKS_DATA)

# ── Procesar cada red ──────────────────────────────────────────────────────
results     <- list()
summary_tbl <- list()
all_diffs   <- list()

for (net in names(files)) {
  cat("\n", strrep("=", 60), "\n")
  cat(" ", net, "\n")
  cat(strrep("=", 60), "\n")

  orig <- read_links(files[[net]]$orig)
  std  <- read_links(files[[net]]$std)

  # Aplicar renombres link-a-link desde Links_data
  net_links   <- links_data |> dplyr::filter(network == net)
  orig_updated <- orig
  rownames(orig_updated) <- seq_len(nrow(orig_updated))  # IDlink empieza en 1

  # Resource renames
  r_renames <- net_links |>
    dplyr::filter(!is.na(resource_original_name)) |>
    dplyr::select(IDlink_original_network, resource_original_name)
  for (i in seq_len(nrow(r_renames))) {
    lid <- r_renames$IDlink_original_network[i]
    if (lid <= nrow(orig_updated))
      orig_updated$resource[lid] <- normalize_name(r_renames$resource_original_name[i])
  }

  # Consumer renames
  c_renames <- net_links |>
    dplyr::filter(!is.na(consumer_original_name)) |>
    dplyr::select(IDlink_original_network, consumer_original_name)
  for (i in seq_len(nrow(c_renames))) {
    lid <- c_renames$IDlink_original_network[i]
    if (lid <= nrow(orig_updated))
      orig_updated$consumer[lid] <- normalize_name(c_renames$consumer_original_name[i])
  }

  # Guardar archivos actualizados
  write.csv(orig_updated, paste0(net, "_links_original_updated.csv"),    row.names = FALSE)
  write.csv(std,          paste0(net, "_links_standardized_updated.csv"), row.names = FALSE)

  # Comparar
  o_set   <- paste(orig_updated$resource, "->", orig_updated$consumer)
  s_set   <- paste(std$resource,           "->", std$consumer)
  added   <- setdiff(s_set, o_set)
  removed <- setdiff(o_set, s_set)
  common  <- intersect(o_set, s_set)

  o_sp <- union(orig_updated$resource, orig_updated$consumer)
  s_sp <- union(std$resource,          std$consumer)
  sp_add <- setdiff(s_sp, o_sp)
  sp_rem <- setdiff(o_sp, s_sp)

  cat(sprintf(
    "\n  Links: orig_upd=%d  std=%d  comunes=%d  +%d  -%d\n",
    length(o_set), length(s_set), length(common), length(added), length(removed)
  ))
  cat(sprintf(
    "  Spp:   orig_upd=%d  std=%d  nuevas=%d  eliminadas=%d\n",
    length(o_sp), length(s_sp), length(sp_add), length(sp_rem)
  ))
  if (length(sp_add) > 0) cat("  Spp nuevas:", paste(sort(sp_add), collapse = ", "), "\n")
  if (length(sp_rem) > 0) cat("  Spp elim.: ", paste(sort(sp_rem), collapse = ", "), "\n")

  # Métricas
  mo <- basic_metrics(orig_updated); ms <- basic_metrics(std)
  cat("\n  Métricas:\n")
  print(rbind(cbind(version = "original_updated", mo),
              cbind(version = "standardized",     ms)), row.names = FALSE)

  summary_tbl[[net]] <- rbind(
    cbind(site = net, version = "original_updated", mo),
    cbind(site = net, version = "standardized",     ms)
  )

  parse_link <- function(links, change) {
    if (length(links) == 0) return(data.frame())
    parts <- strsplit(links, " -> ")
    data.frame(
      site     = net,
      change   = change,
      resource = sapply(parts, `[`, 1),
      consumer = sapply(parts, `[`, 2),
      stringsAsFactors = FALSE
    )
  }
  all_diffs[[net]] <- rbind(parse_link(added, "added"), parse_link(removed, "removed"))

  results[[net]] <- list(orig_updated = orig_updated, std = std,
                         added = added, removed = removed,
                         sp_add = sp_add, sp_rem = sp_rem)
}

# ── Exportar tablas globales ───────────────────────────────────────────────
summary_out <- do.call(rbind, summary_tbl)
diffs_out   <- do.call(rbind, all_diffs)

write.csv(summary_out, "network_metrics_summary_refined.csv",  row.names = FALSE)
write.csv(diffs_out,   "network_differences_refined.csv",      row.names = FALSE)

cat("\n\nArchivos generados:\n")
for (net in names(files)) {
  cat(sprintf("  %s_links_original_updated.csv\n", net))
  cat(sprintf("  %s_links_standardized_updated.csv\n", net))
}
cat("  network_differences_refined.csv\n")
cat("  network_metrics_summary_refined.csv\n")
