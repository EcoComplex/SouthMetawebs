# Functions used by the analysis scripts of this repository.

#' Simulate the metaweb assembly model over a Latin hypercube of parameters
#'
#' @param no_parms Number of parameters to sample (2: m and a; 3: m, a and se).
#' @param no_sims Number of parameter combinations.
#' @param par_ranges Matrix with the range of each parameter in rows (m, a, se).
#' @param A Adjacency matrix of the metaweb.
#' @param model_type 0 discrete time (metaWebNetAssembly), 1 continuous time (metaWebNetAssemblyCT).
#' @param seed Optional random seed.
#'
#' @return A data frame with the parameters and the mean S, L and C of the last 100 time steps.
sim_metaWebAssembly_lhs <- function(no_parms, no_sims, par_ranges, A, model_type=0, seed = NULL) {
  stopifnot(no_parms %in% c(2,3))
  dimA <- nrow(A)
  if (!is.null(seed)) set.seed(seed)
  X <- randomLHS(no_sims, no_parms)

  if(no_parms==2) {
    m <- qunif(X[,1], par_ranges[1,1], par_ranges[1,2])
    e <- qunif(X[,2], par_ranges[2,1], par_ranges[2,2])
    se <- rep(par_ranges[3,1], times=no_sims)
  } else {
    m <- qunif(X[,1], par_ranges[1,1], par_ranges[1,2])
    e <- qunif(X[,2], par_ranges[2,1], par_ranges[2,2])
    se <- qunif(X[,3], par_ranges[3,1], par_ranges[3,2])
  }
  arguments <- data.frame(m = m, e = e, se = se)
  
  # Final time for the simulation
  tf <- 1000
  require(future.apply)
  num_cores <- future::availableCores()
  if (supportsMulticore()) {
    plan(multicore, workers = num_cores)
  } else {
    message("multicore not supported on this OS, falling back to multisession.")
    plan(multisession, workers = num_cores)
  }
  # Use future_lapply for parallelism
  sim_list <- future_lapply(1:no_sims, function(i) {
    # Set a reproducible seed for each worker if needed
    mm <- rep(arguments[i,1], times=dimA)
    aa <- rep(arguments[i,2], times=dimA)
    ee <- rep(arguments[i,3], times=dimA)

    if (model_type == 0) {
      AA <- metaWebNetAssembly(A, mm, aa, ee, tf)
    } else {
      AA <- metaWebNetAssemblyCT(A, mm, aa, ee, tf)
    }
    dfA <- data.frame(S=AA$S[(tf-100):tf], L=as.numeric(AA$L[(tf-100):tf]), T=c((tf-100):tf))
    data.frame(m=mm[1], a=aa[1], se=ee[1], S=mean(dfA$S), L=mean(dfA$L), C=mean(dfA$L)/(mean(dfA$S)*mean(dfA$S)))
  }, future.seed = TRUE)
  plan(sequential)
  sim <- do.call(rbind, sim_list)
  return(sim)
}

#' Fit metaweb assembly model from simulations using ABC and perform a goodness-of-fit test
#'
#' @param emp_S Empirical species richness (S) value for the network.
#' @param emp_C Empirical connectance (C) value for the network.
#' @param web_name Name of the network to fit.
#' @param sims Data frame with simulation results (should include columns S, C, m, a, se, model, model_type, Metaweb).
#' @param tol Tolerance for selecting simulations close to empirical S (default 0.1).
#' @param plot Logical; if TRUE, make and save a plot of the fitting (default TRUE).
#'
#' @return A data frame with ABC rejection results and goodness-of-fit p-values for each model.
#' @export
#'
#' @examples
#' CI_metaWebAssembly_model(emp_S = 20, emp_C = 0.15, web_name = "MyNetwork", sims = simMetaWebAssembly)
CI_metaWebAssembly_model <- function(emp_S, emp_C, web_name, sims, tol=0.1, plot=TRUE){
  # Filter simulations for the given web_name
  sims <- sims %>% filter(Metaweb == web_name)
  if(nrow(sims) == 0) {
    stop(sprintf("No simulations found for web_name: %s", web_name))
  }
  # Fit using the distance to S and C
  sims <- sims %>% mutate(distance = sqrt(((emp_S - S)/emp_S)^2 + ((emp_C - C)/emp_C)^2), min_dist = (distance == min(distance)), type = "Simulations") %>% arrange(distance)
  sel1 <- tibble(S = emp_S, C = emp_C, type = "Empirical")

    gof <- sims %>% group_modify(~ {
    gg <- summary(abc::gfit(target = sel1[1, 1:2], sumstat = .x[, c("S", "C")], nb.replicate = 1000))
    tibble(gof_pvalue = gg$pvalue)
  }) %>% mutate(Network = web_name, fit_type = "S-C")

  rej <- abc::abc(target = sel1[1, 1:2], param = sims[, c("m", "a", "se")], sumstat = sims[, c("S", "C")], tol = 0.05, method = "rejection")
  sel <- rej$unadj.values 
  sel <- inner_join(data.frame(sel),sims)
  if (plot) {
    print(
      ggplot(sel, aes(S, C, color = type)) +
        geom_point(alpha = 0.1) + theme_bw() + stat_ellipse() +
        geom_point(data = sel1, aes(S, C, color = type)) +
        scale_color_viridis_d(name = "", labels = c("Empirical", "Simulations"))
    )
    ggsave(paste0("Figures/Metaweb_fit_", web_name, "_byModel.png"), width = 8, height = 5, units = "in", dpi = 600)
  }
  return(list(sel=sel,gof=gof))
}

#' Simulate Food Webs with Metaweb Assembly Model and Compute Network Metrics
#'
#' This function performs multiple simulations of the metaweb assembly model (`metaWebNetAssemblyCT`)
#' using fitted parameters (`m`, `a`, `se`) sampled from a data frame. It calculates SVD entropy,
#' effective rank, modularity and topological indices for each simulation.
#'
#' @param Adj Adjacency matrix of the metaweb.
#' @param fitted_params A data.frame with columns `m`, `a`, `se`, and `Metaweb` for parameter sampling.
#' @param metaweb_name Character, the name of the metaweb to filter parameters (must match `Metaweb` column).
#' @param nsim Number of simulations to perform (default: 1000).
#' @param final_time Total time steps for the assembly simulation (default: 1000).
#' @param modularity_method Community detection algorithm: "infomap" (default, used
#'   in the analyses) or "spinglass" (used only to compare algorithms).
#'
#' @return A `tibble` with columns `Sim`, `Metaweb`, `S`, `L`, `C`, `LD`, `TLmean`,
#'   `Vulnerability`, `Generality`, `Entropy`, `Rank`, `Modularity`.
#' @importFrom igraph graph_from_adjacency_matrix components induced_subgraph vcount ecount
#' @importFrom multiweb calc_svd_entropy run_infomap
#' @importFrom dplyr bind_rows
#' @export
simulate_metaweb_metrics <- function(Adj, fitted_params, metaweb_name,
                                     nsim = 1000, final_time = 1000,
                                     modularity_method = c("infomap", "spinglass")) {
  modularity_method <- match.arg(modularity_method)

  require(igraph)
  require(dplyr)
  require(future.apply)
  
  param_set <- fitted_params %>% filter(Metaweb == metaweb_name)
  if (nrow(param_set) == 0) stop("No fitted parameters found for this metaweb.")
  
  num_cores <- future::availableCores()
  plan(multisession, workers = num_cores)

  sim_list <- future_lapply(1:nsim, function(i) {
    row <- param_set[sample(1:nrow(param_set), 1), ]
    m <- rep(row$m, nrow(Adj))
    a <- rep(row$a, nrow(Adj))
    se <- rep(row$se, nrow(Adj))
    
    sim <- metaWebNetAssemblyCT(Adj, m = m, e = a, se = se, time = final_time)
    g <- graph_from_adjacency_matrix(sim$A, mode = "directed")
    comps <- components(g)
    g <- induced_subgraph(g, which(comps$membership == which.max(comps$csize)))
    
    S <- vcount(g)
    L <- ecount(g)
    
    if (S < 3 || L == 0) {
      return(tibble(Sim = i, Metaweb = metaweb_name, S = S, L = L,
                    Entropy = NA_real_, Rank = NA_real_, Modularity = NA_real_))
    }
    
    entropy_vals <- tryCatch(multiweb::calc_svd_entropy(g), error = function(e) list(Entropy = NA_real_, Rank = NA_real_))
    mod <- tryCatch({
      comm <- if (modularity_method == "infomap") multiweb::run_infomap(g) else cluster_spinglass(g)
      modularity(comm)
    }, error = function(e) NA_real_)
    
    topo <- multiweb::calc_topological_indices(g)
    
    tibble(Sim = i, Metaweb = metaweb_name, S = S, L = L, C = topo$Connectance,
           LD = topo$LD,
           TLmean = topo$TLmean,
           Vulnerability = topo$Vulnerability,
           Generality = topo$Generality,
           Entropy = entropy_vals$Entropy, 
           Rank = entropy_vals$Rank, 
           Modularity = mod)
  }, future.seed = TRUE)
  
  plan(sequential)
  return(bind_rows(sim_list))
}

#' Faceted density plots of simulated metrics with observed values
#'
#' @param network_info Data frame with observed metrics (must include `site`)
#' @param simulated_metrics Data frame with simulated metrics (must include `Metaweb`)
#' @param metrics Character vector of metrics to plot
#'
#' @return A named list of ggplot objects (one per metric)
#' @export
#'
#' @import ggplot2 dplyr purrr
plot_density_faceted_by_site <- function(network_info,
                                         simulated_metrics,
                                         metrics = c("S", "C", "Modularity")) {
  library(dplyr)
  library(ggplot2)
  library(purrr)
  
  # Prepare simulated data
  sim_data <- simulated_metrics %>%
    rename(site = Metaweb) %>%
    filter(site %in% network_info$site)
  
  # Prepare observed data
  obs_data <- network_info %>%
    select(site, everything()) %>%
    filter(site %in% sim_data$site)
  
  # Determine metrics automatically if needed
  if (is.null(metrics)) {
    numeric_cols <- map_lgl(obs_data, is.numeric)
    candidate_metrics <- names(obs_data)[numeric_cols]
    metrics <- intersect(candidate_metrics, names(sim_data))
    metrics <- setdiff(metrics, c("latitude", "depth_m", "area_km2"))
  }
  
  # Function to build one faceted plot per metric
  plot_one_metric <- function(metric_name) {
    
    obs_metric <- obs_data %>%
      select(site, value = all_of(metric_name))
    
    ggplot(sim_data, aes(x = .data[[metric_name]])) +
      geom_density(fill = "skyblue", alpha = 0.5, linewidth = 0.5) +
      geom_vline(data = obs_metric,
                 aes(xintercept = value),
                 color = "red",
                 linetype = "dashed",
                 linewidth = 0.9) +
      facet_wrap(~ site, scales = "free") +
      labs(x = metric_name,
           y = "Density") +
      theme_bw(base_size = 13) +
      theme(strip.background = element_rect(fill = "gray90"),
            strip.text = element_text(face = "bold"))
  }
  
  # Generate plots
  plots <- map(metrics, plot_one_metric)
  names(plots) <- metrics
  
  plots
}

#' Parallel QSS Computation Across Networks and Self-Damping Values
#'
#' This function evaluates QSS and stability metrics across multiple networks
#' and self-damping values in parallel, ensuring no nested parallelism.
#'
#' @param networks List of igraph networks (same length as sites).
#' @param sites Site names matching networks.
#' @param selfDamping_range Sequence of selfDamping values to try.
#' @param nsim Number of simulations per run (default: 1000).
#' @param negative Max magnitude of negative interactions.
#' @param positive Max magnitude of positive interactions.
#'
#' @return A data frame with QSS, MEing_stable and selfDamping for each site and
#'   self-damping value with QSS > 0.
#' @export
calculate_network_stability_parallel <- function(networks,
                                                 sites,
                                                 selfDamping_range = seq(-1, -10, by = -0.5),
                                                 nsim = 1000,
                                                 negative = -10,
                                                 positive = 1) {
  stopifnot(length(networks) == length(sites))
  
  library(dplyr)
  library(future.apply)
  library(purrr)
  
  # Expand combinations: each row = 1 task
  combos <- expand.grid(
    idx = seq_along(networks),
    selfDamping = selfDamping_range,
    KEEP.OUT.ATTRS = FALSE
  )
  
  combos_list <- split(combos, seq(nrow(combos)))
  
  # Parallel computation across combinations
  results_raw <- future_lapply(combos_list, function(row) {
    i <- row$idx
    sd <- row$selfDamping
    g <- networks[[i]]
    site <- sites[[i]]
    
    eigen_vals <- tryCatch(
      calc_QSS(
        g,
        nsim = nsim,
        ncores = 0,  # ensure no nested parallelism
        negative = negative,
        positive = positive,
        selfDamping = sd,
        returnRaw = TRUE
      ),
      error = function(e) NA_real_
    )
    
    if (all(is.na(eigen_vals))) return(NULL)
    
    qss <- mean(eigen_vals < 0, na.rm = TRUE)
    
    if (qss > 0) {
      MEing_stable <- mean(eigen_vals[eigen_vals < 0], na.rm = TRUE)
      return(tibble(site = site, QSS = qss, MEing_stable = MEing_stable, selfDamping = sd))
    } else {
      return(NULL)
    }
  })
  
  bind_rows(results_raw)
}

#' Simulate QSS and Stability at Given selfDamping Values per Network
#'
#' This function computes QSS and the mean of negative real parts of maximal eigenvalues
#' at fixed self-damping values for a list of networks.
#'
#' @param networks A list of igraph food web networks.
#' @param sites A character vector with names of the networks (same length as `networks`).
#' @param selfDamping_values A numeric vector with selfDamping value for each network.
#' @param nsim Number of simulations per network (default: 1000).
#' @param negative Maximum magnitude of negative interaction strength (default: -10).
#' @param positive Maximum magnitude of positive interaction strength (default: 1).
#'
#' @return A data frame with `site`, `QSS`, `MEing_stable`, and `selfDamping`.
#' @export
simulate_qss_fixed_selfDamping <- function(networks,
                                           sites,
                                           selfDamping_values,
                                           nsim = 1000,
                                           negative = -10,
                                           positive = 1) {
  stopifnot(length(networks) == length(sites), length(sites) == length(selfDamping_values))
  
  library(dplyr)
  library(future.apply)
  
  results <- future_lapply(seq_along(networks), function(i) {
    g <- networks[[i]]
    site <- sites[[i]]
    sd <- selfDamping_values[[i]]
    
    eigen_vals <- tryCatch(
      calc_QSS(
        g,
        nsim = nsim,
        ncores = 0,
        negative = negative,
        positive = positive,
        selfDamping = sd,
        returnRaw = TRUE
      ),
      error = function(e) NA_real_
    )
    
    if (all(is.na(eigen_vals))) return(NULL)
    
    qss <- mean(eigen_vals < 0, na.rm = TRUE)
    meing <- eigen_vals[eigen_vals<0]
    
    tibble(site = site, QSS = qss, MEing_stable = meing, selfDamping = sd)
  }, future.seed = TRUE)
  
  return(bind_rows(results))
}
