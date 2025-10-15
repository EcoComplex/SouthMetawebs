#' Calculate motif counts for observed network and CI for erdos-renyi random networks and Z-scores 
#'
#' @param red igraph network object
#' @param nsim number of simulation to calculate random networks with the same nodes and links
#'
#' @return data.frame with all the results
#' @export
#'
#' @examples
calc_motif_random <- function(red, nsim=1000)
{
    Size <- vcount(red)
    Links <- ecount(red)
    
    redes.r <- lapply(1:nsim, function (x) {
      e <- erdos.renyi.game(Size, Links, type="gnm",directed = TRUE)
      while(components(e)$no>1)
        e <- erdos.renyi.game(Size, Links, type="gnm",directed = TRUE)
      return(e) }
    )
    
    ind <- data.frame()
    require(doFuture)
    registerDoFuture()
    plan(multisession)
    
    ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% 
    {
      mot <- triad_census(redes.r[[i]])
      mot[4] # Exploitative competition
      mot[5] # Apparent competition
      mot[6] # Tri-trophic chain
      mot[9] # Omnivory

     data.frame(explComp=mot[4],apprComp=mot[5],triTroph=mot[6],omnivory=mot[9])
    }
    plan(sequential)
    # 99% confidence interval
    #
    
    qEC <- quantile(ind$explComp,c(0.005,0.995))
    qAC <- quantile(ind$apprComp,c(0.005,0.995))
    qTT <- quantile(ind$triTroph,c(0.005,0.995))
    qOM <- quantile(ind$omnivory,c(0.005,0.995))
    
    # Calculate motif for the original network
    obs <- triad_census(red)
    
    zEC <- (obs[4] - mean(ind$explComp))/sd(ind$explComp)
    zAC <- (obs[5] - mean(ind$apprComp))/sd(ind$apprComp)
    zTT <- (obs[6] - mean(ind$triTroph))/sd(ind$triTroph)
    zOM <- (obs[9] - mean(ind$omnivory))/sd(ind$omnivory)
    
    return(tibble(explComp=obs[4],apprComp=obs[5],triTroph=obs[6],omnivory=obs[9],zEC=zEC,zAC=zAC,zTT=zTT,zOM=zOM,EClow=qEC[1],EChigh=qEC[2],AClow=qAC[1],AChigh=qAC[2],TTlow=qTT[1],TThigh=qTT[2],OMlow=qOM[1],OMhigh=qOM[2]))         
}    

# Calculation of the clustering coefficients and average path for random network simulations
#
#
calc_modularity_random<- function(red, nsim=1000){
  
  t <- calc_topological_indices(red)
    
  redes.r <- lapply(1:nsim, function (x) {
    e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
    while(components(e)$no>1)
      e <- erdos.renyi.game(t$Size, t$Links, type="gnm",directed = TRUE)
    
    return(e) }
    )

  ind <- data.frame()
  require(doFuture)
  registerDoFuture()
  plan(multisession)
  
#   require(doParallel)
#   cn <-detectCores()
# #  cl <- makeCluster(cn,outfile="foreach.log") # Logfile to debug 
#   cl <- makeCluster(cn)
#   registerDoParallel(cl)
  
  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages='igraph') %dopar% 
    {
    m<-calc_modularity_unconnected(redes.r[[i]])
    modl <- m$modularity
    ngrp <- max(m$membership) 
    clus.coef <- transitivity(redes.r[[i]], type="Global")
    cha.path  <- average.path.length(redes.r[[i]])
    data.frame(modularity=modl,ngroups=ngrp,clus.coef=clus.coef,cha.path=cha.path)
  }
  # stopCluster(cl)
  plan(sequential)
  ind <- ind %>% mutate(gamma=t$Clustering/clus.coef,lambda=t$PathLength/cha.path,SWness=gamma/lambda)
  # 99% confidence interval
  #
  qSW <- quantile(ind$SWness,c(0.005,0.995),na.rm = TRUE)
  qmo <- quantile(ind$modularity,c(0.005,0.995))
  qgr <- quantile(ind$ngroups,c(0.005,0.995))
  mcc <- mean(ind$clus.coef)
  mcp <- mean(ind$cha.path)
  mmo <- mean(ind$modularity)
  mgr <- mean(ind$ngroups)
  mSW <- mean(t$Clustering/mcc*mcp/t$PathLength)
  mCI <- 1+(qSW[2]-qSW[1])/2  
  return(list(su=tibble(rndCC=mcc,rndCP=mcp,rndMO=mmo,rndGR=mgr,SWness=mSW,SWnessCI=mCI,MOlow=qmo[1],MOhigh=qmo[2],
                    GRlow=qgr[1],GRhigh=qgr[2]), sim=ind))         
}



#' Title Plot net assembly model S and L average by a moving window to check if equilibrium is reached
#'
#' @param AA output of a net assembly model
#' @param timeW time window used
#' @param fname file name to save the plot 
#' @param emp data.frame with empirical values of S an C
#'
#' @return
#' @export
#'
#' @examples
plot_NetAssemblyModel_eqw <- function(AA,timeW,fname=NULL,emp=NULL){

  df <- data.frame(S=AA$S,L=as.numeric(AA$L),T=c(1:tf)) %>% mutate(C=L/(S*S))
  grandS <- mean(df$S[timeW:nrow(df)])
  grandC <- mean(df$C[timeW:nrow(df)])
  
  df$gr <- rep(1:(nrow(df)/timeW), each = timeW)
  df <- df %>% filter(!is.nan(C)) %>% group_by(gr) %>% summarise(mS=mean(S),sdS=sd(S), mL=mean(L), sdL=sd(L),time=max(T),mC=mean(C),sdC=sd(C))
  if(is.null(fname)){
    print(g1 <- ggplot(df,aes(y=mS,x=time,colour=time))+ theme_bw() + geom_point() + geom_errorbar(aes(ymin=mS-sdS,ymax=mS+sdS)) + scale_color_distiller(palette = "RdYlGn",guide=FALSE)+ geom_hline(yintercept =grandS,linetype=3 ))
    print(g2 <- ggplot(df,aes(y=mC,x=time,colour=time))+ theme_bw() + geom_point() + geom_errorbar(aes(ymin=mC-sdC,ymax=mC+sdC))+ scale_color_distiller(palette = "RdYlGn",guide=FALSE)+ geom_hline(yintercept =grandC,linetype=3 ))
  } else {
    require(cowplot)
    g1 <- ggplot(df,aes(y=mS,x=time,colour=time))+ theme_bw() + geom_point() + geom_errorbar(aes(ymin=mS-sdS,ymax=mS+sdS)) + scale_color_distiller(palette = "RdYlGn",guide=FALSE)+ geom_hline(yintercept =grandS,linetype=3 ) + geom_hline(yintercept =emp$Size,linetype=2 )
    g2 <- ggplot(df,aes(y=mC,x=time,colour=time))+ theme_bw() + geom_point() + geom_errorbar(aes(ymin=mC-sdC,ymax=mC+sdC))+ scale_color_distiller(palette = "RdYlGn",guide=FALSE)+ geom_hline(yintercept =grandC,linetype=3 ) + geom_hline(yintercept =emp$Connectance,linetype=2 )
    g3 <- plot_grid(g1,g2,labels = c("A","B"),align = "h")
    save_plot(fname,g3,base_width=8,base_height=5,dpi=600)
  }
    
  return(list(g1=g1,g2=g2))
}



#' Estimation of z-scores using Meta-Web assembly model as a null 
#'
#' @param webs table with the parameters of the reference network
#' @param web_name name of the reference network
#' @param Adj Adyacency matrix for the meta-web
#' @param mig Migration parameter of the meta-Web assembly model
#' @param ext Exctinction parameter of the meta-Web assembly model
#' @param sec Secondary exctinctions parameter of the meta-Web assembly model
#' @param nsim number of simulations
#' @param final_time number of steps of the simulations 
#'
#' @return
#' @export
#'
#' @examples
calc_modularity_metaWebAssembly<- function(webs, web_name, Adj, mig,ext,sec,nsim=1000,final_time=1000){
  
  t <- webs %>% filter(Network==web_name)

  mig <- rep(mig,nrow(Adj))
  ext <- rep(ext,nrow(Adj))
  sec <- rep(sec,nrow(Adj))
  ind <- data.frame()
  require(doFuture)
  registerDoFuture()
  plan(multisession)
  
  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages=c('meweasmo','igraph'), 
                 .export = c('Adj','ext','mig','final_time')) %dopar% 
  {
    AA <- metaWebNetAssemblyCT(Adj,mig,ext,sec,final_time)
    g <- graph_from_adjacency_matrix( AA$A, mode  = "directed")
    # Select only a connected subgraph graph 
    dg <- components(g)
    g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
    mmm <- infomap(g)
    modl <- mmm$modularity
    ngrp <- length(mmm$csize)
    clus.coef <- transitivity(g, type="Global")
    cha.path  <- average.path.length(g)
    size <- vcount(g)
    links <- ecount(g)
    
    #mmm <- calc_incoherence(g)
    #qss <- calc_QSS(g,10000,ncores=4)

    bind_cols(data.frame(Size=size,Links=links,modularity=modl,ngroups=ngrp,clus.coef=clus.coef,cha.path=cha.path,Q=mmm$Q,mTI=mmm$mTI))#,qss)
  }
  plan(sequential)
  
  ind <- ind %>% mutate(gamma=t$Clustering/clus.coef,lambda=t$PathLength/cha.path,SWness=gamma/lambda)
  # 99% confidence interval
  #
  qSW <- quantile(ind$SWness,c(0.005,0.995), na.rm = TRUE)
  qmo <- quantile(ind$modularity,c(0.005,0.995), na.rm = TRUE)
  qgr <- quantile(ind$ngroups,c(0.005,0.995), na.rm = TRUE)
  mcc <- mean(ind$clus.coef, na.rm = TRUE)
  mcp <- mean(ind$cha.path, na.rm = TRUE)
  mmo <- mean(ind$modularity, na.rm = TRUE)
  mgr <- mean(ind$ngroups, na.rm = TRUE)
  mSW <- mean(t$Clustering/mcc*mcp/t$PathLength, na.rm = TRUE)
  mCI <- 1+(qSW[2]-qSW[1])/2  

  qQ <- quantile(ind$Q,c(0.005,0.995), na.rm = TRUE)
  qTI <- quantile(ind$mTI,c(0.005,0.995), na.rm = TRUE)
  mdlQ <- mean(ind$Q, na.rm = TRUE)
  mdlTI <- mean(ind$mTI, na.rm = TRUE)
  # m <-calc_incoherence(red,ti)
  
  zQ <-  (t$Q- mdlQ)/sd(ind$Q)
  zTI <- (t$mTI - mdlTI)/sd(ind$mTI) # the same as sd(ind$mTI)

  zMO <- (t$Modularity- mmo)/sd(ind$modularity)
  # 99% confidence interval
  #
  # q_qss <- quantile(ind$QSS,c(0.005,0.995),na.rm = TRUE)
  # m_qss <- mean(ind$QSS)
  # q_meing <- quantile(ind$MEing,c(0.005,0.995),na.rm = TRUE)
  # m_meing <- mean(ind$MEing)
  
  # zQSS <- (t$QSS - m_qss)/sd(ind$QSS) 
  # zMEing <- (t$MEing - m_meing)/sd(ind$MEing)

  return(list(su=tibble(mdlCC=mcc,mdlCP=mcp,mdlMO=mmo,mdlGR=mgr,SWness=mSW,SWnessCI=mCI,MOlow=qmo[1],MOhigh=qmo[2],
                    GRlow=qgr[1],GRhigh=qgr[2], mdlQ=mdlQ,mdlTI=mdlTI,Qlow=qQ[1],Qhigh=qQ[2],
                    TIlow=qTI[1],TIhigh=qTI[2],zQ=zQ,zTI=zTI,zMO=zMO
                    #mdlQSS=m_qss,QSSlow=q_qss[1],QSShigh=q_qss[2],
                    #zQSS=zQSS,mdlMEing=m_meing,MEingLow=q_meing[1],MEingHigh=q_meing[2],zMEing=zMEing
                    )
              ,sim=ind))         
}



#' Estimation of QSS z-scores using Meta-Web assembly model as a null 
#'
#' @param red This is the reference network as an igraph object
#' @param Adj Adyacency matrix for the meta-web
#' @param mig Migration parameter of the meta-Web assembly model
#' @param ext Exctinction parameter of the meta-Web assembly model
#' @param nsim number of simulations
#'
#' @return
#' @export
#'
#' @examples
calc_qss_metaWebAssembly<- function(webs, web_name, Adj, mig,ext,sec,nsim=1000,final_time=1000,ncores=0){
  
  t <- webs %>% filter(Network==web_name)
  mig <- rep(mig,nrow(Adj))
  ext <- rep(ext,nrow(Adj))
  sec <- rep(sec,nrow(Adj))
  
  ind <- data.frame()

  require(doFuture)
  registerDoFuture()
  plan(multisession)
  
  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE) %dopar% 
  {
    AA <- metaWebNetAssemblyCT(Adj,mig,ext,sec,final_time)
    g <- graph_from_adjacency_matrix( AA$A, mode  = "directed")
    dg <- components(g)
    g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
    
    size <- vcount(g)
    links <- ecount(g)
    
    # Select only a connected subgraph graph 
    # print(paste("Sim:",i, "Size:", size))
    bind_cols(data.frame(Size=size,Links=links),calc_QSS(g,1000))
  }
  plan(sequential)
  
  # 99% confidence interval
  #
  q_qss <- quantile(ind$QSS,c(0.005,0.995),na.rm = TRUE)
  m_qss <- mean(ind$QSS)
  q_meing <- quantile(ind$MEing,c(0.005,0.995),na.rm = TRUE)
  m_meing <- mean(ind$MEing)
  
  zQSS <- (t$QSS - m_qss)/sd(ind$QSS) 
  zMEing <- (t$MEing - m_meing)/sd(ind$MEing)
  return(list(su=tibble(QSS=t$QSS,mdlQSS=m_qss,QSSlow=q_qss[1],QSShigh=q_qss[2],
                zQSS=zQSS,MEing=t$MEing,mdlMEing=m_meing,MEingLow=q_meing[1],MEingHigh=q_meing[2],zMEing=zMEing),sim=ind))         
}


#' Calculate motif counts for observed network and CI for meta-web assembly model networks and Z-scores 
#'
#' @param red igraph network object
#' @param Adj Adyacency matrix for the meta-web
#' @param mig Migration parameter of the meta-Web assembly model
#' @param ext Exctinction parameter of the meta-Web assembly model
#' @param nsim number of simulation to calculate random networks with the same nodes and links
#'
#' @return data.frame with all the results
#' @export
#'
#' @examples
calc_motif_metaWebAssembly<- function(motdf,fw_name, Adj, mig,ext,sec, nsim=1000,final_time=1000)
{
  mig <- rep(mig,nrow(Adj))
  ext <- rep(ext,nrow(Adj))
  sec <- rep(sec,nrow(Adj))
  
  require(doFuture)
  registerDoFuture()
  plan(multisession)
  
  ind <- data.frame()

  ind <- foreach(i=1:nsim,.combine='rbind',.inorder=FALSE,.packages=c('meweasmo','igraph'), 
                 .export = c('Adj','ext','mig','final_time')) %dopar% 
  {
    AA <- metaWebNetAssemblyCT(Adj,mig,ext,sec,final_time)
    g <- graph_from_adjacency_matrix( AA$A, mode  = "directed")
    # Select only a connected subgraph graph 
    dg <- components(g)
    g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
    
    mot <- triad_census(g)
    mot[4] # Exploitative competition
    mot[5] # Apparent competition
    mot[6] # Tri-trophic chain
    mot[9] # Omnivory
    
    data.frame(explComp=mot[4],apprComp=mot[5],triTroph=mot[6],omnivory=mot[9])
  }
  plan(sequential)
  
  # 99% confidence interval
  #
  
  qEC <- quantile(ind$explComp,c(0.005,0.995))
  qAC <- quantile(ind$apprComp,c(0.005,0.995))
  qTT <- quantile(ind$triTroph,c(0.005,0.995))
  qOM <- quantile(ind$omnivory,c(0.005,0.995))
  
  # Calculate motif for the original network
  obs <- motdf %>% filter(Network==fw_name)
  
  zEC <- (obs$explComp - mean(ind$explComp))/sd(ind$explComp)
  zAC <- (obs$apprComp - mean(ind$apprComp))/sd(ind$apprComp)
  zTT <- (obs$triTroph - mean(ind$triTroph))/sd(ind$triTroph)
  zOM <- (obs$omnivory - mean(ind$omnivory))/sd(ind$omnivory)
  
  return(tibble(explComp=obs$explComp,apprComp=obs$apprComp,triTroph=obs$triTroph,omnivory=obs$omnivory,zEC=zEC,zAC=zAC,zTT=zTT,zOM=zOM,EClow=qEC[1],EChigh=qEC[2],AClow=qAC[1],AChigh=qAC[2],TTlow=qTT[1],TThigh=qTT[2],OMlow=qOM[1],OMhigh=qOM[2]))         
}    



#' Title Plot 5 simulations of net assembly model time series of S and L only the last timeW steps are ploted
#'
#' @param metaW meta-web adjacency matrix 
#' @param m     migration
#' @param q     probability of link
#' @param a     extinction
#' @param timeW time window used
#'
#' @return
#' @export
#'
#' @examples
plot_NetAssemblyModel_sims <- function(metaW,m, q, a, tf,timeW){
  require(viridis)

  if(tf<timeW) stop("timeW parameter must be less than the time of the simulation")
  
  dfA <- data.frame()
  
  for(n in 1:5){
    AA <- metaWebNetAssembly(metaW,m,q,a,tf)
    tdfA <- data.frame(S=AA$S[(tf-timeW):tf],L=as.numeric(AA$L[(tf-timeW):tf]),T=c((tf-timeW):tf))
    tdfA$C <- tdfA$L/(tdfA$S*tdfA$S)
    tdfA$sim <- n
    dfA <- bind_rows(dfA,tdfA)
  }
  gS <- ggplot(dfA, aes(x=T,y=S,colour=sim)) + geom_point() + theme_bw() + geom_hline(yintercept = mean(dfA$S),linetype = 2,colour="grey50") + scale_color_viridis(guide=FALSE)
  print(gS)
  gL <- ggplot(dfA, aes(x=T,y=L,colour=sim)) + geom_point() + theme_bw() + ylab("L") + geom_hline(yintercept = mean(dfA$L),linetype = 2,colour="grey50") + scale_color_viridis(guide=FALSE)
  print(gL)
  gC <- ggplot(dfA, aes(x=T,y=C,colour=sim)) + geom_point() + theme_bw() + ylab("C") + geom_hline(yintercept = mean(dfA$C),linetype = 2,colour="grey50") + scale_color_viridis(guide=FALSE)
  print(gC)
  return(list(gS=gS,gL=gL,gC=gC))
}




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
  #p <- progressor(steps = no_sims )
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
    #p(sprintf("Sim#%d (PID %d)", i, Sys.getpid()))
    dfA <- data.frame(S=AA$S[(tf-100):tf], L=as.numeric(AA$L[(tf-100):tf]), T=c((tf-100):tf))
    data.frame(m=mm[1], a=aa[1], se=ee[1], S=mean(dfA$S), L=mean(dfA$L), C=mean(dfA$L)/(mean(dfA$S)*mean(dfA$S)))
  }, future.seed = TRUE)
  plan(sequential)
  sim <- do.call(rbind, sim_list)
  return(sim)
}


#' Return the names of the basal species: nodes with no incoming links 
#'
#' @param g igraph network
#'
#' @return the basal species 
#' @export
#'
#' @examples
basal_species <- function(g){
  deg <- degree(g, mode="in") # calculate the in-degree: the number of preys
  
  V(g)$indegree <-  deg
  
  basal <- V(g)[indegree==0]
  
  V(g)[basal]$name
  
}


calc_compare_motif <- function(redl, network_name,nsims=1000,rnd_seed=123){
  #
  # Motif
  #
  mot <- triad_census(redl)
  
  # mot[4] = Exploitative competition
  # mot[5] = Apparent competition
  # mot[6] = Tri-trophic chain
  # mot[9] = Omnivory
  # mot[10] = Loop 
  
  freq_all_Motif <- data.frame(Network=network_name,t(mot))
  
  set.seed(rnd_seed)
  motif_ER <- calc_motif_random(redl,nsims) %>% mutate(Network=network_name)
  
  return(list(moter=motif_ER,fallm=freq_all_Motif))
}




#' Fit metaweb assembly model from simulations 
#'
#' @param webs data frame with info about the web we are fitting
#' @param web_name name of the web we are fitting
#' @param sims data frame with simulations
#' @param tol  tolerance of the fitting for the plot
#' @param plot logical if true make and save a plot of the fitting
#'
#' @return
#' @export
#'
#' @examples
fit_metaWebAssembly_model <- function(emp_S, emp_C, web_name, sims, tol=0.1, plot=TRUE){
  # Filter simulations for the given web_name
  sims <- sims %>% filter(Metaweb == web_name)
  orig_tol <- tol
  while (TRUE) {
    sel <- sims %>%
      group_by(model) %>%
      mutate(alpha = m/a) %>%
      filter(S > emp_S * (1 - tol), S < emp_S * (1 + tol)) %>%
      arrange(S)
    if (nrow(sel) > 0) break
    tol <- tol * 2
    warning(sprintf("No simulations within tolerance %.3f for %s, doubling tolerance to %.3f", tol/2, web_name, tol))
  }
  # Fit using the distance to S and C
  sel <- sel %>%
    group_by(model, model_type) %>%
    mutate(
      distance = sqrt(((emp_S - S)/emp_S)^2 + ((emp_C - C)/emp_C)^2),
      min_dist = (distance == min(distance))
    ) %>%
    arrange(distance)
  fittedmw <- sel %>% filter(min_dist) %>% mutate(Network = web_name, fit_type = "S-C")
  if (plot) {
    print(
      ggplot(sel, aes(S, C, color = min_dist)) +
        geom_point(alpha = 0.5) + theme_bw() +
        annotate("point", x = emp_S, y = emp_C, color = "red", size = 3) +
        facet_wrap(~model) +
        scale_color_viridis_d(name = "", labels = c("Simulations", "Empirical", "Fit"))
    )
    ggsave(paste0("Figures/Metaweb_fit_", web_name, "_byModel.png"), width = 8, height = 5, units = "in", dpi = 600)
  }
  return(fittedmw)
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

list_obj_sizes <- function(list_obj=ls(envir=.GlobalEnv)){ 	
  sizes <- sapply(list_obj, function(n) object.size(get(n)), simplify = FALSE) 	
  print(sapply(sizes[order(-as.integer(sizes))], function(s) format(s, unit = 'auto'))) 
} 


#' Simulate and Plot meta web assembly model
#'
#' @param meta igraph metaweb 
#' @param fitted Data frame with the fitted parameters if more than one record select the first
#' @param webst Data frame with the information about networks
#' @param netname Name of the network
#' 
#' @return
#' @export
#'
#' @examples
simulate_plot_metaweb_assembly <- function(meta,fitted,webst,netname){
  # Model with probability of secundary extinction se>0 & se<1
  #
  # Potter
  #
  A <- as_adjacency_matrix(meta,sparse=F)

  tf <- 1000
  set.seed(1110) 
  f <- fitted %>% filter(Network==netname) 
  mm <- rep(f$m[1],nrow(A))
  aa <- rep(f$a[1],nrow(A))
  se <- rep(f$se[1],nrow(A))
  AA <- metaWebNetAssemblyCT(A,mm,aa,se,tf)
  
  figname <- paste0("Figures/",netname,"_metawebSim_avg.png")
  # Running averages
  plot_NetAssemblyModel_eqw(AA,50,figname,webst %>% filter(Network==netname))
  
  # Time series plot
  figname <- paste0("Figures/",netname,"_metawebSim_ts.png")
  plot_NetAssemblyModel(AA,300,figname,webst %>% filter(Network==netname))
  
  return(f[1,])
}


calc_topoRoles_metaWebAssembly <- function(web_name,meanS, Adj, mig,ext,sec,final_time=1000){
  
  mig <- rep(mig,nrow(Adj))
  ext <- rep(ext,nrow(Adj))
  sec <- rep(sec,nrow(Adj))
  
  while(TRUE){
    #message(paste("\nmetaWebNetAssemblyCT", dim(Adj),"trace", sum(diag(Adj)), mig,ext,sec))
    AA <- metaWebNetAssemblyCT(Adj,mig,ext,sec,final_time)
    g <- graph_from_adjacency_matrix( AA$A, mode  = "directed")
    dg <- components(g)
    g <- induced_subgraph(g, which(dg$membership == which.max(dg$csize)))
    #
    # Added this Loop and condition to prevent very small networks when exist a posibility 
    #
    if( vcount(g) > meanS*.8) break  
  }
  modulos<-cluster_spinglass(g)
  
  topo_roles <- incremental_topoRoles(g, web_name)
  hub_conn <- plot_topological_roles(topo_roles,g,modulos) %>% mutate(Network=web_name)

  return(list(tr=topo_roles,hc=hub_conn))
}


#' Simulate Food Webs with Metaweb Assembly Model and Compute Network Metrics
#'
#' This function performs multiple simulations of the metaweb assembly model (`metaWebNetAssemblyCT`)
#' using fitted parameters (`m`, `a`, `se`) sampled from a data frame. It calculates SVD entropy,
#' effective rank, and Infomap modularity for each simulation.
#'
#' @param ig igraph object representing the metaweb.
#' @param fitted_params A data.frame with columns `m`, `a`, `se`, and `Metaweb` for parameter sampling.
#' @param metaweb_name Character, the name of the metaweb to filter parameters (must match `Metaweb` column).
#' @param nsim Number of simulations to perform (default: 1000).
#' @param final_time Total time steps for the assembly simulation (default: 1000).
#'
#' @return A `tibble` with columns: `Sim`, `Metaweb`, `S`, `L`, `Entropy`, `Rank`, `Modularity`.
#' @importFrom igraph graph_from_adjacency_matrix components induced_subgraph vcount ecount
#' @importFrom multiweb calc_svd_entropy run_infomap
#' @importFrom dplyr bind_rows
#' @export
simulate_metaweb_metrics <- function(Adj, fitted_params, metaweb_name,
                                     nsim = 1000, final_time = 1000) {
  
  require(igraph)
  require(dplyr)
  require(future.apply)
  
  param_set <- fitted_params %>% filter(Metaweb == metaweb_name)
  if (nrow(param_set) == 0) stop("No fitted parameters found for this metaweb.")
  
  num_cores <- future::availableCores()
  plan(multisession, workers = num_cores)
  # if (supportsMulticore()) {
  #   plan(multicore, workers = num_cores)
  # } else {
  #   message("multicore not supported on this OS, falling back to multisession.")
  #   plan(multisession, workers = num_cores)
  # }

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
      comm <- multiweb::run_infomap(g)
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

#' Plot Observed vs Simulated Network Metrics
#'
#' @param network_info A data frame with observed metrics including `site`, `S`, `C`, `Entropy`, `Modularity`, `Rank`, etc.
#' @param simulated_metrics A data frame with simulated metrics including `Metaweb`, `S`, `L`, `Entropy`, `Rank`, `Modularity`, etc.
#' @param metrics Character vector of metrics to plot. If NULL, it will use the intersection of numeric columns.
#'
#' @return A named list of ggplot objects comparing each metric across empirical networks.
#' @export
#'
#' @import ggplot2 dplyr
plot_empirical_vs_simulated_metrics <- function(network_info, simulated_metrics, metrics = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # Prepare simulated data
  sim_data <- simulated_metrics %>%
    rename(site = Metaweb) %>%
    filter(site %in% network_info$site)
  
  # Prepare observed data
  obs_data <- network_info %>%
    dplyr::select(site, everything()) %>%
    filter(site %in% sim_data$site)
  
  # Determine metric columns if not provided
  if (is.null(metrics)) {
    numeric_cols <- purrr::map_lgl(obs_data, is.numeric)
    candidate_metrics <- names(obs_data)[numeric_cols]
    metrics <- intersect(candidate_metrics, names(sim_data))
    metrics <- setdiff(metrics, c("latitude", "depth_m", "area_km2"))  # drop geographic or non-network
  }
  
  # Helper for plotting
  plot_metric <- function(metric_name) {
    ggplot(sim_data, aes(x = name, y = .data[[metric_name]])) +
      geom_boxplot(fill = "gray90", color = "gray50") +
      geom_point(data = obs_data, aes(x = site, y = .data[[metric_name]]),
                 color = "red", size = 2.5) +
      labs(y = metric_name, x = NULL) +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 14, face = "bold"))
  }
  
  # Generate plots
  plots <- setNames(lapply(metrics, plot_metric), metrics)
  return(plots)
}



#' Plot Simulated Metric Densities with Observed Value Lines
#'
#' @param network_info A data frame with observed values: must include `site`, `S`, `C`, `entropy`, `modularity`.
#' @param simulated_metrics A tibble with columns: `Metaweb`, `S`, `L`, `Entropy`, `Rank`, `Modularity`.
#' @param metrics Character vector of metric names to plot (default all).
#' @return A list of ggplot objects, one per metric.
#' @export
#'
#' @import ggplot2 dplyr viridisLite
plot_metric_densities_with_empirical <- function(network_info, simulated_metrics,
                                                 metrics = c("S", "C", "Entropy", "Rank", "Modularity")) {
  require(dplyr)
  require(ggplot2)
  require(viridisLite)
  
  # Ensure connectance is available in simulated data
  sim <- simulated_metrics %>%
    mutate(C = L / (S^2)) %>%
    rename(site = Metaweb)
  
  obs <- network_info %>%
    mutate(site = as.character(site)) %>%
    dplyr::select(site, S, C, Entropy, Rank, Modularity)  # Fallback in case Rank is not included
  
  plot_list <- list()
  
  site_levels <- unique(obs$site)
  site_colors <- setNames(viridis(length(site_levels)), site_levels)
  
  for (metric in metrics) {
    if (!metric %in% colnames(sim)) next
    
    sim_metric <- sim %>% filter(!is.na(.data[[metric]]))
    obs_metric <- obs %>% filter(!is.na(.data[[metric]]))
    
    p <- ggplot(sim_metric, aes(x = .data[[metric]], color = site, fill = site)) +
      geom_density(alpha = 0.3) +
      geom_vline(data = obs_metric, aes(xintercept = .data[[metric]], color = site),
                 size = 1, linetype = "solid") +
      scale_color_manual(values = site_colors) +
      scale_fill_manual(values = site_colors) +
      labs(title = paste("Distribution of", metric, "with empirical values"),
           x = metric, y = "Density") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
    
    plot_list[[metric]] <- p
  }
  
  return(plot_list)
}


#' Plot Empirical Values Against Simulated Metric Distributions with CI
#'
#' This function plots empirical network metrics against 95% confidence intervals
#' of simulated values. Sites are ordered by latitude from low to high.
#'
#' @param network_info Data frame with columns: `site`, `latitude`, and metrics (`S`, `C`, `Entropy`, `Rank`, `Modularity`).
#' @param simulated_metrics Data frame with simulated values, must include `Metaweb`, `S`, `L`, `Entropy`, `Rank`, `Modularity`.
#' @param metrics Character vector of metric names to visualize.
#'
#' @return A named list of ggplot2 objects (one per metric).
#' @export
#'
#' @import ggplot2 dplyr tidyr viridisLite
plot_empirical_vs_simulated_ci <- function(network_info,
                                           simulated_metrics,
                                           metrics = c("S", "C", "Entropy", "Rank", "Modularity")) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(viridisLite)
  
  # Compute simulated connectivity
  sim <- simulated_metrics %>%
    mutate(C = L / (S^2)) %>%
    rename(site = Metaweb)
  
  # Prepare empirical data
  obs <- network_info %>%
    mutate(C = C, site = as.character(site)) %>%
    dplyr::select(site, latitude, all_of(metrics))
  
  site_lat <- obs %>% dplyr::select(site, latitude) %>% distinct()
  site_order <- site_lat %>% arrange(latitude) %>% pull(site)
  
  # Color palette
  site_colors <- setNames(viridis(length(site_order)), site_order)
  
  plot_list <- list()
  
  for (metric in metrics) {
    if (!metric %in% names(sim)) next
    
    # Confidence intervals from simulations
    ci_data <- sim %>%
      filter(!is.na(.data[[metric]])) %>%
      group_by(site) %>%
      summarise(
        q2.5 = quantile(.data[[metric]], 0.025, na.rm = TRUE),
        q97.5 = quantile(.data[[metric]], 0.975, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Join with empirical
    obs_metric <- obs %>%
      dplyr::select(site, value = !!sym(metric))
    
    df_plot <- left_join(ci_data, obs_metric, by = "site") %>%
      left_join(site_lat, by = "site") %>%
      mutate(site = factor(site, levels = site_order))
    
    p <- ggplot(df_plot, aes(y = site)) +
      geom_errorbarh(aes(xmin = q2.5, xmax = q97.5), height = 0.25, color = "gray60") +
      geom_point(aes(x = value, color = site), size = 3) +
      scale_color_manual(values = site_colors) +
      labs(
#        title = paste("Empirical vs Simulated CI:", metric),
        x = metric, y = "Site (Ordered by Latitude)"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
    
    plot_list[[metric]] <- p
  }
  
  return(plot_list)
}

#' Plot Network Metrics with Simulated Confidence Intervals + Quantile Regression
#'
#' Fits quantile regressions (tau = 0.25, 0.5, 0.75) and plots regression lines.
#' Returns a dataframe with term-level estimates, SE, and p-values.
#'
#' @param network_info Data frame with empirical metrics (`site`, `latitude`, `area_km2`, etc).
#' @param simulated_metrics Simulated metrics from metaweb simulations (`Metaweb`, `S`, `L`, etc).
#' @param metrics Vector of metric names to plot (default: c("C", "SVDComplexity", "Modularity")).
#' @param xvar Variable for x-axis, either `"latitude"`, `"log_area"`, or `"impact_mean"`.
#' @param fit_model If TRUE, fits quantile regression and returns tidy dataframe (default: TRUE).
#'
#' @return A list with:
#'   - `plots`: Named list of ggplot objects.
#'   - `models`: A tidy dataframe with term estimates, SEs, p-values for each quantile & metric.
#' @export
plot_metric_vs_latitude_ci <- function(network_info,
                                       simulated_metrics,
                                       metrics = c("C", "SVDComplexity", "Modularity"),
                                       xvar = "latitude",
                                       fit_model = TRUE) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggrepel)
  library(viridisLite)
  library(quantreg)
  library(broom)   # for tidy()
  library(purrr)
  
  if ("Metaweb" %in% names(simulated_metrics)) {
    sim <- simulated_metrics %>% rename(site = Metaweb)
  } else {
    sim <- simulated_metrics
  }
  
  obs <- network_info %>%
    mutate(site = as.character(site),
           log_area = log(area_km2)) %>%
    dplyr::select(site, name, C, S, latitude, log_area, depth_m, impact_mean)
  
  x_label <- case_when(
    xvar == "log_area" ~ "Log Area (km²)",
    xvar == "latitude" ~ "Latitude",
    xvar == "impact_mean" ~ "Human impact",
    TRUE ~ xvar
  )
  
  site_order <- obs %>% arrange(latitude) %>% pull(site)
  site_colors <- setNames(viridis(length(site_order)), site_order)
  
  plot_list <- list()
  model_df <- tibble()
  
  taus <- c(0.25, 0.5, 0.75)
  
  for (metric in metrics) {
    if (!metric %in% names(sim)) next
    
    # Confidence intervals from simulations
    ci_data <- sim %>%
      filter(!is.na(.data[[metric]])) %>%
      group_by(site) %>%
      summarise(
        q2.5 = quantile(.data[[metric]], 0.025, na.rm = TRUE),
        q97.5 = quantile(.data[[metric]], 0.975, na.rm = TRUE),
        mean = mean(.data[[metric]], na.rm = TRUE),
        .groups = "drop"
      )
    
    df_plot <- obs %>%
      left_join(ci_data, by = "site") %>%
      mutate(color = site_colors[site],
             xval = case_when(
               xvar == "log_area" ~ log_area,
               xvar == "latitude" ~ latitude,
               xvar == "impact_mean" ~ impact_mean,
               TRUE ~ NA_real_
             ))
    
    # Model data
    df_test <- sim %>%
      dplyr::select(site, value = !!sym(metric)) %>%
      left_join(obs, by = "site")
    
    # Fit models at 3 quantiles
    mods <- purrr::map(taus, function(tau) {
      if(metric == "C") {
        rq(value ~ S + latitude + log_area + impact_mean,
           data = df_test, tau = tau)
      } else {
        rq(value ~ S + C + latitude + log_area + impact_mean,
           data = df_test, tau = tau)
      }
    })
    
    if (fit_model) {
      tidy_mods <- map2_dfr(mods, taus, function(m, tau) {
        broom::tidy(m, se = "boot") %>%
          mutate(tau = tau, Metric = metric)
      })
      model_df <- bind_rows(model_df, tidy_mods)
    }
    
    # Prediction grid
    x_seq <- seq(min(df_plot$xval, na.rm = TRUE), max(df_plot$xval, na.rm = TRUE), length.out = 100)
    newdata <- obs %>%
      summarise(latitude = mean(latitude, na.rm = TRUE),
                log_area = mean(log_area, na.rm = TRUE),
                S = mean(S, na.rm = TRUE),
                C = mean(C, na.rm = TRUE),
                impact_mean = mean(impact_mean, na.rm = TRUE)) %>%
      slice(rep(1, 100)) %>%
      mutate(!!xvar := x_seq)
    
    df_lines <- map2_dfr(mods, taus, function(m, tau) {
      tibble(xval = x_seq,
             pred = predict(m, newdata = newdata),
             tau = tau)
    })
    
    p <- ggplot(df_plot, aes(x = xval, y = mean)) +
      geom_linerange(aes(ymin = q2.5, ymax = q97.5, color = site), linewidth = 1) +
      geom_point(aes(color = site), size = 3) +
      geom_line(data = df_lines, aes(x = xval, y = pred, linetype = factor(tau)),
                linewidth = 0.8, inherit.aes = FALSE) +
      scale_linetype_manual(values = c("dotted", "dashed", "dotdash"),
                            name = "Quantile",
                            labels = c("25%", "50%", "75%")) +
      scale_color_manual(values = site_colors,
                         limits = df_plot$site,
                         labels = df_plot$name) +
      labs(x = x_label, y = metric) +
      theme_bw(base_size = 14) +
      theme(legend.position = "none")
    
    plot_list[[metric]] <- p
  }
  
  return(if (fit_model) list(plots = plot_list, models = model_df) else list(plots = plot_list))
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
#' @return A data frame with best QSS and MEing_stable per site.
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
    
    # Check if the qss has reach the trhreshold
    
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
  
  # Collapse and take best selfDamping per site
  all_results <- bind_rows(results_raw)
  
  # best_results <- all_results %>%
  #   group_by(site) %>%
  #   slice_min(selfDamping, with_ties = FALSE) %>%
  #   ungroup()
  
  return(all_results)
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
#' @param qss_threshold Value above which we stop computing more values (default: 0.1).
#'
#' @return A data frame with QSS, MEing_stable and selfDamping per site.
#' @export
calculate_network_stability_threshold <- function(networks,
                                                 sites,
                                                 selfDamping_range = seq(-1, -10, by = -0.5),
                                                 nsim = 1000,
                                                 negative = -10,
                                                 positive = 1,
                                                 qss_threshold = 0.1) {
  stopifnot(length(networks) == length(sites))
  
  library(dplyr)
  library(future.apply)
  library(purrr)

  # Shared environment to track which sites passed QSS threshold
  qss_reached_sites <- new.env(parent = emptyenv())
  qss_reached_sites$status <- setNames(rep(FALSE, length(sites)), sites)
  
  # Expand combinations: each row = 1 task
  combos <- expand.grid(
    idx = seq_along(networks),
    selfDamping = selfDamping_range,
    KEEP.OUT.ATTRS = FALSE
  )
  
  combos_list <- split(combos, seq(nrow(combos)))
  
  results_raw <- future_lapply(combos_list, function(row) {
    i <- row$idx
    sd <- row$selfDamping
    g <- networks[[i]]
    site <- sites[[i]]
    
    # Skip if site already met QSS threshold
    if (qss_reached_sites$status[[site]]) return(NULL)
    
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
    
    if (qss > qss_threshold) {
      assign("status",
             replace(qss_reached_sites$status, site, TRUE),
             envir = qss_reached_sites)
    }
    
    if (qss > 0) {
      MEing_stable <- mean(eigen_vals[eigen_vals < 0], na.rm = TRUE)
      return(tibble(site = site, QSS = qss, MEing_stable = MEing_stable, selfDamping = sd))
    } else {
      return(NULL)
    }
  })
  
  # Collect results
  all_results <- bind_rows(results_raw)
  
  return(all_results)
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






plot_metric_vs_latitude_ci_gam <- function(network_info,
                                           simulated_metrics,
                                           metrics = c("S", "C", "Rank", "Entropy", "Modularity", "MEing_stable"),
                                           xvar = "latitude",
                                           fit_model = TRUE) {
  library(dplyr); library(ggplot2); library(tidyr); library(ggrepel)
  library(viridisLite); library(mgcv); library(tibble)
  
  # rename Metaweb -> site if present
  if ("Metaweb" %in% names(simulated_metrics)) {
    sim <- simulated_metrics %>% rename(site = Metaweb)
  } else {
    sim <- simulated_metrics
  }
  
  obs <- network_info %>%
    mutate(site = as.character(site),
           log_area = log(area_km2)) %>%
    dplyr::select(site, name, C, S, latitude, log_area, depth_m, impact_mean)
  
  x_label <- dplyr::case_when(
    xvar == "log_area" ~ "Log Area (km²)",
    xvar == "latitude" ~ "Latitude",
    xvar == "impact_mean" ~ "Human impact",
    xvar == "depth_m" ~ "Depth (m)",
    TRUE ~ xvar
  )
  
  site_order <- obs %>% arrange(latitude) %>% pull(site)
  site_colors <- setNames(viridis(length(site_order)), site_order)
  
  plot_list <- list()
  model_list <- list()
  deviance_list <- list()
  
  for (metric in metrics) {
    if (!metric %in% names(sim)) next
    
    # Work on a local copy so we don't mutate sim across metrics
    sim2 <- sim
    
    # If this metric is MEing_stable we transform to positive for Gamma fitting:
    is_transformed <- identical(metric, "MEing_stable")
    if (is_transformed) {
      # create a transformed column with positive values for modelling
      # do not overwrite the original column in sim; create metric_tmp
      sim2 <- sim2 %>% mutate(metric_tmp = - .data[[metric]])
      model_value_name <- "metric_tmp"
    } else {
      sim2 <- sim2 %>% mutate(metric_tmp = .data[[metric]])
      model_value_name <- "metric_tmp"
    }
    
    # compute empirical "intervals" (from the simulations) for plotting
    ci_data <- sim2 %>%
      filter(!is.na(.data[[model_value_name]])) %>%
      group_by(site) %>%
      summarise(
        q2.5 = quantile(.data[[model_value_name]], 0.025, na.rm = TRUE),
        q97.5 = quantile(.data[[model_value_name]], 0.975, na.rm = TRUE),
        mean = mean(.data[[model_value_name]], na.rm = TRUE),
        .groups = "drop"
      )
    # flip back for display if transformed
    if (is_transformed) {
      ci_data <- ci_data %>% mutate(across(c(q2.5, q97.5, mean), ~ - .x))
    }
    
    df_plot <- obs %>%
      left_join(ci_data, by = "site") %>%
      mutate(color = site_colors[site],
             xval = dplyr::case_when(
               xvar == "log_area" ~ log_area,
               xvar == "latitude" ~ latitude,
               xvar == "impact_mean" ~ impact_mean,
               xvar == "depth_m" ~ depth_m,
               xvar == "S" ~ S,
               TRUE ~ NA_real_
             ))
    
    # build model dataset (using the transformed positive column when needed)
    df_test <- sim2 %>%
      dplyr::select(site, value = !!sym(model_value_name)) %>%
      left_join(obs, by = "site") %>%
      drop_na(value, S, latitude, log_area, impact_mean) %>%
      mutate(site = as.factor(site))   # ensure site is a factor for bs="re"
    
    # optionally subsample to limit compute (as you did before)
    if (nrow(df_test) > 1000) {
      df_test <- df_test %>% group_by(site) %>% slice_sample(n = 500) %>% ungroup()
    }
    
    # GAM formula (same for MEing_stable but family differs)
#    if(metric=="C") { 
      formula_gam <- as.formula("value ~ s(S, k = 5) + s(log_area, k = 5) + s(latitude, k = 5) + s(impact_mean, k = 5) + s(depth_m, k=5) + s(site, bs = 're')")
#    } else { 
#      formula_gam <- as.formula("value ~ s(S, k = 5) + s(C, k = 5) + s(log_area, k = 5) + s(latitude, k = 5) + s(impact_mean, k = 5) + s(site, bs = "re")")
#    }
    
    if (is_transformed) {
      # value is positive (we made it so), fit Gamma with log link
      mod_gam <- gam(formula_gam, data = df_test, method = "REML", family = Gamma(link = "log"))
    } else {
      # use Student-t family for robustness (scat())
      mod_gam <- gam(formula_gam, data = df_test, method = "REML", family = scat())
    }
    
    if (fit_model) {
      model_list[[metric]] <- mod_gam
      
      s <- summary(mod_gam)
      # s$s.table exists for smooth terms; if not, create empty tibble
      if (!is.null(s$s.table)) {
        s_tab <- as.data.frame(s$s.table) %>% tibble::rownames_to_column("Term")
        contrib <- s_tab %>%
          mutate(
            Metric = metric,
            contrib_raw = if ("Chi.sq" %in% names(s_tab)) Chi.sq else F,
            PartialDeviance = contrib_raw / sum(contrib_raw, na.rm = TRUE) * 100,
            TotalDevianceExplained = round(100 * s$dev.expl, 1)
          ) %>%
          dplyr::select(Metric, Term, edf, `p-value`, PartialDeviance, TotalDevianceExplained)
      } else {
        contrib <- tibble(Metric = metric, Term = NA_character_, edf = NA_real_, `p-value` = NA_real_,
                          PartialDeviance = NA_real_, TotalDevianceExplained = round(100 * s$dev.expl, 1),
                          Transformed = is_transformed)
      }
      deviance_list[[metric]] <- contrib
    }
    
    # Build newdata for partial prediction (holding others at mean)
    x_seq <- seq(min(df_plot$xval, na.rm = TRUE), max(df_plot$xval, na.rm = TRUE), length.out = 100)
    newdata <- obs %>%
      summarise(
        latitude = mean(latitude, na.rm = TRUE),
        log_area = mean(log_area, na.rm = TRUE),
        S = mean(S, na.rm = TRUE),
        depth_m = mean(depth_m, na.rm = TRUE),
        impact_mean = mean(impact_mean, na.rm = TRUE)
      ) %>%
      slice(rep(1, 100)) %>%
      mutate(!!xvar := x_seq)
    ref_site <- levels(df_test$site)[1]
    
    newdata <- newdata %>%
      mutate(site = factor(ref_site, levels = levels(df_test$site)))
    
    # Predict on link scale, get se on link scale, back-transform using family$linkinv
    pred_obj <- predict(mod_gam, newdata = newdata, se.fit = TRUE, type = "link")
    eta <- pred_obj$fit
    se_eta <- pred_obj$se.fit
    linkinv <- mod_gam$family$linkinv
    
    resp <- linkinv(eta)
    lwr_resp <- linkinv(eta - 2 * se_eta)
    upr_resp <- linkinv(eta + 2 * se_eta)
    
    # if metric was transformed (we fit on -MEing_stable), flip sign back
    if (is_transformed) {
      resp <- -resp
      # careful with upper/lower ordering after sign flip
      lwr_resp <- -lwr_resp
      upr_resp <- -upr_resp
      # ensure lwr <= upr
      tmp_min <- pmin(lwr_resp, upr_resp)
      tmp_max <- pmax(lwr_resp, upr_resp)
      lwr_resp <- tmp_min
      upr_resp <- tmp_max
    }
    
    df_line <- data.frame(
      xval = x_seq,
      pred = resp,
      lwr = lwr_resp,
      upr = upr_resp
    )
    
    p <- ggplot(df_plot, aes(x = xval, y = mean)) +
      geom_linerange(aes(ymin = q2.5, ymax = q97.5, color = site), linewidth = 1) +
      geom_point(aes(color = site), size = 3) +
      geom_line(data = df_line, aes(x = xval, y = pred), linetype = "dashed") +
      geom_ribbon(data = df_line, aes(x = xval, ymin = lwr, ymax = upr), inherit.aes = FALSE,
                  fill = "grey70", alpha = 0.3) +
      scale_color_manual(values = site_colors,
                         limits = df_plot$site,
                         labels = df_plot$name) +
      labs(x = x_label, y = metric) +
      theme_bw(base_size = 14) +
      theme(legend.position = "none")
    
    plot_list[[metric]] <- p
  }
  
  return(list(plots = plot_list,
              models = model_list,
              deviance = bind_rows(deviance_list)))
}


plot_metric_vs_latitude_bayes <- function(network_info,
                                          simulated_metrics,
                                          metrics = c("S", "C", "Rank", "Entropy", "Modularity", "MEing_stable"),
                                          xvar = "latitude",
                                          fit_model = TRUE,
                                          iters = 2000, chains = 4, cores = 4) {
  library(dplyr); library(ggplot2); library(tidyr); library(ggrepel)
  library(viridisLite); library(tibble); library(brms)
  
  if ("Metaweb" %in% names(simulated_metrics)) {
    sim <- simulated_metrics %>% rename(site = Metaweb)
  } else {
    sim <- simulated_metrics
  }
  
  obs <- network_info %>%
    mutate(site = as.character(site),
           log_area = log(area_km2)) %>%
    dplyr::select(site, name, C, S, latitude, log_area, depth_m, impact_mean)
  
  x_label <- dplyr::case_when(
    xvar == "log_area" ~ "Log Area (km²)",
    xvar == "latitude" ~ "Latitude",
    xvar == "impact_mean" ~ "Human impact",
    xvar == "S" ~ "S",
    TRUE ~ xvar
  )
  
  site_order <- obs %>% arrange(latitude) %>% pull(site)
  site_colors <- setNames(viridis(length(site_order)), site_order)
  
  plot_list <- list()
  model_list <- list()
  
  for (metric in metrics) {
    if (!metric %in% names(sim)) next
    
    sim2 <- sim
    is_transformed <- identical(metric, "MEing_stable")
    if (is_transformed) {
      sim2 <- sim2 %>% mutate(metric_tmp = - .data[[metric]])
      model_value_name <- "metric_tmp"
    } else {
      sim2 <- sim2 %>% mutate(metric_tmp = .data[[metric]])
      model_value_name <- "metric_tmp"
    }
    
    ci_data <- sim2 %>%
      filter(!is.na(.data[[model_value_name]])) %>%
      group_by(site) %>%
      summarise(
        q2.5 = quantile(.data[[model_value_name]], 0.025, na.rm = TRUE),
        q97.5 = quantile(.data[[model_value_name]], 0.975, na.rm = TRUE),
        mean = mean(.data[[model_value_name]], na.rm = TRUE),
        .groups = "drop"
      )
    if (is_transformed) {
      ci_data <- ci_data %>% mutate(across(c(q2.5, q97.5, mean), ~ - .x))
    }
    
    df_plot <- obs %>%
      left_join(ci_data, by = "site") %>%
      mutate(color = site_colors[site],
             xval = case_when(
               xvar == "log_area" ~ log_area,
               xvar == "latitude" ~ latitude,
               xvar == "impact_mean" ~ impact_mean,
               TRUE ~ NA_real_
             ))
    
    df_test <- sim2 %>%
      dplyr::select(site, value = !!sym(model_value_name)) %>%
      left_join(obs, by = "site") %>%
      drop_na(value, S, latitude, log_area, impact_mean)
    
    if (nrow(df_test) > 2000) {
      df_test <- df_test %>% group_by(site) %>% slice_sample(n = 200) %>% ungroup()
    }
    
    if (fit_model) {
      # standardize predictors for stability
      df_test <- df_test %>%
        mutate(across(c(value, S, log_area, latitude, impact_mean), scale))
      
      # Student-t likelihood with random intercept by site
      formula <- bf(
        value ~ 1 + s(S) + s(log_area) + s(latitude) + s(impact_mean) + (1 | site),
        family = student()
      )
      
      priors <- c(
        prior(normal(0, 1), class = "Intercept"),
        prior(normal(0, 0.5), class = "b"),
        prior(exponential(1), class = "sigma"),
        prior(exponential(1), class = "sd"),
        prior(exponential(1), class = "nu"))
      
      m <- brm(formula,
               data = df_test,
               prior = priors,
               chains = chains, cores = cores, iter = iters,
               backend = "cmdstanr")  # if installed, faster & more stable
      
      model_list[[metric]] <- m
      
      # Predictions along focal variable
      newdata <- obs %>%
        summarise(
          S = mean(S), log_area = mean(log_area),
          latitude = mean(latitude), impact_mean = mean(impact_mean)
        ) %>%
        slice(rep(1, 100)) %>%
        mutate(!!xvar := seq(min(obs[[xvar]]), max(obs[[xvar]]), length.out = 100))
      
      newdata <- newdata %>% mutate(across(c(S, log_area, latitude, impact_mean), scale))
      
      pred <- posterior_epred(m, newdata = newdata, re_formula = NA)
      mu <- apply(pred, 2, mean)
      pi <- apply(pred, 2, quantile, probs = c(0.025, 0.975))
      
      df_line <- data.frame(
        xval = seq(min(obs[[xvar]]), max(obs[[xvar]]), length.out = 100),
        pred = mu,
        lwr = pi[1,], upr = pi[2,]
      )
      
      if (is_transformed) {
        df_line <- df_line %>% mutate(across(c(pred, lwr, upr), ~ - .x))
      }
      
      p <- ggplot(df_plot, aes(x = xval, y = mean)) +
        geom_linerange(aes(ymin = q2.5, ymax = q97.5, color = site), linewidth = 1) +
        geom_point(aes(color = site), size = 3) +
        geom_line(data = df_line, aes(x = xval, y = pred), linetype = "dashed") +
        geom_ribbon(data = df_line, aes(x = xval, ymin = lwr, ymax = upr),
                    inherit.aes = FALSE, fill = "grey70", alpha = 0.3) +
        scale_color_manual(values = site_colors, limits = df_plot$site, labels = df_plot$name) +
        labs(x = x_label, y = metric) +
        theme_bw(base_size = 14) +
        theme(legend.position = "none")
      
      plot_list[[metric]] <- p
    }

  }
  
  return(list(plots = plot_list, models = model_list))
}

plot_metric_density <- function(metaweb_metrics, network_info, metric_name) {
  library(ggplot2)
  library(dplyr)
  library(rlang)
  
  # Dynamically capture the metric column name
  metric_sym <- sym(metric_name)
  
  # Extract distribution from metaweb replicates
  df_plot <- metaweb_metrics %>% rename(site = Metaweb) %>% 
    dplyr::select(site, !!metric_sym) %>%
    rename(value = !!metric_sym)
  
  # Extract empirical values
  obs_plot <- network_info %>%
    dplyr::select(site, !!metric_sym) %>%
    rename(value = !!metric_sym)
  
  # Build density plot
  ggplot(df_plot, aes(x = value)) +
    geom_density(fill = "skyblue", alpha = 0.5, color = NA, adjust = 1.2) +
    geom_vline(
      data = obs_plot,
      aes(xintercept = value, color = site),
      linewidth = 1
    ) +
    geom_text(
      data = obs_plot,
      aes(x = value, y = 0, label = site, color = site),
      angle = 90, vjust = -0.5, hjust = 0, size = 4
    ) +
    labs(
      x = paste0(metric_name, " (Metaweb-based replicates)"),
      y = "Density",
      title = paste("Distribution of", metric_name, "across metaweb replicates"),
      subtitle = "Vertical lines and labels show empirical network values by site"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}


#' Bayesian hierarchical modeling and visualization of food-web metrics
#'
#' Fits a Student-t Bayesian hierarchical model for each network metric
#' using environmental and structural predictors (S, log_area, latitude,
#' impact_mean, and depth_m) with a random intercept for site.  
#' The function plots the empirical mean and 95% intervals (from metaweb simulations)
#' together with the model posterior mean and 95% credible interval.
#'
#' @param network_info Data frame with observed network metrics and predictors.
#' @param simulated_metrics Data frame with metaweb-based simulated metrics.
#' @param metrics Character vector of metric names to model (default: c("S","C","Rank","Entropy","Modularity","MEing_stable")).
#' @param xvar Predictor for the x-axis in plots (e.g. "latitude", "log_area", "impact_mean").
#' @param fit_model Logical, whether to fit the Bayesian model (default TRUE).
#' @param iters Number of MCMC iterations (default 2000).
#' @param chains Number of MCMC chains (default 4).
#' @param cores Number of CPU cores for parallel computation (default 4).
#' @param max_samples_per_site Maximum samples per site for model fitting (default 500).
#'
#' @return A list with fitted models, posterior samples, and plots for each metric.
#' @export
plot_metric_vs_latitude_bayes <- function(network_info,
                                          simulated_metrics,
                                          metrics = c("S", "C", "Rank", "Entropy", "Modularity", "MEing_stable"),
                                          xvar = "latitude",
                                          fit_model = TRUE,
                                          iters = 2000, chains = 4, cores = 4,
                                          max_samples_per_site = 500) {
  library(dplyr); library(ggplot2); library(tidyr); library(viridisLite)
  library(tibble); require(rethinking)
  
  if ("Metaweb" %in% names(simulated_metrics))
    simulated_metrics <- simulated_metrics %>% rename(site = Metaweb)
  
  obs <- network_info %>%
    mutate(site = as.character(site),
           log_area = log(area_km2)) %>%
    dplyr::select(site, name, C, S, latitude, log_area, depth_m, impact_mean)
  
  x_label <- case_when(
    xvar == "log_area" ~ "Log Area (km²)",
    xvar == "latitude" ~ "Latitude",
    xvar == "impact_mean" ~ "Human impact",
    xvar == "depth_m" ~ "Mean Depth (m)",
    TRUE ~ xvar
  )
  
  site_order <- obs %>% arrange(latitude) %>% pull(site)
  site_colors <- setNames(viridis(length(site_order)), site_order)
  
  out <- list(plots = list(), models = list(), post = list())
  
  for (metric in metrics) {
    if (!metric %in% names(simulated_metrics)) next
    
    sim2 <- simulated_metrics
    is_transformed <- identical(metric, "MEing_stable")
    if (is_transformed) {
      sim2 <- sim2 %>% mutate(metric_tmp = - .data[[metric]])
    } else {
      sim2 <- sim2 %>% mutate(metric_tmp = .data[[metric]])
    }
    
    ci_data <- sim2 %>%
      filter(!is.na(metric_tmp)) %>%
      group_by(site) %>%
      summarise(
        q2.5 = quantile(metric_tmp, 0.025, na.rm = TRUE),
        q97.5 = quantile(metric_tmp, 0.975, na.rm = TRUE),
        mean = mean(metric_tmp, na.rm = TRUE),
        .groups = "drop"
      )
    if (is_transformed)
      ci_data <- ci_data %>% mutate(across(c(q2.5, q97.5, mean), ~ - .x))
    
    df_plot <- obs %>%
      left_join(ci_data, by = "site") %>%
      mutate(color = site_colors[site],
             xval = case_when(
               xvar == "log_area" ~ log_area,
               xvar == "latitude" ~ latitude,
               xvar == "impact_mean" ~ impact_mean,
               xvar == "depth_m" ~ depth_m,
               TRUE ~ NA_real_
             ))
    
    df_test <- sim2 %>%
      dplyr::select(site, value = metric_tmp) %>%
      left_join(obs, by = "site") %>%
      drop_na(value, S, latitude, log_area, impact_mean, depth_m) %>%
      mutate(site = as.factor(site))
    
    if (nrow(df_test) > max_samples_per_site * length(unique(df_test$site))) {
      df_test <- df_test %>% group_by(site) %>% slice_sample(n = max_samples_per_site) %>% ungroup()
    }
    
    if (fit_model) {
      dlist <- list(
        N = nrow(df_test),
        site_id = as.integer(df_test$site),
        n_site = length(unique(df_test$site)),
        value = df_test$value,
        S = standardize(df_test$S),
        log_area = standardize(df_test$log_area),
        latitude = standardize(df_test$latitude),
        impact = standardize(df_test$impact_mean),
        depth = standardize(df_test$depth_m)
      )
      
      m <- ulam(
        alist(
          value ~ dstudent(2, mu, sigma),
          mu <- a[site_id] + bS*S + bA*log_area + bL*latitude + bI*impact + bD*depth,
          a[site_id] ~ normal(a_bar, sigma_site),
          a_bar ~ normal(0, 0.5),
          c(bS, bA, bL, bI, bD) ~ normal(0, 0.2),
          sigma_site ~ exponential(2),
          sigma ~ exponential(2)
        ),
        data = dlist, chains = chains, cores = cores, iter = iters, log_lik = TRUE
      )
      
      post <- extract.samples(m)
      out$models[[metric]] <- m
      out$post[[metric]] <- post
      
      # Build standardized newdata for prediction
      x_seq <- seq(min(df_plot$xval, na.rm = TRUE), max(df_plot$xval, na.rm = TRUE), length.out = 100)
      dnew <- list(
        S = rep(mean(dlist$S), 100),
        log_area = rep(mean(dlist$log_area), 100),
        latitude = rep(mean(dlist$latitude), 100),
        impact = rep(mean(dlist$impact), 100),
        depth = rep(mean(dlist$depth), 100)
      )
      dnew[[xvar]] <- seq(-2, 2, length.out = 100)
      
      npost <- length(post$a_bar)
      mu_post <- matrix(NA, nrow = npost, ncol = length(dnew[[xvar]]))
      for (i in seq_len(npost)) {
        mu_post[i, ] <- post$a_bar[i] +
          post$bS[i]*dnew$S + post$bA[i]*dnew$log_area +
          post$bL[i]*dnew$latitude + post$bI[i]*dnew$impact + post$bD[i]*dnew$depth
      }
      mu_mean <- colMeans(mu_post)
      mu_lwr <- apply(mu_post, 2, quantile, 0.025)
      mu_upr <- apply(mu_post, 2, quantile, 0.975)
      if (is_transformed) {
        mu_mean <- -mu_mean; mu_lwr <- -mu_lwr; mu_upr <- -mu_upr
      }
      
      df_line <- tibble(xval = x_seq, pred = mu_mean, lwr = mu_lwr, upr = mu_upr)
      
      p <- ggplot(df_plot, aes(x = xval, y = mean)) +
        geom_linerange(aes(ymin = q2.5, ymax = q97.5, color = site), linewidth = 1) +
        geom_point(aes(color = site), size = 3) +
        geom_line(data = df_line, aes(x = xval, y = pred), linetype = "dashed") +
        geom_ribbon(data = df_line, aes(x = xval, ymin = lwr, ymax = upr),
                    fill = "grey70", alpha = 0.3, inherit.aes = FALSE) +
        scale_color_manual(values = site_colors) +
        labs(x = x_label, y = metric, title = metric) +
        theme_bw(base_size = 14) + theme(legend.position = "none")
      
      out$plots[[metric]] <- p
    } else {
      p <- ggplot(df_plot, aes(x = xval, y = mean)) +
        geom_linerange(aes(ymin = q2.5, ymax = q97.5, color = site), linewidth = 1) +
        geom_point(aes(color = site), size = 3) +
        scale_color_manual(values = site_colors) +
        labs(x = x_label, y = metric, title = metric) +
        theme_bw(base_size = 14) + theme(legend.position = "none")
      
      out$plots[[metric]] <- p
    }
  }
  
  return(out)
}

# Use stan_glmer from rstanarm for Bayesian hierarchical modeling
#
plot_metric_vs_latitude_bayes1 <- function(network_info,
                                          simulated_metrics,
                                          metrics = c("S", "C", "Rank", "Entropy", "Modularity", "MEing_stable"),
                                          xvar = "latitude",
                                          fit_model = TRUE,
                                          iters = 2000, chains = 4, cores = 4,
                                          max_samples_per_site = 500) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(viridisLite)
  library(rstanarm)
  #library(arm)
  
  if ("Metaweb" %in% names(simulated_metrics))
    simulated_metrics <- simulated_metrics %>% rename(site = Metaweb)
  
  obs <- network_info %>%
    mutate(site = as.character(site),
           log_area = log(area_km2)) %>%
    dplyr::select(site, name, C, S, latitude, log_area, depth_m, impact_mean)
  
  x_label <- case_when(
    xvar == "log_area" ~ "Log Area (km²)",
    xvar == "latitude" ~ "Latitude",
    xvar == "impact_mean" ~ "Human impact",
    xvar == "depth_m" ~ "Depth (m)",
    TRUE ~ xvar
  )
  
  site_order <- obs %>% arrange(latitude) %>% pull(site)
  site_colors <- setNames(viridis(length(site_order)), site_order)
  
  out <- list(plots = list(), models = list(), post = list())
  
  for (metric in metrics) {
    if (!metric %in% names(simulated_metrics)) next
    
    sim2 <- simulated_metrics
    is_transformed <- identical(metric, "MEing_stable")
    if (is_transformed)
      sim2 <- sim2 %>% mutate(metric_tmp = - .data[[metric]])
    else
      sim2 <- sim2 %>% mutate(metric_tmp = .data[[metric]])
    
    ci_data <- sim2 %>%
      filter(!is.na(metric_tmp)) %>%
      group_by(site) %>%
      summarise(q2.5 = quantile(metric_tmp, 0.025),
                q97.5 = quantile(metric_tmp, 0.975),
                mean = mean(metric_tmp),
                .groups = "drop")
    
    if (is_transformed)
      ci_data <- ci_data %>% mutate(across(c(q2.5, q97.5, mean), ~ - .x))
    
    df_plot <- obs %>%
      left_join(ci_data, by = "site") %>%
      mutate(color = site_colors[site],
             xval = case_when(
               xvar == "log_area" ~ log_area,
               xvar == "latitude" ~ latitude,
               xvar == "impact_mean" ~ impact_mean,
               xvar == "depth_m" ~ depth_m,
               TRUE ~ NA_real_
             ))
    
    df_test <- sim2 %>%
      dplyr::select(site, value = metric_tmp) %>%
      left_join(obs, by = "site") %>%
      drop_na(value, S, latitude, log_area, impact_mean, depth_m) %>%
      mutate(site = factor(site))
    
    # Optional subsampling
    if (nrow(df_test) > max_samples_per_site * length(unique(df_test$site))) {
      df_test <- df_test %>%
        group_by(site) %>%
        slice_sample(n = max_samples_per_site) %>%
        ungroup()
    }
    
    if (fit_model) {
      df_std <- df_test %>%
        mutate(
          S_z = standardize(S),
          log_area_z = standardize(log_area),
          latitude_z = standardize(latitude),
          impact_z = standardize(impact_mean),
          depth_z = standardize(depth_m)
        )
      
      formula <- as.formula(
        "value ~ S_z + log_area_z + latitude_z + impact_z + depth_z + (1 | site)"
      )
      
      m <- stan_glmer(
        formula,
        data = df_std,
        family = gaussian(),#(link= "log"),
        chains = chains, iter = iters, cores = cores,
        prior = normal(0, 0.5),
        prior_intercept = normal(0, 1),
        prior_aux = exponential(2)
      )
      
      out$models[[metric]] <- m
      
      # Prediction line
      newdata <- tibble(
        S_z = mean(df_std$S_z),
        log_area_z = mean(df_std$log_area_z),
        latitude_z = mean(df_std$latitude_z),
        impact_z = mean(df_std$impact_z),
        depth_z = mean(df_std$depth_z)
      )
      x_seq <- seq(min(df_plot$xval, na.rm = TRUE),
                   max(df_plot$xval, na.rm = TRUE),
                   length.out = 100)
      newdata <- newdata[rep(1, length(x_seq)), ]
      newdata[[paste0(xvar, "_z")]] <- standardize(x_seq)
      
      mu_post <- posterior_epred(m, newdata = newdata)
      mu_mean <- colMeans(mu_post)
      mu_lwr <- apply(mu_post, 2, quantile, 0.025)
      mu_upr <- apply(mu_post, 2, quantile, 0.975)
      if (is_transformed) {
        mu_mean <- -mu_mean; mu_lwr <- -mu_lwr; mu_upr <- -mu_upr
      }
      
      df_line <- tibble(xval = x_seq, pred = mu_mean, lwr = mu_lwr, upr = mu_upr)
      
      p <- ggplot(df_plot, aes(x = xval, y = mean)) +
        geom_linerange(aes(ymin = q2.5, ymax = q97.5, color = site), linewidth = 1) +
        geom_point(aes(color = site), size = 3) +
        geom_line(data = df_line, aes(x = xval, y = pred), linetype = "dashed") +
        geom_ribbon(data = df_line, aes(x = xval, ymin = lwr, ymax = upr),
                    fill = "grey70", alpha = 0.3, inherit.aes = FALSE) +
        scale_color_manual(values = site_colors) +
        labs(x = x_label, y = metric, title = metric) +
        theme_bw(base_size = 14) +
        theme(legend.position = "none")
      
      out$plots[[metric]] <- p
    } else {
      p <- ggplot(df_plot, aes(x = xval, y = mean)) +
        geom_linerange(aes(ymin = q2.5, ymax = q97.5, color = site), linewidth = 1) +
        geom_point(aes(color = site), size = 3) +
        scale_color_manual(values = site_colors) +
        labs(x = x_label, y = metric, title = metric) +
        theme_bw(base_size = 14) +
        theme(legend.position = "none")
      
      out$plots[[metric]] <- p
    }
  }
  
  return(out)
}
