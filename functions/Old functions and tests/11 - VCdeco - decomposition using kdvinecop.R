##########################################################################
#### 12 - VCdeco -  decomposition using kdvinecop
##########################################################################

##########################################################################
### Function to read out edges
get_edges <- function(fit){
  trees <- length(fit)-3
  d <- nrow(fit$matrix)
  edges <- NULL
  pos_tree <- NULL
  pos_edge <- NULL
  for(i in 1:trees){
    for(j in 1:(d-i)){
      edges <- c(edges,fit[[i]][[j]][[2]])  
      pos_tree <- c(pos_tree,i)
      pos_edge <- c(pos_edge,j)
    }
  }
  return(list(name=edges,tree=pos_tree,edge=pos_edge))
}


##########################################################################
### Function to get all copulas
get_copulas <- function(fit){
  trees <- sum(is.element(names(fit),paste0("T",1:100)))
  d <- nrow(fit$matrix)
  copulas <- vector(mode="list",sum(1:d-1))
  k=1
  for(i in 1:trees){
    for(j in 1:(d-i)){
      copulas[[k]] <- fit[[i]][[j]]$c
      k=k+1
    }
  }
  copulas
}


##########################################################################
### Simulate earnings
incSim <- function(data0, data1,
                   formula, jit,
                   usim,
                   group.marginals, #Marginals of which group (0/1) to be used?
                   tau=seq(5,95,5)/100,
                   multiply=1                # multiply simulated incomes by factor 
){
  
  # Retrieve input from VCopObj
  varnames <- VCopObj$varnames
  n <- nrow(usim)
  d <- ncol(usim)
  
  # Compute weighted quantiles
  dataSim <- lapply(1:d, function(x){
                            if(group.marginals[x]==0){
                                 wtd.quantile(data0[,x+1],  weights=data0[,"(weights)"], probs=usim[,x])
                                 }else{
                                 wtd.quantile(data1[,x+1],  weights=data1[,"(weights)"], probs=usim[,x])
                                 }
  })
  dataSim <- do.call( "cbind", dataSim)
 
  # Round jittered data
  dataSim[,jit] <- round(dataSim[,jit])
  
  # Compute income
  dataSim <- as.data.frame(cbind(rep(0,n),dataSim))
  names(dataSim) <- names(data0)[1:(d+1)]
  y <- rowSums(model.matrix(formula,dataSim))*multiply
  
  # Export quantiles and other distributional stats
  res <- list(quants=quantile(y, probs=tau),
              stats=wtd.sum(y, rep(1,n)),
              ysim=y)
  return(res)           
}

##########################################################################
### Decompose differences and export results
##########################################################################

VC_deco_sim <- function(VCfit0, VCfit1, 
                        formula, jit, 
                        data0, data1,
                        reference=1,           # if reference==1, copula of group 1 is used to analyze the impact of marginals
                        n=1000,                # number of simulated observations
                        tau=seq(5,95,5)/100,   # selected quantiles
                        multiply=(1/40),       # multiply simulated incomes by factor 
                        quasi=TRUE,            # quasi random numbers or pseudo random numbers? 
                        print_info=TRUE,       # print info?
                        copula_detail=TRUE,    # decompose detailed copula
                        export_usim=TRUE       # export simulated quasi copula data
){   
  
  # Create list of vine copula densities
  allcopulas <- c(get_copulas(VCfit0),
                  get_copulas(VCfit1))
  copula <- VCfit0
  
  # Dimension
  d <- nrow(VCfit0$matrix)
  edges0 <- get_edges(VCfit0)
  edges1 <- get_edges(VCfit1)
  #common_edges <- intersect(edges0[which(nchar(edges0)==3)],edges1[which(nchar(edges1)==3)])
  #pedges0 <- match(common_edges,edges0)
  #pedges1 <- match(common_edges,edges1)
  if(reference==1){
    edgesref <- edges1[[1]]
    edgescom <- edges0[[1]]
  }else{
    edgesref <- edges0[[1]]
    edgescom <- edges1[[1]]
  }
  # Selection only undconditional copulas (first d-1 elements)
  pos_edges <- match(edgescom[1:(d-1)],edgesref[1:(d-1)])#match(edges1[[1]][which(nchar(edges0[[1]])==3)],edges0[[1]][which(nchar(edges1[[1]])==3)])
  pos_edges <- pos_edges[which(is.na(pos_edges)==FALSE)]
  common_edges <- edgesref[pos_edges]
  
  # Create matrix to estimate counterfactual distributions  
  nedges <- sum(1:(d-1))
  ncop_detail <- ifelse(length(common_edges)==(d-1),d-1,length(common_edges))
  
  # Set reference group
  r <- ifelse(reference==1,1,0)
  c <- ifelse(reference==1,0,1)
  
  # Create estimation matrix
  group.marginals <- matrix(rep(r,d*d), ncol=d, byrow=TRUE)  
  diag(group.marginals) <- rep(c,d)
  group.marginals <- rbind(rep(r,d),
                           group.marginals,
                           matrix(rep(c,(2+ncop_detail)*d), ncol=d))
  group.copula <- matrix(rep(c,ncop_detail*nedges), ncol=nedges, byrow=TRUE)
  j=1
  for(i in pos_edges){
    group.copula[j,i] <- r
    j=j+1
  }
  group.copula <- rbind(matrix(rep(r,(2+d)*nedges), ncol=nedges),
                        rep(c,nedges),
                        group.copula)
  
  
  # Simulate distributions and compute statistics 
  quants <- NULL
  stats <- NULL
  
  # Define which distribution to estimate; estimate copula decomposition only if requested
  nseq <- ifelse(copula_detail==TRUE,nrow(group.marginals), nrow(group.marginals)-ncop_detail)
  iseq <- 1:nseq
  ysim <- matrix(NA,n,nseq)
  usim_all <- vector(mode="list",length=2+ifelse(copula_detail==TRUE,ncop_detail,0))
  k=1
  
  if(print_info==TRUE){cat("\nSimulating counterfactual distributions...")}
  for(i in iseq){
    
    if(print_info==TRUE){cat((paste0(i,"...")))}
    
    gm <- group.marginals[i,]
    gc <- group.copula[i,]
    
    # Simulate (quasi-)copula data if required data have not already been simulated
    if(i==1|i>1&identical(gc,group.copula[i-ifelse(i==1,0,1),])==FALSE){
    
      # Construct copula from different groups
      copula.select.r <-  edgesref[is.element(gc,r)]
      if(length(copula.select.r)==nedges){
        copula.select.c <-  character()
        }else{
        copula.select.c <- setdiff(edgescom,copula.select.r)
       }
      if(r==0){
        copula.select <- c(which(is.element(edges0[[1]],copula.select.r)),
                           which(is.element(edges1[[1]],copula.select.c))+nedges)
        copula.tree <-   c(edges0[[2]][which(is.element(edges0[[1]],copula.select.r))],
                           edges1[[2]][which(is.element(edges1[[1]],copula.select.c))])
        copula.edge <-  c(edges0[[3]][which(is.element(edges0[[1]],copula.select.r))],
                          edges1[[3]][which(is.element(edges1[[1]],copula.select.c))])
      }else{
        copula.select <- c(which(is.element(edges0[[1]],copula.select.c)),
                           which(is.element(edges1[[1]],copula.select.r))+nedges)
        copula.tree <-   c(edges0[[2]][which(is.element(edges0[[1]],copula.select.c))],
                           edges1[[2]][which(is.element(edges1[[1]],copula.select.r))])
        copula.edge <-  c(edges0[[3]][which(is.element(edges0[[1]],copula.select.c))],
                          edges1[[3]][which(is.element(edges1[[1]],copula.select.r))])
      }
      copulas <- allcopulas[copula.select]
      
      # Assign copulas to vine copula object
      for(j in 1:length(copula.edge)){
       copula[[copula.tree[j]]][[copula.edge[j]]]$c<- copulas[[j]]
      }
      
      # Simulate quasi-compula observations from new copula constructions
      usim <- rkdevinecop(n,copula, quasi=quasi)
      usim_all[[k]] <- usim
      k=k+1
    }
    
    quants_stats <- incSim(data0,data1,
                           formula,jit,
                           usim=usim,
                           group.marginals=gm,
                           tau=tau,
                           multiply=multiply)
    
    quants <- rbind(quants,quants_stats$quants)
    stats <- rbind(stats,quants_stats$stats)
    ysim[,i] <- quants_stats$ysim
  }
  
  # Define matrix with differences to be computed
  diffagg <- matrix(c(1,    1+d+2,
                      1,    1+d+1,
                      1+d+1,1+d+2), ncol=2, byrow=TRUE) 
  diffmarg <- matrix(c(rep(1,d),
                         2:(d+1)),  ncol=2, byrow=FALSE) 
  diffcop <- matrix(c((1+d+3):(1+d+3+ncop_detail-1),
                         rep(1+d+2,ncop_detail)),  ncol=2, byrow=FALSE) 
  diffmatrix <- rbind(diffagg,diffmarg,diffcop)
  
  if(reference==0){
    diffmatrix <- diffmatrix[,c(2,1)]
  }
  diffmatrix <- diffmatrix[iseq, ]
  
  # Compute difference between empirical and simulated "observed" distribution
  quants_diff <- log(quants[diffmatrix[,1],]) - log(quants[diffmatrix[,2],])
  stats_diff <- stats[diffmatrix[,1],] - stats[diffmatrix[,2],]
  
  # Compute differences
  quants_deco <- log(quants[diffmatrix[,1],]) - log(quants[diffmatrix[,2],])
  stats_deco <- stats[diffmatrix[,1],] - stats[diffmatrix[,2],]
  
  # Compute interactions: Marginals
  quants_deco <- rbind(quants_deco,quants_deco[2,]-colSums(quants_deco[4:7,]))
  stats_deco <- rbind(stats_deco,stats_deco[2,]-colSums(stats_deco[4:7,]))
  
  # Compute interactions: Copula
  if(copula_detail==TRUE){
    quants_deco <- rbind(quants_deco,quants_deco[3,]-colSums(quants_deco[8:13,]))
    stats_deco <- rbind(stats_deco,stats_deco[3,]-colSums(stats_deco[8:13,]))
  }
  
  # Retrieve names of variables and names of edges/copulas
  varnames <- VCopObj$varnames
  edgesnames <-  VC_arrange_var_names(VCopObj)$vcopula
  
  # Name results
  diffnames1 <- c("Observed","Marginals","Copula",paste0("M",1:4),paste0("C",1:6),"M Interaction", "C Interaction")
  diffnames2 <- c("1 Observed","2 Marginals","3 Copula",paste0(1:4," ",varnames),paste0(1:6, " ",edgesnames),"5 Interaction", "7 Interaction")
  diffcat1 <- c("Observed","Marginal","Copula",rep("Marginal",4),rep("Copula",6),"Marginal","Copula")
  diffcat2 <- c(rep("aggregate",3), rep("detailed",12))
  diffnames <- data.frame(names1=diffnames1,names2=diffnames2,cat1=diffcat1,cat2=diffcat2)
  
  if(copula_detail==FALSE){
    diffnames <- diffnames[-c(8:13,15),]
    diffnames[,1] <- as.character(diffnames[,1])
  }   
  
  ## Prepare main results for export
  quants_deco <- as.data.frame(quants_deco) 
  stats_deco  <- as.data.frame(stats_deco) 
  
  names(quants_deco) <- tau
  names(stats_deco) <- names(quants_stats$stats)
  
  quants_deco <- cbind(diffnames,quants_deco)
  stats_deco <- cbind(diffnames,stats_deco)
  
  ## Check fit, 
  # i.e. estimate statistics of observed marginals and compare them 
  # to simulated statistics
  
  outcome_variable <- as.character(VCopObj$formula)[2]
  y0 <- VCopObj$data0[,outcome_variable]
  y1 <- VCopObj$data1[,outcome_variable]
  w0 <- VCopObj$data0[,"(weights)"]
  w1 <- VCopObj$data1[,"(weights)"]
  
  quantsObs0 <- wtd.quantile(y0,weights=w0, probs=tau, na.rm=TRUE)
  quantsObs1 <- wtd.quantile(y1,weights=w1, probs=tau, na.rm=TRUE)
  statsObs0 <- wtd.sum(y0,w0)
  statsObs1 <- wtd.sum(y1,w1)
  
  # Compute difference between empirical and simulated "observed" distribution
  quants_sim_vs_obs <- data.frame(tau=tau,group1=log(quants[1,])-log(quantsObs1),group0=log(quants[7,])-log(quantsObs0))
  quants_sim_vs_obs$DiffInDiff <- quants_sim_vs_obs$group1 - quants_sim_vs_obs$group0
  stats_sim_vs_obs <- data.frame(statistic=names(statsObs1),group1=stats[1,]-statsObs1,group0=stats[7,]-statsObs0)
  
  # Export simulated y
  ysim <- as.data.frame(ysim[,c(1,7)])
  names(ysim) <- c("group1","group0")
  
  ## Export usim?
  if(export_usim==FALSE){
    usim_all <- NULL
    ysim <- NULL
  }
  
  ## Return results 
  return(list(quants_deco=quants_deco,
              stats_deco=stats_deco,
              quants=quants,
              stats=stats,
              quants_sim_vs_obs=quants_sim_vs_obs,
              stats_sim_vs_obs=stats_sim_vs_obs,
              usim_all=usim_all,
              ysim=ysim,
              par=list(reference=reference,
                       multiply=multiply,
                       n=n, 
                       quasi = quasi,
                       copula_detail = copula_detail)))
}
