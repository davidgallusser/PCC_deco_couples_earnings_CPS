##########################################################################
#### 2 - VCdeco -  decomposition using rvinecopulib
##########################################################################

##########################################################################
### Function to get all copulas
get_edge_names1 <- function(fit){
  c1 <- summary(fit)$conditioned
  c2 <- summary(fit)$conditioning
  vn <- fit$names
  names <- mapply(function(y,x) paste0(paste0(y[1],",",y[2]),
                                       ifelse(length(x)==0,
                                              "",
                                              paste0(";",paste0(x,collapse=",")))), 
                  y=c1,x=c2)
  names_s <- mapply(function(y,x) paste0(paste0(sort(y)[1],",",sort(y)[2]),
                                         ifelse(length(x)==0,
                                                "",
                                                paste0(";",paste0(sort(x),collapse=",")))), 
                    y=c1,x=c2)
  
  st <- fit$structure$struct_array
  edges <- unlist(lapply(st, function(x) 1:length(x)))
  trees <- unlist(mapply(function(y,x) rep(x,length(y)), y= st,x=1:length(st)))
  varnames <-  mapply(function(y,x) paste0(paste0(vn[sort(y)[1]],",",vn[sort(y)[2]]),
                                           ifelse(length(x)==0,
                                                  "",
                                                  paste0(";",paste0(vn[sort(x)],collapse=",")))), 
                      y=c1,x=c2)
  list(names=names,names_sorted=names_s,trees=trees,edges=edges, varnames=varnames)
}


##########################################################################
### Function to get names of edges (rvinecopulib:::get_graph() is similar)
get_edge_names <- function(fit){
    vn <- fit$names
    M <- get_matrix(fit)
    d <- get_structure(fit)$d
    t <- get_structure(fit)$trunc_lvl 
    cd <- cg <- NULL
    cdlist <- cglist <- list()
    csd <- csg <- NULL
    vnd <- vng <- NULL
    vnx <- vny <- NULL
    tree <- edge <- NULL
    ti <- 0
    k=1
    for(i in 1:t){
        e <- d-i
        ti <- ti + 1
        ei <- 0
        for(j in 1:e){
          cd1 <- c(M[d+1-j,j],M[i,j])
          if(is.element(0,cd1)==FALSE){
          cg1 <- if((i-1)<1){NA}else{M[(i-1):1,j]}
          ei <- ei + 1
          tree <- c(tree,ti)
          edge <- c(edge,ei)
          cdlist[[k]] <- cd1
          cglist[[k]] <- cg1
          k=k+1
          }
      }
    }
    
cd  <- unlist(lapply(cdlist,function(x) paste0(x,collapse=",")))
csd <- unlist(lapply(cdlist,function(x) paste0(sort(x),collapse=",")))
cg  <- unlist(lapply(cglist,function(x) paste0(x,collapse=",")))
csg <- unlist(lapply(cglist,function(x) paste0(sort(x),collapse=",")))
vnd <- unlist(lapply(cdlist,function(x) paste0(vn[sort(x)],collapse=",")))
vng <- unlist(lapply(cglist,function(x) paste0(vn[sort(x)],collapse=",")))

names_u <- names_s <- varnames <- rep(NA,length(cd))
  for(i in 1:length(cd)){
      names_u[i] <- paste0(cd[i],ifelse(is.na(cg[i])|cg[i]==""|cg[i]=="NA","",paste0(";",cg[i])))
      names_s[i] <- paste0(csd[i],ifelse(is.na(csg[i])|csg[i]=="","",paste0(";",csg[i])))
      varnames[i] <- paste0(vnd[i],ifelse(is.na(vng[i])|vng[i]=="","",paste0(";",vng[i])))
   }

list(edgesnames=names_u,edgesnames_sorted=names_s,
     tree=tree,edge=edge,
     varnames=varnames,
     vars_conditioned=cdlist,
     vars_conditioning=cglist)
}


##########################################################################
### Function to get all copulas
get_copulas_in_list <- function(fit){
  x <- get_all_pair_copulas(fit)
  d <- fit$structure$d
  copulas <- x[[1]]
  if(d>2){
  for(i in 2:length(x)){
    copulas <- c(copulas,x[[i]])
  }
  }
  copulas
}


##########################################################################
### Simulate earnings
incSim <- function(data0,data1,
                   formula, # jit, jit_bd,
                   usim,
                   group.marginals, #Marginals of which group (0/1) to be used?
                   tau=seq(5,95,5)/100,
                   multiply=1                # multiply simulated incomes by factor 
){
  
  # Retrieve dimensions
  n <- nrow(usim)
  d <- ncol(usim)
  
  # Compute weighted quantiles
  dataSim <- lapply(1:d, function(x){
                            if(group.marginals[x]==0){
                                 wtd.quantile(data0[,x+1],
                                              weights=data0[,"(weights)"],
                                              probs=usim[,x])
                                 }else{
                                 wtd.quantile(data1[,x+1],
                                              weights=data1[,"(weights)"],
                                              probs=usim[,x])
                                 }
  })
  dataSim <- do.call( "cbind", dataSim)
 
  # Round jittered data
  #dataSim[,jit] <- round(dataSim[,jit])
  #round_jittered <- floor(abs(log(jit_bd,base=10)))
  #dataSim[,jit] <- do.call("cbind",lapply(jit,function(x) round(dataSim[,x],round_jittered[x])))
  
  # Simulate income
  dataSim <- as.data.frame(cbind(rep(0,n),dataSim))
  names(dataSim) <- names(data0)[1:(d+1)]
  y <- rowSums(model.matrix(formula,dataSim))*multiply
  
  # Export quantiles and other distributional stats
  res <- list(quants=quantile(y, probs=tau),
              stats=wtd.sum(y, rep(1,n), log_trans=TRUE),
              ysim=y,
              xsim=dataSim)
  return(res)           
}

##########################################################################
### Decompose differences and export results
##########################################################################

VC_deco_sim <- function(VCfit0, VCfit1,
                        formula, #jit, jit_bd,
                        data0, data1,
                        reference=1,           # if reference==1, copula of group 1 is used to analyze the impact of marginals
                        n=1000,                # number of simulated observations
                        tau=seq(5,95,5)/100,   # selected quantiles
                        multiply=(1/40),       # multiply simulated incomes by factor 
                        quasi=TRUE,            # quasi random numbers or pseudo random numbers? 
                        detailed_deco=TRUE,    # detailed decomposition?  (i.e. isolate influence of single pair copulas?)
                        print_info=TRUE,       # print info?
                        export_usim=TRUE,      # export simulated quasi copula data
                        import_usim=list(),    # import simulated quasi copula data
                        cores=cores
){   
  
  # Create list of vine copula densities
  copulas <- list(VCfit0,VCfit1)
  allcopulas_as_list <- list(get_copulas_in_list(VCfit0),
                             get_copulas_in_list(VCfit1))

  
  # Dimension
  d <- VCfit0$structure$d
  edges0 <- get_edge_names(VCfit0)
  edges1 <- get_edge_names(VCfit1)
  #common_edges <- intersect(edges0[which(nchar(edges0)==3)],edges1[which(nchar(edges1)==3)])
  #pedges0 <- match(common_edges,edges0)
  #pedges1 <- match(common_edges,edges1)
  if(reference==1){
    edgesref <- edges1[[1]]
    edgesref_s <- edges1[[2]]
    edgescom <- edges0[[1]]
    edgescom_s <- edges0[[2]]
    edgesref_all <- edges1
    edgescom_all <- edges0
  }else{
    edgesref <- edges0[[1]]
    edgesref_s <- edges0[[2]]
    edgescom <- edges1[[1]]
    edgescom_s <- edges1[[2]]
    edgesref_all <- edges0
    edgescom_all <- edges1
  }
  # Selection only undconditional copulas (first d-1 elements)
  #sel <- 1:(d-1)
  sel <- 1:length(edgescom)
  pos_edges_all <- na.exclude(match(edgescom_s,edgesref_s))
  common_edges <- edgesref[pos_edges_all]
  common_edges_vn <-  edgesref_all[[5]][pos_edges_all]
  #permuted_edgesref <- edgesref!=edgesref_s
  #permuted_common_edges <- permuted_edgesref[pos_edges_all]
  
  #sel <- 1:length(edgescom)
  #pos_edges <- na.exclude(match(edgescom[sel],edgesref[sel]))#match(edges1[[1]][which(nchar(edges0[[1]])==3)],edges0[[1]][which(nchar(edges1[[1]])==3)])
  #pos_edges_s <- na.exclude(match(edgescom_s[sel],edgesref[sel]))
  #pos_edges_all <- unique(c(pos_edges,pos_edges_s))
  #permuted_common_edges <- c(rep(FALSE,length(pos_edges)),rep(TRUE,length(setdiff(pos_edges_s,pos_edges))))
  #common_edges <- edgesref[pos_edges_all]
  #common_edges_vn <- edgesref_all[[5]][pos_edges_all]
 
  # Create matrix to estimate counterfactual distributions  
  nedges <- sum(1:(d-1))
  n_common_copulas <- length(common_edges)#min(c(d-1,length(common_edges)))
  n_detailed_copulas <- ifelse(nedges>1,n_common_copulas+sum(1:(n_common_copulas-1)),1)
  n_detailed_marginals <- d+nedges
  
  # Set reference group
  r <- ifelse(reference==1,1,0)
  c <- ifelse(reference==1,0,1)
  
  # Create estimation matrix
  group.marginals <- matrix(rep(r,d*n_detailed_marginals), ncol=d, byrow=TRUE)  
  diag(group.marginals) <- rep(c,d)
  k=d+1
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      group.marginals[k,c(i,j)] <- rep(c,2)
      k=k+1
    }
  }
  group.marginals <- rbind(rep(r,d),
                           group.marginals,
                           matrix(rep(c,(2+n_detailed_copulas)*d), ncol=d))
  
  group.copula <- matrix(rep(c,n_detailed_copulas*nedges),
                         ncol=nedges, byrow=TRUE)
  j=1
  for(i in pos_edges_all){
    group.copula[j,i] <- r
    j=j+1
  }
  if(n_common_copulas>1){
      k=n_common_copulas+1
      for(i in 1:(n_common_copulas-1)){
        for(j in (i+1):n_common_copulas){
         group.copula[k,pos_edges_all[c(i,j)]] <- rep(r,2)
          k=k+1
        }
      }
  }
  group.copula <- rbind(matrix(rep(r,(2+n_detailed_marginals)*nedges), ncol=nedges),
                        rep(c,nedges),
                        group.copula)
  
  #group.copula.permuted <-  matrix(FALSE,nrow(group.copula),ncol(group.copula))
  #for(j in which(permuted_common_edges==TRUE)){
  #  sel <- setdiff(which(group.copula[,pos_edges_all[j]]==r),1:(n_detailed_marginals+2))
  #  group.copula.permuted[sel,pos_edges_all[j]] <- TRUE
  #  #for(i in pos_edges_s){
  #  #      group.copula.permuted[d+3+j,i] <- TRUE
  #  #      }
  #}  
  
  ### Simulate distributions and compute statistics 
  # Only aggregate decomposition?
  if(detailed_deco){
      iseq <- 1:nrow(group.marginals)
  }else{
      iseq <- c(1,1+d+nedges+1:2)
  }
  ysim <- matrix(NA,n,length(iseq))
  xsim <- vector(mode="list",length=length(iseq))
  usim_all <- vector(mode="list",length=2+ifelse(detailed_deco,n_common_copulas,0))
  k=1
  l=1
  
  quants <- NULL
  stats <- NULL
  
  if(print_info==TRUE){cat("\nSimulating counterfactual distributions...")}
  for(i in iseq){
    
    if(print_info==TRUE){cat((paste0(i,"...")))}
    
    gm <- group.marginals[i,]
    gc <- group.copula[i,]
    #gc_p <- group.copula.permuted[i,]
    
    # Simulate (quasi-)copula data if required data have not already been simulated
    if(i==1|i>1&identical(gc,group.copula[i-ifelse(i==1,0,1),])==FALSE){
    
      ## Construct counterfactual copulas 
      # Which pair-copulas of reference group are needed?
      copula.select.r <- edgesref_s[is.element(gc,r)] # selection of sorted names
      copula.select.rr <- edgesref[is.element(gc,r)]  # selection of actual names
      
      if(length(copula.select.r)==nedges){
        # If all edges of reference copula are selected,
        # select entire vine copula of reference group
        cg <- r+1
        sc_select <- numeric()   #selection of pair-copulas of comparison group
        sr_select <- 1:nedges    #selection of pair-copulas of reference group
        tree_replace <- edgesref_all[[3]] #tree of pair-copulas
        edge_replace <- edgesref_all[[4]] #edge index of pair-copulas
        permuted <- rep(FALSE,nedges) #permute pair-copula
        }else{
        # If not all edges of reference copula are selected,
        # select vine copula of comparison group
        cg <- c+1
      
        # Which pair-copulas of comparison group are used?
        copula.select.c <- setdiff(edgescom_s,copula.select.r)
        
        # Position of selected pair-copulas from  reference and 
        # comparison group in vine copula of comparison group
        sc_replace <- match(copula.select.c,edgescom_s)
        sr_replace <- match(copula.select.r,edgescom_s)
     
        # Position of selected pair-copulas respective vine copula
        sc_select <- which(is.element(edgescom_s,copula.select.c))
        sr_select <- which(is.element(edgesref_s,copula.select.r))
        
        # Trees and edge index of pair-copulas in vine copula of comparison group
        # where selected pair-copulas are "plug in".
        tree_replace <- c(edgescom_all[[3]][sc_replace],edgescom_all[[3]][sr_replace])
        edge_replace <- c(edgescom_all[[4]][sc_replace],edgescom_all[[4]][sr_replace])
        
        # Do selected pair-copulas need permutation before "plugging-in"?
        copula.replace.c <- edgescom[sr_replace]
        permuted <- c(rep(FALSE,length(copula.select.c)),copula.replace.c!=copula.select.rr)
        }
      
        # Contruct copulas
        pair_copulas_replace <- c(allcopulas_as_list[[c+1]][sc_select],
                                  allcopulas_as_list[[r+1]][sr_select])
        
        # Assign copulas to vine copula object
        for(j in 1:length(edge_replace)){
          pair_copula <- pair_copulas_replace[[j]]
          if(permuted[j]){pair_copula$parameters <- t(pair_copula$parameters)}
          copulas[[cg]]$pair_copulas[[tree_replace[j]]][[edge_replace[j]]] <- pair_copula
        }
        
        # Simulate quasi-compula observations from new copula constructions
        #usim <- rvinecop(n, copulas[[cg]], qrng = quasi, cores = cores)
        usim <- rvinecopulib:::vinecop_sim_cpp(copulas[[cg]], n, quasi, cores, rvinecopulib:::get_seeds())
        colnames(usim) <- copulas[[cg]]$names[1:d]
        usim_all[[k]] <- usim
        gc_name <- paste(gc, collapse="-")
        names(usim_all)[k] <- gc_name
        k=k+1
     
    }
    
    # Use imported simulated copula data if required
    if(length(import_usim)!=0&gc_name%in%names(import_usim)){
      usim <- import_usim[[gc_name]]
    }
    
    quants_stats <- incSim(data0,data1,
                           formula, #jit, jit_bd,
                           usim=usim,
                           group.marginals=gm,
                           tau=tau,
                           multiply=multiply)
    
    quants <- rbind(quants,quants_stats$quants)
    stats <- rbind(stats,quants_stats$stats)
    ysim[,l] <- quants_stats$ysim
    xsim[[l]] <- quants_stats$xsim
    l=l+1
  }
  
  ## Define matrix with differences to be computed
  # Direct differences
  if(detailed_deco){
    diffagg <- matrix(c(1,    1+n_detailed_marginals+2,
                         1,    1+n_detailed_marginals+1,
                         1+n_detailed_marginals+1,1+n_detailed_marginals+2), ncol=2, byrow=TRUE) 
    diffmarg <- matrix(c(rep(1,d),
                       2:(d+1)),  ncol=2, byrow=FALSE) 
    diffcop <- matrix(c((1+n_detailed_marginals+3):(1+n_detailed_marginals+3+n_common_copulas-1),
                      rep(1+n_detailed_marginals+2,n_common_copulas)),  ncol=2, byrow=FALSE) 
    diffmatrix <- rbind(diffagg,diffmarg,diffcop)
  }else{
    diffmatrix <- diffagg <- matrix(c(1,3,1,2,2,3),3,2, byrow=TRUE)
  }
  if(reference==0){diffmatrix <- diffmatrix[,2:1]}

  ## Compute direct effects
  quants_deco <- log(quants[diffmatrix[,1],]) - log(quants[diffmatrix[,2],])
  stats_deco <- stats[diffmatrix[,1],] - stats[diffmatrix[,2],]
  
  ## Compute detailed deco
  if(detailed_deco){
  ## Compute first order interactions: Marginals
  diffinteractions_marg <-  matrix(1,nedges,4)
  k=1
  for(i in 1:(d-1)){
    for(j in (i+1):d){
      diffinteractions_marg[k,2:4] <- c(i,j,d+k)+1
      k=k+1
    }
  }
  if(reference==0){diffinteractions_marg <-  diffinteractions_marg[, 4:1]}
  quants_deco <- rbind(quants_deco,
                       log(quants[diffinteractions_marg[,1],]) - log(quants[diffinteractions_marg[,4],])
                       - log(quants[diffinteractions_marg[,1],]) + log(quants[diffinteractions_marg[,2],])
                       - log(quants[diffinteractions_marg[,1],]) + log(quants[diffinteractions_marg[,3],]))
  stats_deco <- rbind(stats_deco,
                      stats[diffinteractions_marg[,1],] - stats[diffinteractions_marg[,4],]
                      - stats[diffinteractions_marg[,1],] + stats[diffinteractions_marg[,2],]
                      - stats[diffinteractions_marg[,1],] + stats[diffinteractions_marg[,3],])
  
  ## Compute first order interactions of copula if more than one copula
  if(n_common_copulas>1){
    diffinteractions_cop <-  matrix(1+n_detailed_marginals+2,sum(1:(n_common_copulas-1)),4)
    k=1
    for(i in 1:(n_common_copulas-1)){
      for(j in (i+1):n_common_copulas){
        diffinteractions_cop[k,3:1] <- c(i,j,n_common_copulas+k)+(1+n_detailed_marginals+2)
        k=k+1
      }
    }
    if(reference==0){diffinteractions_cop <-  diffinteractions_cop[, 4:1]}
    quants_deco <- rbind(quants_deco,
                         log(quants[diffinteractions_cop[,1],]) - log(quants[diffinteractions_cop[,4],])
                         -log(quants[diffinteractions_cop[,2],]) + log(quants[diffinteractions_cop[,4],])
                         -log(quants[diffinteractions_cop[,3],]) + log(quants[diffinteractions_cop[,4],]))
    stats_deco <- rbind(stats_deco,
                        stats[diffinteractions_cop[,1],] - stats[diffinteractions_cop[,4],] 
                        -stats[diffinteractions_cop[,2],] + stats[diffinteractions_cop[,4],] 
                        -stats[diffinteractions_cop[,3],] + stats[diffinteractions_cop[,4],])
  }
  
  ## Compute higher order interactions
  # Compute interactions order > 1: Marginals
  sel <- c(4:(3+d),(3+d+n_common_copulas+1):(3+d+n_common_copulas+nedges))
  quants_deco <- rbind(quants_deco,quants_deco[2,]-colSums(quants_deco[sel,]))
  stats_deco <- rbind(stats_deco,stats_deco[2,]-colSums(stats_deco[sel,]))
  
  # Compute sum of interactions order > 0: Marginals
  sel <- 4:(3+d)
  quants_deco <- rbind(quants_deco,quants_deco[2,]-colSums(quants_deco[sel,]))
  stats_deco <- rbind(stats_deco,stats_deco[2,]-colSums(stats_deco[sel,]))

  ## Copula higher order interactions if available
  if(n_common_copulas>1){
  # Compute interactions order > 1: Copula
  sel <- c((3+d+1):(3+d+n_common_copulas),(3+d+n_common_copulas+nedges+1):(3+d+n_common_copulas+nedges+sum(1:(n_common_copulas-1))))
  quants_deco <- rbind(quants_deco,
                       quants_deco[3,]-colSums(quants_deco[sel,]))
  stats_deco <- rbind(stats_deco,
                      stats_deco[3,]-colSums(stats_deco[sel,]))
    
  # Compute sum interactions order > 0: Copula
  sel <- (3+d+1):(3+d+n_common_copulas)
    quants_deco <- rbind(quants_deco,
                         quants_deco[3,]-colSums(quants_deco[sel,]))
    stats_deco <- rbind(stats_deco,
                        stats_deco[3,]-colSums(stats_deco[sel,]))
  }else{
    quants_deco <- rbind(quants_deco,quants_deco[3,]-quants_deco[(3+d+1),],0)
    stats_deco <- rbind(stats_deco,stats_deco[3,]-stats_deco[(3+d+1),],0)
  }
 
  } #end detailed deco
  
  # Retrieve names of variables and names of edges/copulas
  varnames <- colnames(usim)
  
  # Name results
  diffnames1 <- c("Total difference","Marginals","Copula")
  diffnames2 <- c("1 Total difference","2 Marginals","3 Copula")
  diffcat1 <- c("Observed","Marginal","Copula")
  diffcat2 <- c(rep("aggregate",3))
  
  if(detailed_deco){
  combM <- combn(d,2)
  if(n_common_copulas>1){combC <- combn(n_common_copulas,2)}
  diffnames1 <- c(diffnames1,
                  paste0("M ",1:d),
                  paste0("C ",1:n_common_copulas),
                  paste0("M ",combM[1,],"&",combM[2,]))
  if(n_common_copulas>1){diffnames1 <- c(diffnames1,paste0("C ",combC[1,],"&",combC[2,]))}
  diffnames1  <- c(diffnames1, "M Interactions > 1", "M Interactions > 0", "C Interactions > 1", "C Interactions > 0")
  diffnames2 <- c(diffnames2,
                  paste0(1:d," ", varnames),
                  paste0(1:n_common_copulas, " ",common_edges_vn),
                  paste0((d+1):(d+ncol(combM))," ",varnames[combM[1,]]," & ",varnames[combM[2,]]))
  if(n_common_copulas>1){diffnames2 <- c(diffnames2,paste0((n_common_copulas+1):(n_common_copulas+ncol(combC))," ",
                                                           common_edges_vn[combC[1,]]," & ",common_edges_vn[combC[2,]]))}
  diffnames2 <- c(diffnames2, 
                  paste0((d+ncol(combM)+1)," M Interactions > 1"),
                  paste0((d+1)," M Interactions > 0"),
                  paste0(ifelse(n_common_copulas>1,(n_common_copulas+ncol(combC)+1),2)," C Interactions > 1"),
                  paste0(ifelse(n_common_copulas>1,(n_common_copulas+1),2)," C Interactions > 0"))
  
  diffcat1 <-   c(diffcat1,
                  rep("Marginal",d),
                  rep("Copula",n_common_copulas),
                  rep("Marginal", ncol(combM)),
                  rep("Copula",ifelse(n_common_copulas>1,ncol(combC),0)),
                  "Marginal","Marginal", "Copula", "Copula")
  diffcat2 <- c(diffcat2,
                rep("direct",d+n_common_copulas),
                rep("interaction",nedges+sum(1:(n_common_copulas-1))),
                c("interaction","direct","interaction","direct"))
  
  }
  diffnames <- data.frame(names1=diffnames1,names2=diffnames2,cat1=diffcat1,cat2=diffcat2)
  
  #diffnames1 <- c("Observed","Marginals","Copula",paste0("M",1:d),paste0("C",1:n_common_copulas),"M Interaction", "C Interaction")
  #diffnames2 <- c("1 Observed","2 Marginals","3 Copula",paste0(1:d," ",varnames),paste0(1:n_common_copulas, " ",common_edges_vn),paste0(d+1," Interaction"), paste0(n_common_copulas+1," Interaction"))
  #diffcat1 <- c("Observed","Marginal","Copula",rep("Marginal",d),rep("Copula",n_common_copulas),"Marginal","Copula")
  #diffcat2 <- c(rep("aggregate",3), rep("detailed",nrow(group.marginals)-1))

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
  y0 <- data0[,1]
  y1 <- data1[,1]
  w0 <- data0[,"(weights)"]
  w1 <- data1[,"(weights)"]
  quantsObs0 <- wtd.quantile(y0,weights=w0, probs=tau, na.rm=TRUE)
  quantsObs1 <- wtd.quantile(y1,weights=w1, probs=tau, na.rm=TRUE)
  statsObs0 <- wtd.sum(y0,w0,log_trans=TRUE)
  statsObs1 <- wtd.sum(y1,w1,log_trans=TRUE)
  
  # Compute differences between empirical and simulated distribution
  quants_sim_vs_obs <- data.frame(tau=tau,
                                  group1=log(quants[diffagg[1,1],])-log(quantsObs1),
                                  group0=log(quants[diffagg[1,2],])-log(quantsObs0)
                                  )
  quants_sim_vs_obs$DiffInDiff <- quants_sim_vs_obs$group1 - quants_sim_vs_obs$group0
  stats_sim_vs_obs <- data.frame(statistic=names(statsObs1),
                                 group1=stats[diffagg[1,1],]-statsObs1,
                                 group0=stats[diffagg[1,2],]-statsObs0)
  rmse <- c(group0=sqrt(mean(sum(quants_sim_vs_obs$group0)^2)),
            group1=sqrt(mean(sum(quants_sim_vs_obs$group1)^2)))
  
  # Export simulated y
  sel <- c(diffagg[1,1],diffagg[1,2])
  ysim <- as.data.frame(ysim[,sel])
  xsim <- xsim[sel]   
  names(ysim) <- names(xsim) <- c("group1","group0")
  
  ## Export usim?
  if(export_usim==FALSE){
    usim_all <- NULL
    ysim <- NULL
    xsim <- NULL
  }
  
  ## Return results 
  return(list(quants_deco=quants_deco,
              stats_deco=stats_deco,
              quants=quants,
              stats=stats,
              quants_sim_vs_obs=quants_sim_vs_obs,
              stats_sim_vs_obs=stats_sim_vs_obs,
              rmse=rmse,
              usim_all=usim_all,
              ysim=ysim,
              xsim=xsim,
              par=list(reference=reference,
                       multiply=multiply,
                       tau=tau,
                       n=n, 
                       quasi = quasi,
                       detailed_deco = detailed_deco)))
}
