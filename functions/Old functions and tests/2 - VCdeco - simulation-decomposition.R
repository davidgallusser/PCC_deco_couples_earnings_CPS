##########################################################################
#### 2 - VCdeco - simulation-decomposition
##########################################################################

#### The following functions simulate counterfactual distributions based
#### on the estimated vine copula construction and the estimated marginal
#### distributions. The counterfactuals are then used to decomposed observed
#### differences into the contribution of the marginals and the copula,
#### respectively. 

##########################################################################
#### Sample from a non-parametric copula estimate
##########################################################################

## D-Vine structure, four variables, following structure:
#T1: v1   <->   v2  <->  v3  <-> v4
#T2:  v1,v2  <->  v2,v3   <->  v3,v4
#T3:   v1,v3;v2  <->  v2,v4;v3
#             v1,v4;v2,v3

# Corresponding matrix
# Matrix <- matrix(c(1, 0, 0, 0,
#                    4, 2, 0, 0,
#                    3, 4, 3, 0,
#                    2, 3, 4, 4), byrow=TRUE, ncol=4)

## For details of simulation algorithm see Czado (2019: p. 127 and p. 134) 
## Function is analogous for u_sim <- rkdecop(1000, npedge11)

rKdeVineCop <- function(n,obj,quasi=TRUE){
  
  # Retrieve copula kernel density estimates 
  CT1E1 <- obj$copulas$CT1E1
  CT1E2 <- obj$copulas$CT1E2 
  CT1E3 <- obj$copulas$CT1E3
  CT2E1 <- obj$copulas$CT2E1
  CT2E2 <- obj$copulas$CT2E2
  CT3E1 <- obj$copulas$CT3E1
  type <- obj$type
  varnames <- obj$varnames
  varorder <- obj$order
 
  ### Step 0
  # Sample w_1...w_4 \sim iid unif(0,1)
  # Quasi-random numbers or pseudo-random numbers?
  if(quasi==TRUE){
    w <- as.data.frame(qrng::ghalton(n, d = 4, method="generalized"))
  }else{
    w <- as.data.frame(matrix(runif(n*4),ncol=4, byrow=TRUE))
  }
  
  ### Simulate from D-Vine
  if(type=="DVine"){

  # V and V^2 are 4x4 upper-triangle matrices with the following entries: 
  ## Step/column 1
  # *** v11 = v11^2 = u1 = w1  ***
  v11 = w[,1]
  v11.2 = w[,1]
  
  ## Step/column 2 
  # v22 = w2 = C_{2|1}(u2|u1)
  v22 = w[,2]
  # v12 = u2 = C^-1_{2|1}(w2|u1) = h^-1_{2|1}(v22|v11^2)
  v12 = hkdecop(cbind(v11.2,v22), CT1E1, cond.var=1, inverse=TRUE)
  # *** v12^2 = u2 ***
  v12.2 = v12
  # v22^2 = C_{1|2}(u1|u2) = h_{1|2}(u1|u2)=h_{1|2}(v11^2|v12) 
  v22.2 = hkdecop(cbind(v11.2,v12), CT1E1, cond.var=2)
  
  ## Step/column 3
  # v33 = w3 = C(u3|u1:2)
  v33 = w[,3]
  # v23 = C_{3|2}(u3|u2) = h_{3|1:2}^-1(w3|C(u1|u2)) = h^-1_{3|1:2}(v33|v22^2)
  v23 =  hkdecop(cbind(v22.2,v33), CT2E1, cond.var=1, inverse=TRUE)
  # v33^2 = C_{1|2:3}(u1|u2:3) = h_{1|3;2}(C_{1|2}(u|2)|C_(3|2)(u3|u2)) = h_{1|3;2}(v22^2|v23)
  v33.2 = hkdecop(cbind(v22.2,v23), CT2E1, cond.var=2)
  # *** v13 = u3 = C^-1_{3|1:2}(w3|u1:2) = h^-1_{3|2}(h^-1_{3|1:2}(w3|C_{1|2}(u1|u2))|u2) = h^-1_{3|2}(v23|v12.2) ***
  v13 = hkdecop(cbind(v12.2,v23), CT1E2, cond.var=1, inverse=TRUE)
  # v13^2 = v13
  v13.2 = v13
  # v23^2 = C_{2|3}(u2|u3) = h_{2|3}(v12^2|v13)
  v23.2 =  hkdecop(cbind(v12.2,v13), CT1E2, cond.var=2)
  
  ## Step/column 4
  # v44 = w4 = C_{4|1:3}(u4|u1:3)
  v44 = w[,4]
  # v34 = C_{4|2:3}(u4|u2:3) = h^-1_{4|1;2:3}(w4|C_{1|2:3(u1|u2:3))})=h^-1_{4|1;2:3}(v44|v33^2)
  v34 = hkdecop(cbind(v33.2,v44), CT3E1, cond.var=1, inverse = TRUE)
  # v24 = C_{4|3}(u4|u3) = h^-1{4|2;3}(C_{4|23}(u4|u2:3)|C{2|3}(u2|u3))=h^-1_{4|2;3}(v34|v23^2)
  v24 = hkdecop(cbind(v23.2,v34), CT2E2,  cond.var=1, inverse = TRUE) 
  # *** v14 = u4 = C^-1_{4|1:3}(w4|u1:3) )= h^-1_{4|3}(C_{4|3}(u4|u3)|u3) = h^-1_{4|3}(v24|v13^2) ***
  v14 = hkdecop(cbind(v13.2,v24), CT1E3, cond.var=1,  inverse = TRUE) 
  
  }else{
  ### Simulate from C-vine
    
  # V is 4x4 upper-triangle matrices with the following entries:  
  #v1
  v11 = w[,1]    
  #v2
  v12 = hkdecop(w[,c(1,2)], CT1E1, cond.var=1, inverse = TRUE)
  #v3
  v23 = hkdecop(w[,c(2,3)], CT2E1, cond.var=1, inverse = TRUE)
  v13 = hkdecop(cbind(w[,1],v23), CT1E2, cond.var=1, inverse = TRUE)
  #v4
  v34 = hkdecop(w[,c(3,4)], CT3E1, cond.var=1, inverse = TRUE)
  v24 = hkdecop(cbind(w[,2],v34), CT2E2, cond.var=1, inverse = TRUE)
  v14 = hkdecop(cbind(w[,1],v24), CT1E3, cond.var=1, inverse = TRUE)
  }
  
  # Finally, save simulated copula data: 
  SimObs <- as.data.frame(cbind(v11,v12,v13,v14))
  names(SimObs) <- varnames 
  SimObs <- SimObs[,order(varorder)]
  
  return(SimObs)
}

##########################################################################
### Estimate density of vine copula
##########################################################################

dKdeVineCop <- function(u,obj){
  
  # Retrieve copula kernel density estimates 
  C12 <- obj$copulas[[1]]
  C23 <- obj$copulas[[2]]  
  C34 <- obj$copulas[[3]]  
  C13.2 <- obj$copulas[[4]]  
  C24.3 <- obj$copulas[[5]]  
  C14.23 <- obj$copulas[[6]] 
  #varnam <- obj[[7]][1:4] 
  
  # If is.vector, create matrix
  if(is.vector(u)){u <- matrix(u,ncol=length(u))}
  
  ### Tree 1
  # Edge 1.1: u1,u2
  d12 <- dkdecop(u[,c(1,2)],C12)
  # Edge 1.2: u2,u3
  d23 <- dkdecop(u[,c(2,3)],C23)
  # Edge 1.3: u3,u4
  d34 <- dkdecop(u[,c(3,4)],C34)
  
  ### Tree 2
  # Edge 2.1: u1,u3; u2 
  h1.2 <- hkdecop(u[,c(1,2)], C12, cond.var=2) # Hfunc1:= h_2.1(u3|u1) and Hfunc2:= h_1.2(u1|u3) where index 1,2 
  h3.2 <- hkdecop(u[,c(2,3)], C23, cond.var=1) 
  
  d13.2 <- dkdecop(cbind(h1.2,h3.2),C13.2)
  
  # Edge 2.2: u2,u4; u3 
  h2.3 <- hkdecop(u[,c(2,3)], C23, cond.var=2)
  h4.3 <- hkdecop(u[,c(3,4)], C34, cond.var=1) 
  
  d24.3 <- dkdecop(cbind(h2.3,h4.3), C24.3)
  
  ### Tree 3
  # Edge 3.1:  u1,u4; u2, u3 
  h1.23 <- hkdecop(cbind(h1.2,h3.2), C13.2, cond.var=2) 
  h4.23 <- hkdecop(cbind(h2.3,h4.3), C24.3, cond.var=1) 
  
  d14.23 <- dkdecop(cbind(h1.23,h4.23), C14.23)
  
  ### Compute and return estimated copula density 
  d <- d24.3*d14.23*d13.2*d12*d23*d34
  return(d)
}


##########################################################################
### Simulate earnings
##########################################################################

incSim <- function(VCopObj,
                   usim,
                   group.marginals=rep(1,4), #Marginals of which group (0/1) to be used?
                   tau=seq(5,95,5)/100,
                   multiply=1                # multiply simulated incomes by factor 
                   ){
  
  # Retrieve input from VCopObj
  varnames <- VCopObj$varnames
  formula <- VCopObj$formula
  jit <- VCopObj$jit
  n <- nrow(usim)
  
  # Data sets
  data0 <- VCopObj$data0
  data1 <- VCopObj$data1
  
  # Compute weighted quantiles
  dataSim <- as.data.frame(matrix(rep(NA,n*5),ncol=5))
  names(dataSim) <-  c(names(VCopObj$data0)[1],varnames)
  dataSim[,1]  <- 0
  dataSim[,2]  <- if(group.marginals[1]==0){wtd.quantile(data0[,varnames[1]], weights=data0[,"(weights)"], probs=usim[,varnames[1]])}else{wtd.quantile(data1[,varnames[1]], weights=data1[,"(weights)"], probs=usim[,varnames[1]])}
  dataSim[,3]  <- if(group.marginals[2]==0){wtd.quantile(data0[,varnames[2]], weights=data0[,"(weights)"], probs=usim[,varnames[2]])}else{wtd.quantile(data1[,varnames[2]], weights=data1[,"(weights)"], probs=usim[,varnames[2]])}
  dataSim[,4]  <- if(group.marginals[3]==0){wtd.quantile(data0[,varnames[3]], weights=data0[,"(weights)"], probs=usim[,varnames[3]])}else{wtd.quantile(data1[,varnames[3]], weights=data1[,"(weights)"], probs=usim[,varnames[3]])}
  dataSim[,5]  <- if(group.marginals[4]==0){wtd.quantile(data0[,varnames[4]], weights=data0[,"(weights)"], probs=usim[,varnames[4]])}else{wtd.quantile(data1[,varnames[4]], weights=data1[,"(weights)"], probs=usim[,varnames[4]])}
  
  # Round jittered data
  dataSim[,jit+1] <- round(dataSim[,jit+1])
  
  # Compute income
  dataSim[,1] <- rowSums(model.matrix(formula,dataSim))*multiply

  # Export quantiles and other distributional stats
  res <- list(quants=quantile(dataSim[,1], probs=tau),
              stats=wtd.sum(dataSim[,1], rep(1,n)),
               ysim=dataSim[,1])
  return(res)           
}

##########################################################################
### Decompose differences and export results
##########################################################################

VC_deco_sim <- function(VCopObj,
                        reference=1,           # if reference==1, copula of group 1 is used to analyze the impact of marginals
                        n=1000,                # number of simulated observations
                        tau=seq(5,95,5)/100,   # selected quantiles
                        multiply=(1/40),       # multiply simulated incomes by factor 
                        quasi=TRUE,            # quasi random numbers or pseudo random numbers? 
                        print_info=TRUE,       # print info?
                        copula_detail=TRUE,    # decompose detailed copula
                        export_usim=TRUE       # export simulated quasi copula data
                        ){   
  
  # Define statistics matrix  
  r <- ifelse(reference==1,1,0)
  c <- ifelse(reference==1,0,1)
  
  group.marginals <- matrix(rep(r,4*4), ncol=4, byrow=TRUE)  
  diag(group.marginals) <- rep(c,4)
  group.marginals <- rbind(rep(r,4),
                           group.marginals,
                           matrix(rep(c,8*4), ncol=4))
  
  group.marg.inter <- matrix(c(c,c,r,r,
                               c,r,c,c,
                               c,r,r,c,
                               r,c,c,r,
                               r,c,r,c,
                               r,r,c,c), byrow=TRUE, ncol=4)
  
  group.copula <- matrix(rep(c,6*6), ncol=6, byrow=TRUE)  
  diag(group.copula) <- rep(r,6)
  group.copula <- rbind(matrix(rep(r,6*6), ncol=6),
                        rep(c,6),
                        group.copula)
  group.cop.inter <- matrix(c(r,r,c,c,c,c,
                              r,c,r,c,c,c,
                              r,c,c,r,c,c,
                              c,r,r,c,c,c,
                              c,r,c,r,c,c,
                              c,c,r,r,c,c), byrow=TRUE, ncol=6)
  
  # Simulate distributions and compute statistics 
  quants <- NULL
  stats <- NULL

  if(copula_detail==TRUE){
    iseq <- 1:13
    ysim <- matrix(NA,n,13)
  }else{
    iseq <- 1:7  
    ysim <- matrix(NA,n,7)
  }
  usim_all <- vector(mode="list",length=length(iseq)-5)
  k=1
  
  if(print_info==TRUE){cat("\nSimulating counterfactual distributions...")}
  for(i in iseq){
    
    if(print_info==TRUE){cat((paste0(i,"...")))}
    
    gm <- group.marginals[i,]
    gc <- group.copula[i,]
    
    # Simulate (quasi-)copula data if required data have not already been simulated
    if(i==1|i>1&identical(gc,group.copula[i-ifelse(i==1,0,1),])==FALSE){
      
      # Construct copula from different groups
      copula.select <- 1:6+6*gc
      copula <- c(VCopObj$copula0$copulas,VCopObj$copula1$copulas)
      copula <- copula[copula.select]
      copula <- list(copulas=copula,
                     type=VCopObj$copula0$type,
                     varnames=VCopObj$copula0$varnames,
                     order=VCopObj$copula0$order)

      # Simulate quasi-compula data 
      usim <- rKdeVineCop(n,copula,quasi=quasi)
      usim_all[[k]] <- usim
      k=k+1
    }
    
    quants_stats <- incSim(VCopObj,
                           usim=usim,
                           group.marginals=gm,
                           tau=tau,
                           multiply=multiply)
    
    quants <- rbind(quants,quants_stats$quants)
    stats <- rbind(stats,quants_stats$stats)
    ysim[,i] <- quants_stats$ysim
  }
  
  # Define matrix with differences to be computed
  if(reference==1){
    diffmatrix <- cbind(c(1,1,6,1,1,1,1,8:13),
                        c(7,6,7,2,3,4,5,rep(7,6)))
  }else{
    diffmatrix <- cbind(c(7,6,7,2,3,4,5,rep(7,6)),
                        c(1,1,6,1,1,1,1,8:13))
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


##########################################################################
### Get conditional copula of simulated and observed distribution
##########################################################################

CondCop_check_fitted <- function(VCopObj,variable,q=10){
  
  VCopObj <- VCopObj$VCfit
  deco_res <-  VCopObj$VCdeco
  
  outcome_variable <- as.character(VCopObj$formula)[2]
  y0 <- VCopObj$data0[,outcome_variable]
  y1 <- VCopObj$data1[,outcome_variable]
  w0 <- VCopObj$data0[,"(weights)"]
  w1 <- VCopObj$data1[,"(weights)"]
  
  ysim <- deco_res$ysim
  
  l0 <- share0(round(VCopObj$data0[,variable],0),w0)
  l1 <- share0(round(VCopObj$data1[,variable],0),w1)
  
  Obs0 <- CondCop(rank.wt(y0,w0), VCopObj$u0[,variable],u=l0,q=q)
  Obs1 <- CondCop(rank.wt(y1,w1), VCopObj$u1[,variable],u=l1,q=q)
  Sim0 <- CondCop(rank.wt(ysim[,2]), deco_res$usim_all[[2]][,variable],u=l0,q=q)
  Sim1 <- CondCop(rank.wt(ysim[,1]), deco_res$usim_all[[1]][,variable],u=l1,q=q)
  
  res <- data.frame(tau=rep((0:(q-1))/q,2),
                    delta=c(Sim1-Obs1,Sim0-Obs0),
                    group=rep(c("Group 1","Group 0"),each=length(Obs1))                      )
  plot <- ggplot(res, aes(tau,delta, color=group)) + 
                geom_line() + geom_point() + 
                facet_grid(~group) + geom_hline(yintercept = 0, color="lightgrey")
  return(plot)
}

##########################################################################
### Check fitted bivariate dependencies
##########################################################################

VC_check_fitted_dependencies <- function(VCopObj){
  
  deco_res <- VCopObj$VCdeco
  VCopObj <- VCopObj$VCfit

  order <- VCopObj$order
  type <- VCopObj$type
  usim_all <- deco_res$usim_all

  vn <- VCopObj$varnames[order]
  if(type=="DVine"){
    vA <- c(1,2,3,1,2,1)
    vB <- c(2,3,4,3,4,4)
    sltri <- c(2,7,12,3,8,4)
  }else{
    vA <- c(1,1,1,2,2,3)
    vB <- c(2,3,4,3,4,4)
    sltri <- which(lower.tri(matrix(1:16, ncol=4)))
  }
  vcomb <- paste0(1:6,rep(" ",6),vn[vA],rep("-",6),vn[vB],c(rep("",3),rep(" (cond.)",3)))
  
  rho0obs <- cov.wt(VCopObj$u0[,order], wt=VCopObj$data0$`(weights)`, cor=TRUE)$cor[sltri] 
  rho1obs <- cov.wt(VCopObj$u1[,order], wt=VCopObj$data1$`(weights)`, cor=TRUE)$cor[sltri] 
  
  rho0sim <- cov.wt(usim_all[[2]][,order], wt=rep(1,nrow(usim_all[[2]])), cor=TRUE)$cor[sltri]
  rho1sim <- cov.wt(usim_all[[1]][,order], wt=rep(1,nrow(usim_all[[1]])), cor=TRUE)$cor[sltri]
  
  tau0obs  <- TauMatrix(VCopObj$u0[,order], weights=VCopObj$data0$`(weights)`)[sltri] 
  tau1obs  <- TauMatrix(VCopObj$u1[,order], weights=VCopObj$data1$`(weights)`)[sltri] 
  
  tau0sim  <- TauMatrix(usim_all[[2]][,order])[sltri] 
  tau1sim  <- TauMatrix(usim_all[[1]][,order])[sltri] 
  
  res <- data.frame(var=rep(vcomb,8),
                    observed=rep(c("obs","sim","obs","sim"),each=12),
                    group=rep(rep(c("group 0","group 1"),each=6),4),
                    stat=rep(c("rho","tau"),each=24),
                    rho=c(rho0obs,rho1obs, rho0sim,rho1sim,tau0obs, tau1obs, tau0sim,tau1sim)
                    )
  return(res)
}



