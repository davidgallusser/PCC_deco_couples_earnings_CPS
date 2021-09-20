##########################################################################
#### 25 - VCdeco - bootstrap inference using rvinecopulib
##########################################################################

#### The following code offers a bootstrap function for the vine copula
#### decomposition
require(parallel)
require(pbapply)

##########################################################################
VC_deco_bs <- function(VCfit, deco_res, it=100, 
                       cores=1, psu=NULL, stratum=NULL, alpha=0.05){
  # Store simulation parameter
  par <- deco_res$aggregate$par

  # Resample 
  s <- VC_resample(VCfit, it=it, psu=psu, stratum=stratum)

  # Reestimate (parallelized if requested)
  if(cores==1){  
    #res <- lapply(s, function(x) VC_deco_for_bs(VCopObj,x))
    bs_res <- pblapply(s, function(x) VC_deco_for_bs(VCfit,x,par))
  }else{
    # Set up clusters on every core
    cl <- makeCluster(cores)
  
    # Run all libraries on Cluster
    clusterEvalQ(cl, {
      library("Hmisc")
      library("rvinecopulib")
      library("sampleSelection")
      library("Formula")
      library("survey")
    })
    
    # Names of functions to be exported
    foos <- as.vector(lsf.str(,as.environment(".GlobalEnv")))
    
    # Export results and user written function to clusters
    clusterExport(cl=cl, 
                  varlist=c("VCfit","par",foos),
                  envir=environment()
    )

    # Run parLapply on 
    #res <- parLapply(cl, s, function(x) VC_deco_for_bs(VCopObj,x))
    bs_res <- pblapply(s, function(x) VC_deco_for_bs(VCfit,x,par), cl=cl)
    
    # Stop cluster
    stopCluster(cl)
  }  
  
  ## Estimate bs standard deviation of decomposition results
  se <- list()
  if(length(bs_res[[1]])==1){
    bs_res <- list(aggregate=lapply(bs_res, function(x) x[[1]]))
  }else{
    bs_res <- list(aggregate=lapply(bs_res, function(x) x[[1]]),
                     composition_effect=lapply(bs_res, function(x) x[[2]]),
                     wage_structure_effect=lapply(bs_res, function(x) x[[3]]))
  }
  
  for(i in 1:length(bs_res)){
  res <- bs_res[[i]]
  
  # Quantiles
  quants_deco <- do.call("rbind", lapply(res, function(x) x[[1]]))
  #quants_deco <- lapply(split( quants_deco, quants_deco[,1]), function(x) data.frame(x[1,1:4],t(sqrt(diag(cov(x[,5:ncol(x)]))))))
  quants_deco_se <- lapply(split(quants_deco, quants_deco[,1]), function(x) data.frame(x[1,1:4],t(robust_bs_se(x[,5:ncol(x)]))))
  quants_deco_se <- do.call("rbind", quants_deco_se)
  names(quants_deco_se) <- names(res[[1]][[1]])

  # Other stats 
  stats_deco <- do.call("rbind", lapply(res, function(x) x[[2]]))
  #stats_deco <- lapply(split(stats_deco, stats_deco[,1]), function(x) data.frame(x[1,1:4],t(sqrt(diag(cov(x[,5:ncol(x)]))))))
  stats_deco_se <- lapply(split(stats_deco, stats_deco[,1]), function(x) data.frame(x[1,1:4],t(robust_bs_se(x[,5:ncol(x)]))))
  stats_deco_se <- do.call("rbind", stats_deco_se)
  names(stats_deco_se) <- names(res[[1]][[2]])
  
  ## Estimate bs standard deviation (i.e. square of variance) of check of fit
  # Quantiles
  quants_sim_vs_obs <- do.call("rbind", lapply(res, function(x) x[[5]]))
  #quants_sim_vs_obs <- lapply(split( quants_sim_vs_obs, quants_sim_vs_obs[,1]), function(x) data.frame(x[1,1],t(sqrt(diag(cov(x[,2:ncol(x)]))))))
  quants_sim_vs_obs_se <- lapply(split( quants_sim_vs_obs, quants_sim_vs_obs[,1]), function(x) data.frame(x[1,1],t(robust_bs_se(x[,2:ncol(x)]))))
  quants_sim_vs_obs_se <- do.call("rbind", quants_sim_vs_obs_se)
  names(quants_sim_vs_obs_se) <- names(res[[1]][[5]])
  
  # Other stats 
  stats_sim_vs_obs <- do.call("rbind", lapply(res, function(x) x[[6]]))
  #stats_sim_vs_obs <- lapply(split(stats_sim_vs_obs, stats_sim_vs_obs[,1]), function(x) data.frame(x[1,1],t(sqrt(diag(cov(x[,2:ncol(x)]))))))
  stats_sim_vs_obs_se <- lapply(split(stats_sim_vs_obs, stats_sim_vs_obs[,1]), function(x) data.frame(x[1,1],t(robust_bs_se(x[,2:ncol(x)]))))
  stats_sim_vs_obs_se <- do.call("rbind", stats_sim_vs_obs_se)
  names(stats_sim_vs_obs_se) <- names(res[[1]][[6]])
  
  ## Estimate bs quantiles of KS limiting distribution for uniform confidence bands
  # Quantile treatment effects
  quants_deco_est <- deco_res[[i]]$quants_deco
  nqs <-  ncol(quants_deco_est)-4
  quants_deco_critval <- cbind(quants_deco,quants_deco_est[rep(1:nrow(quants_deco_est),it),5:ncol(quants_deco_est)])
  quants_deco_critval <- lapply(split(quants_deco_critval,quants_deco_critval[,1]),
                                function(x) critval_uniform_ci(x[,(4+1):(4+nqs)],x[1,(4+nqs+1):(4+2*nqs)],alpha))
  quants_deco_critval <- cbind(quants_deco_se[,1:4],as.data.frame(do.call("rbind",  quants_deco_critval )))
  names(quants_deco_critval)[5] <- "t"
  
  # Diference between simulated and observed quantiles
  quants_sim_vs_obs_est_1 <- t(deco_res[[i]]$quants_sim_vs_obs)[2,]
  quants_sim_vs_obs_est_0 <- t(deco_res[[i]]$quants_sim_vs_obs)[3,]
  quants_sim_vs_obs_est_dd <- t(deco_res[[i]]$quants_sim_vs_obs)[4,]
  quants_sim_vs_obs_bs_1 <- do.call("rbind", lapply(res, function(x) t(x[[5]])[2,]))
  quants_sim_vs_obs_bs_0 <- do.call("rbind", lapply(res, function(x) t(x[[5]])[3,]))
  quants_sim_vs_obs_bs_dd <- do.call("rbind", lapply(res, function(x) t(x[[5]])[4,]))
  
  
  quants_sim_vs_obs_critval <- data.frame(group=c("group1","group0","DiffInDiff"),
                                          t=c(critval_uniform_ci(quants_sim_vs_obs_bs_1,quants_sim_vs_obs_est_1,alpha),
                                              critval_uniform_ci(quants_sim_vs_obs_bs_0,quants_sim_vs_obs_est_0,alpha),
                                              critval_uniform_ci(quants_sim_vs_obs_bs_dd,quants_sim_vs_obs_est_dd,alpha)))
  
  ### P-values of Kolmogorov-Smirnov and Cramèr-von-Mises-tests
  # Diference between simulated and observed quantiles
  quants_sim_vs_obs_KS_CMS_test <- c(KS_CMS_test(quants_sim_vs_obs_bs_1,quants_sim_vs_obs_est_1,quants_sim_vs_obs_bs_1),
                                     KS_CMS_test(quants_sim_vs_obs_bs_0,quants_sim_vs_obs_est_0,quants_sim_vs_obs_bs_0))
  names(quants_sim_vs_obs_KS_CMS_test) <- c("KS group 1","CMS group 1","KS group 0", "CMS group 0")
  
  # Return results
  se[[i]] <-  list(quants_deco_se=quants_deco_se,
              stats_deco_se=stats_deco_se,
              quants_sim_vs_obs_se=quants_sim_vs_obs_se,
              stats_sim_vs_obs_se=stats_sim_vs_obs_se,
              quants_deco_critval=quants_deco_critval,
              quants_sim_vs_obs_critval=quants_sim_vs_obs_critval,
              quants_sim_vs_obs_KS_CMS_test=quants_sim_vs_obs_KS_CMS_test
              )
  names(se)[i] <- names(deco_res)[i]
  }
  return(se)
}

##########################################################################
# Resample function
VC_resample <- function(VCfit,it,psu,stratum){
  
  n0 <- nrow(VCfit$data0)
  n1 <- nrow(VCfit$data1)

  w0 <- VCfit$data0[,"(weights)"]
  w1 <- VCfit$data1[,"(weights)"]
  
  if(is.null(psu)&is.null(stratum)){
  #Naive n-out-of-m-bootstrap over entire sample (with constant adjustment of all weights) 
  s <- lapply(1:it, function(x){
                                s0 <- sample(1:n0,n0,replace=TRUE)
                                s1 <- sample(1:n1,n1,replace=TRUE)
                                list(s0=s0,
                                     s1=s1,
                                     w0=sum(w0, na.rm=TRUE)/sum(w0[s0], na.rm=TRUE),
                                     w1=sum(w1, na.rm=TRUE)/sum(w1[s1], na.rm=TRUE))}
              )
  
  }else{
  if(is.null(stratum)){
    stratum <- list(1,1)
  }else{
    stratum <- list(VCfit$data0[,stratum],VCfit$data1[,stratum])
  }
  if(is.null(psu)){
    psu <- list(1:length(w0),1:length(w1))
    }else{
    psu <- list(VCfit$data0[,psu],VCfit$data1[,psu])
  }
  # N-out-of-m-bootstrap within single strata, computing a constant adjustment factor 
  # for sample weights within each stratum. 
  s <- lapply(1:it, function(x) resample_within_strata(stratum[[1]],
                                                       stratum[[2]],
                                                       psu[[1]],
                                                       psu[[2]],
                                                       w0,
                                                       w1))
  } 
  
  return(s)
}
  

##########################################################################
# Function for weight adjusted resampling within strata 

resample_within_strata <-  function(stratum0,stratum1,psu0,psu1,w0,w1){

# Resample PSU within every stratum
s0 <- lapply(split(psu0, stratum0), function(x) sample(unique(x),length(unique(x)),replace=TRUE))
s1 <- lapply(split(psu1, stratum1), function(x) sample(unique(x),length(unique(x)),replace=TRUE))

# Search position of single observations from a all sampled PSU
s0 <- lapply(s0,function(x) unlist(lapply(x, function(y) which(psu0 %in% y))))
s1 <- lapply(s1,function(x) unlist(lapply(x, function(y) which(psu1 %in% y))))

## Adjust weights 
# Divide sum of resampled weights by stratum by sum of observed weights by stratum
wadj0 <- do.call("c",lapply(s0, function(x) sum(w0[x])))/
         do.call("c",lapply(split(w0,stratum0), sum))
n0 <- unlist(lapply(s0, length))
wadj0 <- unlist(sapply(1:length(wadj0), function(x) rep(wadj0[x],n0[x])))

# Divide sum of resampled weights by stratum by sum of observed weights by stratum
wadj1 <- do.call("c",lapply(s1, function(x) sum(w1[x])))/
         do.call("c",lapply(split(w1,stratum1), sum))
n1 <- unlist(lapply(s1, length))
wadj1 <- unlist(sapply(1:length(wadj1), function(x) rep(wadj1[x],n1[x])))

# Return results
return(list(s0=unlist(s0),
            s1=unlist(s1),
            wadj0=wadj0,
            wadj1=wadj1))

}

##########################################################################
# Entire decomposition procedure in one function for bootstrap 

VC_deco_for_bs <- function(VCfit,s,par){
  
  # Decomposition with control variables? 
  counterw <- is.null(VCfit$VCfitC)==FALSE
  
  # Draw bootstrap sample
  data0 <- VCfit$data0[s[[1]],]  
  data1 <- VCfit$data1[s[[2]],] 
  
  # Get dimension, vartypes and structure
  var_types <- VCfit$VCfit0$var_types
  structure0 <- VCfit$VCfit0$structure
  structure1 <- VCfit$VCfit1$structure
  if(counterw){ 
  structureC <- VCfit$VCfitC$structure
  } 
  controls <- VCfit$VCfit0$controls
  formula <- VCfit$formula
  frw <- VCfit$frw
  d <- length(all.vars(formula))-1
  impute <- VCfit$impute
  
  # Re-Impute missing values if required
  if(is.null(impute)==FALSE){
    for(i in 1:length(impute)){
      data0$`(iw)` <- data0[,impute[[i]][[3]]]*s[[3]]
      data0 <- impute_wages(impute[[i]][[1]],
                            impute[[i]][[2]],
                            data=data0, weights = "(iw)"  , group=NULL, method="ml",
                            print_info=FALSE, short_info=FALSE)
      data1$`(iw)` <- data1[,impute[[i]][[3]]]*s[[4]]
      data1 <- impute_wages(impute[[i]][[1]],
                            impute[[i]][[2]],
                            data=data1, weights = "(iw)", group=NULL, method="ml",
                            print_info=FALSE, short_info=FALSE)
    }
  }
  
  # Perform reweighting if control variables used
  if(counterw){
    data0$`(group)` <- "group0"
    data1$`(group)` <- "group1"
    dataC <- rbind(data0,data1)
    dataC$`(group)` <- as.factor(dataC$`(group)`)
    invisible(capture.output(rw <- dfl_deco(frw,dataC,
                             weights=`(weights)`, group=`(group)`, 
                             reference="group1", stats=FALSE)$psi))
    sel <- which(dataC[,"(group)"]=="group1")
    dataC <- dataC[sel,]
    dataC[,"(weights)"] <- dataC[,"(weights)"] *rw[sel,1] 
  }
    
  # Compute pseudo-copula obs for bootstrapped sample (i.e. F(x))
  u0 <- apply(data0[,2:(d+1)], 2, function(x) rank.wt(x,data0[,"(weights)"], ties_method = "max"))
  u1 <- apply(data1[,2:(d+1)], 2, function(x) rank.wt(x,data1[,"(weights)"], ties_method = "max"))
  if(counterw){ 
    uC <- apply(dataC[,2:(d+1)], 2, function(x) rank.wt(x,dataC[,"(weights)"], ties_method = "max"))
  }
  
  if(any(var_types=="d")){
    # Compute left-sided limit of discrete/integer variables (i.e. F(x^-)=F(x-1))
    sel <- which(var_types=="d")+1
    u0minus <- apply(as.data.frame(data0[,sel]), 2, function(x) rank.wt(x,data0[,"(weights)"], ties_method = "max", left_limit_CDF = TRUE))
    u1minus <- apply(as.data.frame(data1[,sel]), 2, function(x) rank.wt(x,data1[,"(weights)"], ties_method = "max", left_limit_CDF = TRUE))
    colnames(u0minus) <- colnames(u1minus) <- paste0(colnames(u0minus),"_minus")
    u0 <- cbind(u0,u0minus)
    u1 <- cbind(u1,u1minus)
    if(counterw){ 
      uCminus <- apply(as.data.frame(dataC[,sel]), 2, function(x) rank.wt(x,dataC[,"(weights)"], ties_method = "max", left_limit_CDF = TRUE))
      colnames(uCminus) <-  colnames(u1minus)
      uC <- cbind(uC,uCminus)  
    }
  }
  
  ## Re-estimate copula with adjusted weights
  VCfit0 <- vinecop(u0, var_types=var_types, 
                    family_set = controls$family_set, 
                    par_method=controls$par_method,
                    nonpar_method = controls$nonpar_method,
                    structure=structure0, mult=controls$mult, 
                    weights=data0[,"(weights)"]*s[[3]], tree_crit=controls$tree_crit,
                    cores=1, show_trace=FALSE)
  
  VCfit1 <- vinecop(u1, var_types=var_types,
                    family_set = controls$family_set, 
                    par_method=controls$par_method,
                    nonpar_method = controls$nonpar_method,
                    structure=structure1, mult=controls$mult, 
                    weights=data1[,"(weights)"]*s[[4]], tree_crit=controls$tree_crit,
                    cores=1, show_trace=FALSE)
  
  if(counterw){ 
  VCfitC <- vinecop(uC, var_types=var_types, 
                      family_set = controls$family_set, 
                      par_method=controls$par_method,
                      nonpar_method = controls$nonpar_method,
                      structure=structureC, mult=controls$mult, 
                      weights=dataC[,"(weights)"], tree_crit=controls$tree_crit,
                      cores=1, show_trace=FALSE)
  }

  ## Re-perform decomposition
  deco_res <- VC_deco_sim(VCfit0, VCfit1, formula=formula,# jit=jit, jit_bd=jit_bd,
                        data0=data0, data1=data1, 
                        reference=par$reference,  n=par$n, tau=par$tau,
                        multiply = par$multiply, quasi = par$quasi, 
                        detailed_deco = par$detailed_deco, 
                        export_usim=FALSE, print_info = FALSE, cores=1)
  if(counterw){ 
    # "Composition" part
    import_usim <-  deco_res$usim_all["1-1-1-1-1-1"]
    deco_res_X <- VC_deco_sim(VCfitC, VCfit1,
                              formula=formula,# jit=jit, jit_bd=jit_bd,
                              data0=dataC, data1=data1,
                              reference=par$reference,  n=par$n, tau=par$tau,
                              multiply = par$multiply, quasi = par$quasi, 
                              detailed_deco = par$detailed_deco, 
                              export_usim=FALSE, print_info = FALSE, cores=1,
                              import_usim=import_usim)
    # "Wage Structure" part
    import_usim <-  deco_res$usim_all["0-0-0-0-0-0"]
    deco_res_WS <- VC_deco_sim(VCfit0, VCfitC,
                               formula=formula,# jit=jit, jit_bd=jit_bd,
                               data0=data0, data1=dataC,
                               reference=par$reference,  n=par$n, tau=par$tau,
                               multiply = par$multiply, quasi = par$quasi, 
                               detailed_deco = par$detailed_deco, 
                               export_usim=FALSE, print_info = FALSE, cores=1,
                               import_usim=import_usim)
    
    deco_res <- list(aggregate=deco_res,
                     composition_effect=deco_res_X,
                     wage_structure_effect=deco_res_WS)
  } else{
    deco_res <- list(aggregate=deco_res)
  }

  
  return(deco_res)
  
}

##########################################################################
# Compute robust bs s.e.
robust_bs_se <- function(x){
  if(is.null(ncol(x))){x <- as.matrix(x)}
  apply(x,2,function(x) (quantile(x,0.75, na.rm=TRUE) - quantile(x,0.25,na.rm=TRUE))/(qnorm(0.75)-qnorm(0.25)))
}

##########################################################################
# Compute bs cirtical value for uniform confidence interval 
# by finding quantile of bootstrapped Kolmogorov-Smirnov test distribution
# (slightly adapted code from Chen et al. (2017), Counterfactual::VarainceEval()

critval_uniform_ci <- function(quants_deco_bs,quants_deco,alpha=0.05){
  Vuqf <- robust_bs_se(quants_deco_bs)^2
  nqs <- length(quants_deco)
  reps <- nrow(quants_deco_bs)
  Kuqf <- sqrt((quants_deco_bs - 
                kronecker(matrix(as.numeric(quants_deco),1,nqs),
                          matrix(1, reps, 1)))^2
               /kronecker(matrix(Vuqf, 1, nqs),
                          matrix(1, reps, 1)))
  Kuqfsel <- Kuqf[, apply(Kuqf, 2, function(x) all(is.finite(x)))]
  Kmaxuqf <- apply(Kuqfsel, 1, max)
  Kalpha <- quantile(Kmaxuqf, 1 - alpha, na.rm = TRUE, names = FALSE)
  return(Kalpha)
}

##########################################################################
# Compute p-values of Kolmogorov-Smirnov- and Cram?r-von-Mises-type tests 

KS_CMS_test <- function(boot_test_numerator, 
                        obs_test_numerator,
                        variable_for_variance){
  Vuqf <- robust_bs_se(variable_for_variance)^2
  Vuqf = Vuqf + 1e-09
  reps <- nrow(boot_test_numerator)
  nqs = ncol(boot_test_numerator)
  Kuqf = sqrt(boot_test_numerator^2/
                kronecker(matrix(Vuqf, 1, nqs), matrix(1, reps, 1)))
  Kmaxuqf = apply(Kuqf, 1, max)
  KSstat = max(sqrt(obs_test_numerator^2/Vuqf))
  Kuqf2 = Kuqf^2
  Kmeanuqf = rowMeans(Kuqf2)
  CMSstat = mean(obs_test_numerator^2/Vuqf)
  testboot = c(mean(Kmaxuqf > KSstat), mean(Kmeanuqf > CMSstat))
  return(testboot)
}

#library(Counterfactual)
#counterfactual()
#Counterfactual::VarianceEval()
#Counterfactual::InferenceTestingEval()
#Counterfactual::TestingEval()


stats_marginals_est_for_bs <- function(marginal_stats,s){
  
  marginal_stats$data0 <- marginal_stats$data0[s[[1]],]
  marginal_stats$data1 <- marginal_stats$data1[s[[2]],]
  marginal_stats$data0[,"(weights)"] <- marginal_stats$data0[,"(weights)"]*s[[3]]
  marginal_stats$data1[,"(weights)"] <- marginal_stats$data1[,"(weights)"]*s[[4]]
  
  res <- stats_marginals_est(marginal_stats$data0,
                             marginal_stats$data1,
                             tau=marginal_stats$tau,
                             kendall=marginal_stats$kendall,
                             trimgini=marginal_stats$trimgini,
                             trimcor=marginal_stats$trimcor
                             )
  
  return(res)
}

##########################################################################
# Function to bootstrap marginal statistics

stats_marginals_bs <- function(marginal_stats, it=100, 
                         ncore=1, psu=NULL, stratum=NULL, alpha=0.05){
  
  # Resample 
  s <- VC_resample(marginal_stats, it=it, psu=psu, stratum=stratum)
  
  # Restimate (parallelized if requested)
  if(ncore==1){  
    #res <- lapply(s, function(x) VC_deco_for_bs(VCopObj,x))
    res <- pblapply(s, function(x) stats_marginals_est_for_bs(marginal_stats,x))
    
  }else{
    # Set up clusters on every core
    cl <- makeCluster(ncore)
    
    # Run all libraries on Cluster
    clusterEvalQ(cl, {
      library("Hmisc")
      library("VineCopula")   
    })
    
    # Names of functions to be exported
    foos <- as.vector(lsf.str(,as.environment(".GlobalEnv")))
    
    # Export results and user written function to clusters
    clusterExport(cl=cl, 
                  varlist=c("marginal_stats",foos),
                  envir=environment()
    )
    
    # Run parLapply on 
    #res <- parLapply(cl, s, function(x) VC_deco_for_bs(VCopObj,x))
    res <- pblapply(s, function(x) stats_marginals_est_for_bs(marginal_stats,x), cl=cl)
    
    # Stop cluster
    stopCluster(cl)
  }  
  
  ## Estimate bs standard deviation of decomposition results
  # Quantiles
  quants <- do.call("rbind", lapply(res, function(x) x[[1]]))
  quants_se <- lapply(split(quants , quants[,c("tau","group")]), function(x) data.frame(x[1,1:2],t(robust_bs_se(x[,3:7]))))
  quants_se <- do.call("rbind", quants_se)
  
  # Other stats 
  stats <- do.call("rbind", lapply(res, function(x) x[[2]]))
  stats_se <- lapply(split(stats, stats[,c("stat","group")]), function(x) data.frame(x[1,1:2],t(robust_bs_se(x[,3:7]))))
  stats_se <- do.call("rbind", stats_se)

  # Corr measures
  corrs <- do.call("rbind", lapply(res, function(x) x[[3]]))
  corrs_se <- lapply(split(corrs, corrs[,c("type","vars","group")]), function(x) data.frame(x[1,c(1,2,4)],t(robust_bs_se(x[,3]))))
  corrs_se <- do.call("rbind", corrs_se)
  names(corrs_se)[4] <- "rho"
  
  # Return results
  return(list(quants_se=quants_se,
              stats_se=stats_se,
              corrs_se=corrs_se))
}

