##########################################################################
#### 1 - VCdeco - estimation using rvinecopulib
##########################################################################

##########################################################################
#### Load library
# require(rvinecopulib)

##########################################################################
#### New wrapper function
VC_deco <- function(formula,data,group=numeric(),weights=numeric(),
                    ## Select best initial trees
                    find_best_tree=FALSE,  #
                    ## Estimation
                    #jit=NA,               # which variables should be jittered?
                    #jit_bd=0.5,           # jittering bandwidth
                    var_types = "c",       #
                    family_set="tll",      # specify parametric or nonparaemtric copula modell 
                    structure=NA,          # copula matrix; if NA, algorithm searches for best tree
                    par_method = "mle",    # "mle" or "itau"  (MM) est. for parametric estimation
                    nonpar_method="quadratic", #est. method nonparametric models: "quadratic"= std. trans estimator with local likelihood approx. of order two
                    mult=1,                # bandwidth multiplier
                    print_info=TRUE,       # print progress 
                    group0=NULL,           # define group 0
                    tree_crit="rho",       # Criterion for tree selection
                    ## Decomposition
                    deco=TRUE,             # perform decomposition (if FALSE only copulas are fitted)
                    detailed_deco=TRUE,    # perform detailed copula decomposition (i.e. isolate influence of single pair copulas)
                    reference=1,           # copula of reference group (c(0,1)) is used to analyze the impact of marginals
                    n=1000,                # number simulated observations
                    tau=seq(5,95,5)/100,   # quantiles to estimated
                    multiply=(1/40),       # multiply simulated incomes by factor 
                    quasi=TRUE,            # quasi-random numbers or pseudo-random numbers?
                    export_usim=TRUE,      # export simulated quasi copula data
                    ## Imputation   
                    impute=NULL,           # list of models with variables to impute
                    impute_bs_only=TRUE,   # (re-)impute only when bootstrapping?
                    ## Bootstrap
                    bootstrap=FALSE,       # boostrap 
                    it=100,                # boostrap iterations
                    psu=NULL,              # psu 
                    stratum=NULL,          # stratum   
                    alpha=0.05,            # significance for uniform ci
                    cores=1                # number of cores
){
  
  #---------------------------------------------------------------
  # 1 Get data and parameters
  
  #Use match.call to create data.frame
  #mf = match.call()
  #m = match(c("formula", "data", "weights","group"), names(mf), 0)
  #mf = mf[c(1, m)]
  #mf$drop.unused.levels = TRUE
  #mf[[1]] = as.name("model.frame")
  
  # Check if control variables included
  f <- as.Formula(formula)
  formula_select <- formula(f, collapse=TRUE)
  counterw <- length(f)[2]==2
  if(counterw){
    frw <- formula(f,rhs=2)
  }else{
    frw <- NULL
  }
  f <- formula(f,rhs=1)
 
  #Add psu and stratum to formula in order to integrate it to model frame
  if(is.null(psu)==TRUE&is.null(stratum)==FALSE){
    #mf$formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)], collapse=" ")," +",stratum)) 
    formula_select <- as.formula(paste0(paste0(as.character(formula_select)[c(2,1,3)], collapse=" ")," +",stratum))
  }else if(is.null(psu)==FALSE&is.null(stratum)==TRUE){
    #mf$formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)], collapse=" ")," +",psu))   
    formula_select <- as.formula(paste0(paste0(as.character(formula_select)[c(2,1,3)], collapse=" ")," +",psu))   
  }else if(is.null(psu)==FALSE&is.null(stratum)==FALSE){
    #mf$formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)], collapse=" ")," +",psu," +", stratum))   
    formula_select <- as.formula(paste0(paste0(as.character(formula_select)[c(2,1,3)], collapse=" ")," +",psu," +", stratum))   
  }
  
  # Add independent variables used for imputations to "selection" formula
  if(is.null(impute)==FALSE){
    l <- length(impute)
    selw <- seq(3,l*3,3)
    selm <- setdiff(1:(3*l),selw)
    vars_to_add <- unique(unlist(lapply(unlist(impute)[selm],function(x) as.character(x)[c(2,3)])))
    vars_to_add <- c(vars_to_add, unique(unlist(unlist(impute)[selw])))
    vars_to_add <- paste(vars_to_add, collapse=" + ")
    formula_select <- as.formula(paste0(paste(as.character(formula_select)[c(2,1,3)], collapse=" "), " + ", vars_to_add))
  }
  
  # Evaluate model.frame
  # mf = eval.parent(mf)
  # mt = attr(mf, "terms")
  
  #Retrieve weights
  #if(is.null(model.weights(mf))){
  #  mf$`(weights)` <- rep(1,nrow(mf))
  #}
  
  # Check group variable
  if(length(group)==0|is.numeric(group)&length(group)!=nrow(data)){
    stop("group must be speficied or have length of data")
  }else{
    data$group <- data[,group]
  }
  
  #Specify weigths if empty
  if(length(weights)==0){
    data$weights <- rep(1,nrow(data))
  }else{
    data$weights <- data[,weights]
  }

  #---------------------------------------------------------------
  # 2 Impute missing values if required
  if(is.null(impute)==FALSE&impute_bs_only==FALSE){
    for(i in 1:length(impute)){
      if(print_info){cat("\n Imputing", as.character(impute[[i]][[2]])[2], "\n")}
      data <- impute_wages(impute[[i]][[1]],
                           impute[[i]][[2]],
                           data=data, weights = impute[[i]][[3]], group=group, method="ml",
                           print_info=print_info)
    }
  }
  
  #formula <- EARNINGS ~ WAGE_m:ANNUALHOURS_m + WAGE_f:ANNUALHOURS_f
  #formula_select <- EARNINGS ~ WAGE_m:ANNUALHOURS_m + WAGE_f:ANNUALHOURS_f
  #tail(get_all_vars(formula_select,cdata))
  #data <- cdata
  #mf <- get_all_vars(formula_select,cdata,weights="ASECWTH_f", group="YEAR_m")
  #head(get_all_vars(formula_select,cdata,weights="ASECWTH_f", group="YEAR_m"))
  
  #---------------------------------------------------------------
  # 3 Split sample and create copula data 
  # Keep only required variables
  #mf <- model.frame(formula_select,data,weights=weights, group=data$`(group)`)
  mf <- get_all_vars(formula_select,data,weights=weights,group=group)

  # Rename weights and group variable
  names(mf)[is.element(names(mf),c("weights","group"))] <- c("(weights)","(group)")
  
  # Select group 0
  if(is.null(group0)){
    if(is.factor(mf[,"(group)"])){
      sel0 <- levels(mf[,"(group)"])[1]
      sel1 <- levels(mf[,"(group)"])[2]
    }else{
      sel0 <- sort(unique(mf[,"(group)"]))[1]
      sel1 <- sort(unique(mf[,"(group)"]))[2]
    }
  }else{
    if(length(which(mf[,"(group)"]==group0))==0){
    stop("No valid group defined!")
    }else{
    sel0 <- group0
    sel1 <- base::setdiff(unique(unique(mf[,"(group)"])),sel0)
    }
  }
  
  # Adjust formula
  formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)],collapse=" ")," - 1"))
  
  #drop observations with zero weight
  if(print_info==TRUE){
    cat(paste0(c("\nDropped because of 0 household weight:", length(which(mf[,"(weights)"]==0&mf[,"(group)"]==sel0)),"in group 0 and",
                 length(which(mf[,"(weights)"]==0&mf[,"(group)"]!=sel0)) ,"in group 1.\n")))
  }
  mf <- mf[which(mf[,"(weights)"]!=0),]
  
  # Select 
  data0 <- mf[which(mf[,"(group)"]==sel0),]
  data1 <- mf[which(mf[,"(group)"]==sel1),] 
  
  #---------------------------------------------------------------
  # 4 Perform reweighting if control variables used
  if(counterw){
    if(print_info){cat("Estimating weights for counterfactuals... \n")}
    invisible(capture.output(rw <- dfl_deco(frw,mf,
                             weights=`(weights)`, group=`(group)`, 
                             reference=sel1, stats=FALSE)$psi))
    s <- which(mf[,"(group)"]==sel1)
    dataC <- mf[s,]
    dataC[,"(weights)"] <- dataC[,"(weights)"]*rw[s] 
  }
  
  #---------------------------------------------------------------
  # 5 Create copula data 
  # get dimension of copula and varnames variables
  d <- max(length(all.vars(f))-1,2) #ncol(data0)-3
  l <- length(all.vars(f))
  varnames <- all.vars(f)[(l-d+1):l] 
  
  # Add jittering to discrete data
  #data0[,jit+1] <- apply(data0[,jit+1], 2, function(x) add_jitter(x))
  #data1[,jit+1] <- apply(data1[,jit+1], 2, function(x) add_jitter(x))
  #if(is.null(jit)==FALSE){
  #  if(length(jit_bd)!=d){
  #    jit_bd_o <- jit_bd
  #    jit_bd <- rep(0,d)
  #    jit_bd[jit] <- jit_bd_o
  #  }
  #  data0[,jit+1] <- do.call("cbind",lapply(jit,function(x) add_jitter(data0[,x+1],bd=jit_bd[x])))
  #  data1[,jit+1] <- do.call("cbind",lapply(jit,function(x) add_jitter(data1[,x+1],bd=jit_bd[x])))
  #}
  #Define var_types
  if(length(var_types)==1){
    var_types <- rep(var_types,d)
  }
  
  # Compute pseudo-copula obs (i.e. F(x))
  s <- ifelse(l==d,0,1)
  u0 <- apply(data0[,(1:d)+s], 2, function(x) rank.wt(x,data0[,"(weights)"], ties_method = "max"))
  u1 <- apply(data1[,(1:d)+s], 2, function(x) rank.wt(x,data1[,"(weights)"], ties_method = "max"))
  if(counterw){ 
  uC <- apply(dataC[,(1:d)+s], 2, function(x) rank.wt(x,dataC[,"(weights)"], ties_method = "max"))
  }
  
  if(any(var_types=="d")){
  # Compute left-sided limit of discrete/integer variables, i.e. F(x^-)=F(x-1)
  sel <- which(var_types=="d")+s
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
  
  #6 --------------------------------------------------------------
  # find best tree if required 
  if(find_best_tree){
  best_trees <- VC_find_best_tree(u0, u1, 
                               w0=data0[,"(weights)"], w1=data1[,"(weights)"],
                               nbest=12, var_types=var_types, tree_crit=tree_crit,
                               family_set=family_set)
  return(best_trees)
  }
  
  #2 ---------------------------------------------------------------
  # Fit model 
  show_trace <- print_info
  estcores <- ifelse(bootstrap,1,cores)
  if(class(structure)[1]=="list"&class(structure)[1]!="rvine_structure"){
    structure0 <- structure[[1]]
    structure1 <- structure[[2]]
    structureC <- structure[[2]]
  }else{
    structure0 <- structure
    structure1 <- structure
    structureC <- structure
  }

  if(print_info){cat("Density estimation group 0... \n")}
  VCfit0 <- vinecop(u0, var_types=var_types, family_set = family_set, 
                    nonpar_method = nonpar_method,
                    par_method = par_method,
                    structure=structure0, mult=mult, 
                    weights=data0[,"(weights)"], tree_crit=tree_crit,
                    cores=estcores, show_trace=show_trace)

  if(print_info){cat("\nDensity estimation group 1... \n")}
  VCfit1 <- rvinecopulib::vinecop(u1, var_types=var_types, family_set = family_set, 
                    nonpar_method = nonpar_method,
                    par_method = par_method,
                    structure=structure1, mult=mult, 
                    weights=data1[,"(weights)"], tree_crit=tree_crit,
                    cores=estcores, show_trace=show_trace)


  if(counterw){ 
  if(print_info){cat("\nDensity estimation reweighted group 1... \n")}
  VCfitC <- rvinecopulib::vinecop(uC, var_types=var_types, family_set = family_set, 
                      nonpar_method = nonpar_method,
                      structure=structureC, mult=mult, 
                      weights=dataC[,"(weights)"], tree_crit=tree_crit,
                      cores=estcores, show_trace=show_trace)
  }else{
    VCfitC=NULL
    dataC=NULL
    uC=NULL
  }
  
  VCfit <- list(VCfit0=VCfit0, VCfit1=VCfit1, VCfitC=VCfitC,
                formula=f, frw=frw, #jit=jit, 
                data0=data0, data1=data1,dataC=dataC,
                u0=u0, u1=u1, uC=uC,
                impute=impute)
  
  
  #7 ---------------------------------------------------------------
  # Perform decomposition
  if(deco){
    if(print_info){cat("\nAggregate decomposition... \n")}
    deco_res <- VC_deco_sim(VCfit0, VCfit1,
                            formula=f,# jit=jit, jit_bd=jit_bd,
                            data0=data0, data1=data1,
                            reference=reference,  n=n, tau=tau,
                            multiply = multiply, quasi = quasi, 
                            detailed_deco = detailed_deco, print_info = print_info,
                            cores=estcores)
    if(counterw){ 
    if(print_info){cat("\nComputing composition effect... \n")}
    # "Composition" part
    import_usim <-  deco_res$usim_all["1-1-1-1-1-1"]
    deco_res_X <- VC_deco_sim(VCfitC, VCfit1,
                              formula=f,# jit=jit, jit_bd=jit_bd,
                              data0=dataC, data1=data1,
                              reference=reference,  n=n, tau=tau,
                              multiply = multiply, quasi = quasi, 
                              detailed_deco = detailed_deco, print_info = print_info,
                              import_usim=import_usim, 
                              cores=estcores)
    # "Wage Structure" part
    if(print_info){cat("\nComputing wage structure effect... \n")}
    import_usim <-  deco_res$usim_all["0-0-0-0-0-0"]
    deco_res_WS <- VC_deco_sim(VCfit0, VCfitC,
                               formula=f,# jit=jit, jit_bd=jit_bd,
                               data0=data0, data1=dataC,
                               reference=reference,  n=n, tau=tau,
                               multiply = multiply, quasi = quasi, 
                               detailed_deco = detailed_deco, print_info = print_info,
                               import_usim=import_usim, 
                               cores=estcores)
    deco_res <- list(aggregate=deco_res,
                     composition_effect=deco_res_X,
                     wage_structure_effect=deco_res_WS)
    }else{
    deco_res <- list(aggregate=deco_res)
    }
    
    
  }else{
    deco_res <- NULL  
  }
  
  #8 ---------------------------------------------------------------
  # Boostrap
  if(bootstrap){
    cat("\nBootstrapping...\n")
    VCdeco_se <- VC_deco_bs(VCfit, deco_res, it=it, cores=cores,
                            psu=psu, stratum=stratum, alpha=alpha)  
  }else{
    VCdeco_se <- NULL  
  }
  
  #9 ---------------------------------------------------------------
  # Return results
  return(list(VCfit=VCfit,
              VCdeco=deco_res,
              VCdeco_se=VCdeco_se))
}

###########################################################################
#### Find "best" first trees
VC_find_best_tree <- function(u0, u1, w0, w1, nbest=5,
                         var_types, tree_crit, family_set){
  
  # Get varnames and dimension
  d <- length(var_types)#ncol(u0)
  
  # Estimate weights
  cat("Estimate weights (", tree_crit, ")...\n")
  if(tree_crit=="tau"){
    Ktau0 <- VineCopula::TauMatrix(u0[,1:d], weights=w0)
    Ktau1 <- VineCopula::TauMatrix(u1[,1:d], weights=w1)
  }else if(tree_crit=="rho"){
    Ktau0 <- cov.wt(u0[,1:d],w=w0, cor=TRUE)$cor
    Ktau1 <- cov.wt(u1[,1:d],w=w1, cor=TRUE)$cor  
  }else{
    Ktau0 <- bicop_AICs(u0,w=w0, var_types=var_types)
    Ktau1 <- bicop_AICs(u1,w=w1, var_types=var_types)
  }
  
  # Find best CVine combinations
  a <- rep(1:d)
  b <- lapply(a, function(x) setdiff(a,x))
  ordersC <- matrix(c(a,unlist(b)),d,d,byrow=TRUE)
  EdgesByKtau0 <- t(apply(ordersC,1,function(x) gettaus(x,Ktau0,"CVine")))
  EdgesByKtau1 <- t(apply(ordersC,1,function(x) gettaus(x,Ktau1,"CVine")))
  TauSumOverEdges0 <- rowSums(abs(EdgesByKtau0))
  TauSumOverEdges1 <- rowSums(abs(EdgesByKtau1))
  os0 <- order(TauSumOverEdges0, decreasing =TRUE)
  os1 <- order(TauSumOverEdges1, decreasing =TRUE)
  best0C <- ordersC[os0,]
  best1C <- ordersC[os1,]
  besttau0C <- EdgesByKtau0[os0,]
  besttau1C <- EdgesByKtau1[os1,] 
  
  # Find unique DVine combinations
  orders <- expand.grid(as.data.frame(matrix(rep(1:d,d),d,d)))
  orders <- orders[which(apply(orders,1,function(x) length(unique(x))==d)),]
  ordersD <- orders[1:(nrow(orders)/2),]
  EdgesByKtau0 <- t(apply(ordersD,1,function(x) gettaus(x,Ktau0,"DVine")))
  EdgesByKtau1 <- t(apply(ordersD,1,function(x) gettaus(x,Ktau1,"DVine")))
  TauSumOverEdges0 <- rowSums(abs(EdgesByKtau0))
  TauSumOverEdges1 <- rowSums(abs(EdgesByKtau1))
  os0 <- order(TauSumOverEdges0, decreasing =TRUE)
  os1 <- order(TauSumOverEdges1, decreasing =TRUE)
  best0D <- ordersD[os0,]
  best1D <- ordersD[os1,]
  besttau0D <- EdgesByKtau0[os0,]
  besttau1D <- EdgesByKtau1[os1,] 
  
  maxd <- min(nbest,nrow(best0D))
  maxc <- min(nbest,nrow(best0C))
  best0D <- best0D[c(1:maxd),]
  best1D <- best1D[c(1:maxd),]
  best0C <- best0C[c(1:maxc),]
  best1C <- best1C[c(1:maxc),]
  
  sD <- setdiff(1:maxd,match(apply(best0D ,1,function(x) paste0(x,collapse="-")), 
               apply(best1D,1,function(x) paste0(x,collapse="-"))))
  sC <- setdiff(1:maxc,match(apply(best0C ,1,function(x) paste0(x,collapse="-")), 
                             apply(best1C,1,function(x) paste0(x,collapse="-"))))
      
  bestD <- rbind(best0D,best1D[sD,])
  bestC <- rbind(best0C,best1C[sC,])  
                    
  dvine <- apply(bestD,1, function(x) rvinecopulib::dvine_structure(rev(x),trunc_lvl = 1))
  cvine <- apply(bestC,1, function(x) rvinecopulib::cvine_structure(rev(x),trunc_lvl = 1))
  
  # Fit full models
  cat("Fit models...\n")
  dvine_fits <- pbapply::pblapply(dvine, function(x) list(fit0=rvinecopulib::vinecop(u0, var_types=var_types, family_set=family_set,
                                                             weights=w0, tree_crit="rho", structure=x),
                                                          fit1=rvinecopulib::vinecop(u1, var_types=var_types, family_set=family_set,
                                                             weights=w1, tree_crit="rho", structure=x)))
  
  cvine_fits <- pbapply::pblapply(cvine, function(x) list(fit0=rvinecopulib::vinecop(u0, var_types=var_types, family_set=family_set,
                                                                       weights=w0, tree_crit="rho", structure=x),
                                                          fit1=rvinecopulib::vinecop(u1, var_types=var_types, family_set=family_set,
                                                                weights=w1, tree_crit="rho", structure=x)))
  
  
  # Get ICs and export
  dvine_IC <- do.call("rbind", lapply(dvine_fits, function(x) data.frame(AIC0=AIC(x[[1]]),
                                                          BIC0=BIC(x[[1]]),
                                                          AIC1=AIC(x[[2]]),
                                                          BIC1=BIC(x[[2]]))))
  
  vine_names <- do.call("c",lapply(dvine, function(x) paste0("D-vine ",paste0(x$order, collapse="-"), collapse= " ")))
  dvine_IC <- cbind(vine_names,dvine_IC)
  
  cvine_IC <- do.call("rbind", lapply(cvine_fits, function(x) data.frame(AIC0=AIC(x[[1]]),
                                                                         BIC0=BIC(x[[1]]),
                                                                         AIC1=AIC(x[[2]]),
                                                                         BIC1=BIC(x[[2]]))))
  vine_names <- do.call("c",lapply(cvine, function(x) paste0("C-vine ",paste(paste(x$order[length(x$order)], x$order[1:(length(x$order)-1)], sep="-"), collapse=", "), collapse= " ")))
  cvine_IC <- cbind(vine_names,cvine_IC)
  vine_IC <- rbind(dvine_IC, cvine_IC)
  
  return(list(vine_IC=vine_IC, dvine=dvine, cvine=cvine, Ktau0=Ktau0, Ktau0=Ktau1))

} 

############################################################################
#### Function to get tau of first tree DVine-CVine orders
gettaus <- function(x,tau,type="DVine"){
  res <- NULL
  if(type=="DVine"){
    for(i in 2:length(x)){
        res <-c(res,tau[x[i-1],x[i]])
    }
  }else{
    for(i in 2:length(x)){
       res <- c(res,tau[x[1],x[i]])
    }    
  }
  res
}

############################################################################
#### Best AIC 
bicop_AICs <- function(u,w,var_types){
  d<- length(var_types)
  # If we deal with continous variables, we need to add left-limit cdf
  add <- NULL
  j=1
  for(i in 1:d){
    add <- c(add,ifelse(var_types[i]=="d",d+j,NA))
    j=j+ifelse(var_types[i]=="d",1,0)
  }
  # Loop with bivariate copulas
  AICs <- diag(rep(0,d))
  for(j in 1:(d-1)){
    for(i in (j+1):d){
     sel <- c(i,j,na.exclude(add[c(i,j)]))
     AICs[i,j] <- AIC(rvinecopulib::bicop(u[,sel],var_types=var_types[c(i,j)], family_set="tll", weights=w))
    }
  }
  AICs[upper.tri(AICs)] <- AICs[lower.tri(AICs)]
  AICs
}

############################################################################
#### Weigthed partial correlation function 
wtd.pcor <- function(x,z=1,w){
  d <- ncol(x)
  res <- do.call("cbind",lapply(setdiff(1:d,z), function(y) resid(lm(x[,y] ~ x[,z], weights=w))))
  cov.wt(res,w=w, cor=TRUE)$cor
}


##########################################################################
#### Jittering function
# Jitter dist. like in Example 1 of Nagler (2018: 5)
# unifom with theta=0, alternative: theta = 0.3
add_jitter <- function(x, theta=0, nu=5, bd=0.5, quasi=TRUE) {
  n <- length(x)
  if(quasi){
    u <- qunif(qrng::ghalton(n, d = 1, method="generalized"),-bd,bd)[sample(n)]
    b <- qbeta(qrng::ghalton(n, d = 1, method="generalized"),nu,nu)[sample(n)]
  }else{
    u <- runif(n,-bd,bd)
    b <- rbeta(n,nu,nu)
  }
  x <- x+u+theta*(b-0.5)
  return(x)
}





