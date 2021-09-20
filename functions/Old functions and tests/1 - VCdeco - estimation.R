##########################################################################
#### 1 - VCdeco - estimation
##########################################################################

#### The following functions estimate a non-parametric vine copula 
#### simulataneously to the data of two groups

require("VineCopula")    # Vine copula modeling
require("kdecopula")     # Copula density estimation

##########################################################################
#### Costum jittering function
# Jitter dist. like in Example 1 of Nagler (2018: 5)
# unifom with theta=0, alternative: theta = 0.3
add_jitter <- function(x, theta=0, nu=5, bd=0.5) {
  n <- length(x)
  x <- round(x)+runif(n,-bd,bd)+theta*(rbeta(n,nu,nu)-0.5)
  return(x)
}

##########################################################################
#### Two auxiliary functions 

# Create vector
create_vector <- function(x){
                  if(length(x)!=6){
                    x <- rep(x[1],6)
                  }
                  return(x)
                  }

# Extracts coefficients from a copula matrix and stores it to vector
extract_copulamatrix <- function(x){
  if(is.matrix(x)){
    x <- x[lower.tri(x)][c(3,5,6,2,4,1)]
  }else if(length(x)==1){
    x <- rep(x,6)
  }
  return(x)
}

# Create D vine copula matrix
create_copulamatrix <- function(x){
  x <- diag(x)
  k=1
  for(i in nrow(x):2){
    x[i,1:(i-1)] <- diag(x)[-(k)]
    k <- c(k,k+1)
  }
  return(x)
}


##########################################################################
#### Fit an non-parametric vine copula
## D-Vine or C-Vine structure, four variables 

# Matrix <- matrix(c(1, 0, 0, 0,
#                    4, 2, 0, 0,
#                    3, 4, 3, 0,
#                    2, 3, 4, 4), byrow=TRUE, ncol=4)

# Matrix <- matrix(c(4, 0, 0, 0,
#                    3, 3, 0, 0,
#                    2, 2, 2, 0,
#                    1, 1, 1, 1), byrow=TRUE, ncol=4)

# Non-parametric kernel density estimation
kdeVineCop <- function(u, w=NULL,
                       order=1:4,
                       type="DVine",
                       varnames=NULL,
                       bw=NA, mult=1, method="TLL2nn", knots=30,
                       #wadj=0,
                       print_info=TRUE){
  
  n <- nrow(u)
  if(is.null(w)){w <- rep(1,n)}  
  
  #nadj <- n
  # Check for weights
  #if(is.null(w)==FALSE){
    # Adjust for weights if there a weights specified 
    
    #w <- round(w/min(w),wadj)*10^wadj
    #w <- unlist(sapply(1:length(w), function(x) rep(x,w[x])))
    
    # Multiply observation according to their sample weights
    #u <- u[w,]
    
    # Number of actual observations
    #nadj <- length(w)
    
    #if(print_info==TRUE){cat("Observations:",n,"unadjusted, ", nadj, "adjusted.\n" )}
  #} 
  
  # Set up parameter matrices
  bw <- create_vector(bw)
  mult <- create_vector(mult)#*(nadj/n)
  method <- create_vector(method)
  knots <- create_vector(knots)
  
  ##########################################################################
  ### Reorder in order to get right D-Vine structure
  
  #if(is.matrix(order)){
  #  order <- diag(order)
  #}else{
  #Matrix <- create_copulamatrix(order)
  #}
  
  # Rearrange imput matrix
  u <- u[,order]
  
  #  Retrieve names
  if(is.null(varnames)){
    varnames <- colnames(u)
  }

  ##########################################################################
  ### C- or D-vine? Selection of variables
  
  if(type=="DVine"){
  #########################################################################
  ### Fit D-Vine if required
  #
  #T1: v1   <->   v2  <->  v3  <-> v4
  #T2:  v1,v2  <->  v2,v3   <->  v3,v4
  #T3:   v1,v3;v2  <->  v2,v4;v3
  #             v1,v4;v2,v3  
    
  ##########################################################################
  ### DVine Tree 1: 
  if(print_info==TRUE){cat("Tree 1...")}

  # Edge 1.1: u1,u2
  CT1E1 <- wtd.kdecop(u[,c(1,2)], weights=w, method=method[1], mult=mult[1], knots=knots[1], bw=bw[1])
  # Edge 1.2: u2,u3
  CT1E2 <- wtd.kdecop(u[,c(2,3)], weights=w, method=method[2], mult=mult[2], knots=knots[2], bw=bw[2])
  # Edge 1.3: u3,u4
  CT1E3 <- wtd.kdecop(u[,c(3,4)], weights=w, method=method[3], mult=mult[3], knots=knots[3], bw=bw[3])
  
  ##########################################################################
  ### DVine Tree 2: 
  
  if(print_info==TRUE){cat("2...")}
  
  # Edge 2.1: u1,u3; u2 
  h1c <- hkdecop(u[,c(1,2)], CT1E1, cond.var=2) # Hfunc1:= h_2.1(u3|u1) and Hfunc2:= h_1.2(u1|u3) where index 1,2 
  h3c <- hkdecop(u[,c(2,3)], CT1E2, cond.var=1) 
  CT2E1 <- wtd.kdecop(as.copuladata(cbind(h1c,h3c)), weights=w, method=method[4], mult=mult[4], knots=knots[4], bw=bw[4])
  
  # Edge 2.2: u2,u4; u3 
  h2c <- hkdecop(u[,c(2,3)], CT1E2, cond.var=2)
  h4c <- hkdecop(u[,c(3,4)], CT1E3, cond.var=1)
  CT2E2 <- wtd.kdecop(as.copuladata(cbind(h2c,h4c)), weights=w, method=method[5], mult=mult[5], knots=knots[5], bw=bw[5])
  
  ##########################################################################
  ### DVine Tree 3: 
  
  if(print_info==TRUE){cat("3...")}
  
  # Edge 3.1:  u1,u4; u2, u3 
  h1cc <- hkdecop(cbind(h1c,h3c), CT2E1, cond.var=2) 
  h4cc <- hkdecop(cbind(h2c,h4c), CT2E2, cond.var=1) 
  CT3E1 <-  wtd.kdecop(as.copuladata(cbind(h1cc,h4cc)), weights=w, method=method[6], mult=mult[6], knots=knots[6], bw=bw[6])
  
  }else{
  #########################################################################
  ### Fit C-Vine if required
  #
  #T1: v1   <->   v2   |  v1  <->  v3  |  v1  <-> v4
  #T2: v2;v1 <->  v3;v1 | v2;v1 <->  v4;v1
  #T3: v3;v2,v1  <->  v4;v2,v1
    
  #########################################################################
  ### C-vine Tree 1: 
  if(print_info==TRUE){cat("Tree 1...")}
    
  # Edge 1.1: u1,u2
  CT1E1 <- wtd.kdecop(u[,c(1,2)], weights=w, method=method[1], mult=mult[1], knots=knots[1], bw=bw[1])
  # Edge 1.2: u2,u3
  CT1E2 <- wtd.kdecop(u[,c(1,3)], weights=w, method=method[2], mult=mult[2], knots=knots[2], bw=bw[2])
  # Edge 1.3: u3,u4
  CT1E3 <- wtd.kdecop(u[,c(1,4)], weights=w, method=method[3], mult=mult[3], knots=knots[3], bw=bw[3])
    
  ##########################################################################
  ### CVine Tree 2: 
  
  if(print_info==TRUE){cat("2...")}
  
  # Edge 2.1: u1,u3; u2 
  h2c <- hkdecop(u[,c(1,2)], CT1E1, cond.var=1) # Hfunc1:= h_2.1(u3|u1) and Hfunc2:= h_1.2(u1|u3) where index 1,2 
  h3c <- hkdecop(u[,c(1,3)], CT1E2, cond.var=1)
  CT2E1 <- wtd.kdecop(as.copuladata(cbind(h2c,h3c)), weights=w, method=method[4], mult=mult[4], knots=knots[4], bw=bw[4])
  
  # Edge 2.2: u2,u4; u3 
  h4c <- hkdecop(u[,c(1,4)], CT1E3, cond.var=1)
  CT2E2 <- wtd.kdecop(as.copuladata(cbind(h2c,h4c)), weights=w, method=method[5], mult=mult[5], knots=knots[5], bw=bw[5])
  
  ##########################################################################
  ### CVine Tree 3: 
  
  if(print_info==TRUE){cat("3...")}
  
  # Edge 3.1:  u1,u4; u2, u3 
  h3cc <- hkdecop(cbind(h2c,h3c), CT2E1, cond.var=1) 
  h4cc <- hkdecop(cbind(h2c,h4c), CT2E2, cond.var=1) 
  CT3E1 <-  wtd.kdecop(as.copuladata(cbind(h3cc,h4cc)), weights=w, method=method[6], mult=mult[6], knots=knots[6], bw=bw[6])
  
  }
  
  ##########################################################################
  ### Export results
  copulas <- list(CT1E1=CT1E1, CT1E2=CT1E2, CT1E3=CT1E3,
                  CT2E1=CT2E1, CT2E2=CT2E2, 
                  CT3E1=CT3E1)
  
  cAIC <- sum(unlist(lapply(copulas[1:6], function(x) x$info$cAIC)))
  
  res <- list(copulas=copulas,
              type=type,
              varnames=varnames,
              order=order,
              #Matrix=Matrix,
              cAIC=cAIC
  )
  return(res)
  
}

############################################################################
#### New wrapper function

VC_deco <- function(formula,data,group,weights,
                    ## Estimation
                    jit=c(1,3),            # adds jitter to variables
                    order=1:4,             # Order of variables; determines D-Vine or C-Vine structure
                    type="DVine",          # D-Vine or C-Vine structure?
                    mult=1,                # bandwidth multiplier
                    print_info=TRUE,       # print progress info
                    group0=NULL,           # define group 0
                    ## Find best tree
                    find_tree=FALSE,       # find optimal vine copula construction
                    nbest=5,               # number of vcc to be estimated
                    treecrit="tau",        # Use Kendall's tau oder Spearma's rho to select order?
                    ## Decomposition
                    deco=TRUE,             # performs deco
                    reference=1,           # if reference==1, copula of group 1 is used to analyze the impact of marginals
                    n=1000,                # number of simulated observations
                    tau=seq(5,95,5)/100,   # quantiles to estimated
                    multiply=(1/40),       # multiply simulated incomes by factor 
                    quasi=TRUE,            # quasi random numbers or pseudo random numbers?
                    copula_detail=TRUE,    # decompose detailed copula
                    export_usim=TRUE,      # export simulated quasi copula data
                    ## Bootstrap
                    bootstrap=FALSE,       # boostrap 
                    it=100,                # boostrap iterations
                    psu=NULL,              # psu 
                    stratum=NULL,          # stratum   
                    alpha=0.05,            # significance for uniform ci
                    ncore=1                # number of core for bootstrap
  ){
  
  #---------------------------------------------------------------
  ### 1 Prepare data
  #Use match.call to create data.frame
  mf = match.call()
  m = match(c("formula", "data", "weights","group"), names(mf), 0)
  mf = mf[c(1, m)]
  mf$drop.unused.levels = TRUE
  mf[[1]] = as.name("model.frame")
  
  #Add psu and stratum to formula in order to integrate it to model frame
  if(is.null(psu)==TRUE&is.null(stratum)==FALSE){
  mf$formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)], collapse=" ")," +",stratum)) 
  }else if(is.null(psu)==FALSE&is.null(stratum)==TRUE){
  mf$formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)], collapse=" ")," +",psu))   
  }else if(is.null(psu)==FALSE&is.null(stratum)==FALSE){
  mf$formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)], collapse=" ")," +",psu," +", stratum))   
  }
  
  # Evaluate model.frame
  mf = eval.parent(mf)
  mt = attr(mf, "terms")
  
  # Adjust formula
  formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)],collapse=" ")," - 1"))
  
  #Retrieve weights
  if(is.null(model.weights(mf))){
    mf$`(weights)` <- rep(1,nrow(mf))
  }
  
  # Select group 0
  if(is.null(group0)|length(which(mf[,"(group)"]==group0))){
    if(is.factor(mf[,"(group)"])){
      sel0 <- levels(mf[,"(group)"])[1]
    }else{
      sel0 <- sort(unique(mf[,"(group)"]))[1]
    }
  }else{
     sel0 <- group0
  }

  #drop observations with zero weight
  if(print_info==TRUE){
    cat(paste0(c("\nDropped because of 0 household weight:", length(which(mf[,"(weights)"]==0&mf[,"(group)"]==sel0)),"in group 0 and",
                 length(which(mf[,"(weights)"]==0&mf[,"(group)"]!=sel0)) ,"in group 1.\n")))
  }
  mf <- mf[which(mf[,"(weights)"]!=0),]
  
  # Select 
  data0 <- mf[which(mf[,"(group)"]==sel0),]
  data1 <- mf[which(mf[,"(group)"]!=sel0),] 
  
  #2 ---------------------------------------------------------------
  # Either estimate prespecified VCC or find optimal order
  if(find_tree==FALSE){
    VCfit <- VC_fit(data0,data1,jit,order,type, mult, print_info)
  }else{
    cat("Find best tree... \n")
    VCfit <- VC_find_tree(data0[,2:5],data1[,2:5],
                          data0[,"(weights)"],data1[,"(weights)"],
                          treecrit=treecrit,
                          jit=jit,nbest=nbest)
    #print_best_order(VCfit)  
  }
  
  VCfit <- c(VCfit,list(formula=formula))
  
  #3 ---------------------------------------------------------------
  # Perform decomposition
  if(deco){
    VCdeco <- VC_deco_sim(VCfit, reference=reference, 
                          multiply = multiply, tau=tau, print_info = print_info,
                          n=n, quasi = quasi, copula_detail = copula_detail)
  }else{
    VCdeco <- NULL  
  }
  
  #4 ---------------------------------------------------------------
  # Boostrap
  if(bootstrap){
    cat("\nBootstrapping...\n")
    VCdeco_se <- VC_deco_bs(VCfit,VCdeco,it=it, ncore=ncore,
                            psu=psu, stratum=stratum, alpha=alpha)  
  }else{
    VCdeco_se <- NULL  
  }
  
  #5 ---------------------------------------------------------------
  # Return results
  return(list(VCfit=VCfit,
              VCdeco=VCdeco,
              VCdeco_se=VCdeco_se))
}


############################################################################
#### Wrapper function to fit copulas in both groups
VC_fit <- function(data0,data1,jit,order,type, mult, print_info){

  # Add jittering to discrete data
  data0[,jit+1] <- apply(data0[,jit+1], 2, function(x) add_jitter(x))
  data1[,jit+1] <- apply(data1[,jit+1], 2, function(x) add_jitter(x))
  
  # Compute weighted percentile ranks
  u0 <- apply(data0[,2:5], 2, function(x) rank.wt(x,data0[,"(weights)"]))
  u1 <- apply(data1[,2:5], 2, function(x) rank.wt(x,data1[,"(weights)"]))
  
  # Store varnames
  varnames <- names(data0)[2:5]
  
  # Make sure to select only first 4 elements of order vector 
  order <- order[1:4]
  
  ### Find tree or estimate predetermined order?
  ## Estimate
  if(print_info==TRUE){cat("Copula fitting group 0...\n")}
    copula0 <- kdeVineCop(u0, w=data0[,"(weights)"], order=order, type=type, mult=mult)
    
    if(print_info==TRUE){cat("\nCopula fitting group 1...\n")}
    copula1 <- kdeVineCop(u1, w=data1[,"(weights)"], order=order, type=type, mult=mult)
    
    res <- list(copula0=copula0,
                copula1=copula1,
                order=order,
                type=type,
                data0=data0,data1=data1,
                u0=u0,u1=u1,
                varnames=varnames, jit=jit)
}

############################################################################
#### Wrapper function to fit copulas for both groups
#VC_deco_fit <- function(formula,data,group,weights,
#                        jit=c(1,3),       #adds jitter to variables
#                        order=1:4,        #Order of variables; determines D-Vine or C-Vine structure
#                        type="DVine",     #D-Vine or C-Vine structure?
#                        #reference=NULL,   #define reference group
#                        mult=1,           #bandwidth multiplier
#                        print_info=TRUE,  #print progress info
#                        find_tree=FALSE,  #find optimal bandwidth
#                        nbest=5           #number of orders/models to be estimated
#                        ){
#  
#  #Use match.call to create data.frame
#  mf = match.call()
#  m = match(c("formula", "data", "weights", "na.action","group"), names(mf), 0)
#  mf = mf[c(1, m)]
#  mf$drop.unused.levels = TRUE
#  
#  mf[[1]] = as.name("model.frame")
#  mf = eval.parent(mf)
#  mt = attr(mf, "terms")
#  
#  # Adjust formula
#  formula <- as.formula(paste0(paste0(as.character(formula)[c(2,1,3)],collapse=" ")," - 1"))
#  
#  #Retrieve outcome variable
#  #dep = model.response(mf, "numeric")
#  
#  #Retrieve weights
#  if(is.null(model.weights(mf))){
#    mf$`(weights)` <- rep(1,nrow(mf))
#  }
#  
#  # Set reference group 
#  if(is.null(reference)==FALSE){
#    sel0 <- reference
#  } else if(is.factor(mf[,"(group)"])){
#    sel0 <- levels(mf[,"(group)"])[1]
#  }else{
#    sel0 <- sort(unique(mf[,"(group)"]))[1]
#  }
#  
#  #drop observations with zero weight
#  if(print_info==TRUE){
#    cat(paste0(c("\nDropped because of 0 household weight:", length(which(mf[,"(weights)"]==0&mf[,"(group)"]==sel0)),"in group 0 and",
#               length(which(mf[,"(weights)"]==0&mf[,"(group)"]!=sel0)) ,"in group 1.\n")))
#  }
#  mf <- mf[which(mf[,"(weights)"]!=0),]
#  
#  # Select 
#  data0 <- mf[which(mf[,"(group)"]==sel0),]
#  data1 <- mf[which(mf[,"(group)"]!=sel0),]
#  
#  # Add jittering to discrete data
#  data0[,jit+1] <- apply(data0[,jit+1], 2, function(x) add_jitter(x))
#  data1[,jit+1] <- apply(data1[,jit+1], 2, function(x) add_jitter(x))
#  
#  # Compute weighted percentile ranks
#  u0 <- apply(data0[,2:5], 2, function(x) rank.wt(x,data0[,"(weights)"]))
#  u1 <- apply(data1[,2:5], 2, function(x) rank.wt(x,data1[,"(weights)"]))
#  
#  # Store varnames
#  varnames <- names(data0)[2:5]
#  
#  # Make sure to select only first four elements of order vector 
#  order <- order[1:4]
#  
#  ### Find tree or estimate predetermined order?
#  if(find_tree==FALSE){
#  
#  ## Estimate
#  if(print_info==TRUE){cat("Copula fitting group 0...\n")}
#  copula0 <- kdeVineCop(u0, w=data0[,"(weights)"], order=order, type=type, mult=mult)
#  
#  if(print_info==TRUE){cat("\nCopula fitting group 1...\n")}
#  copula1 <- kdeVineCop(u1, w=data1[,"(weights)"], order=order, type=type, mult=mult)
#  
#  # Compute average Log Likelihood
#  # logl0 <- sum(log(dKdeVineCop(u0, copula0))*data0[,"(weights)"])/sum(data0[,"(weights)"])
#  # logl1 <- sum(log(dKdeVineCop(u1, copula1))*data1[,"(weights)"])/sum(data1[,"(weights)"])
#  
#  res <- list(copula0=copula0,
#              copula1=copula1,
#              order=order,
#              type=type,
#              data0=data0,
#              data1=data1,
#              u0=u0,
#              u1=u1,
#              formula=formula,
#              varnames=varnames,
#              jit=jit
#              )
#  
#  }else{
#  cat("Finding best tree... \n")
#  res <- VC_find_tree(data0[,2:5],data1[,2:5],data0[,"(weights)"],data1[,"(weights)"],jit=jit,nbest=nbest)
#  res <- c(res, list(formula=formula))
#  }

#  # Print tree search if requested
#  if(find_tree==TRUE){print_best_order(res)}
  
#  ## Return results
#  return(res)
#} 



############################################################################
#### Function to find optimal tree
VC_find_tree <- function(data0, data1, 
                         w0=NULL, w1=NULL,
                         jit=c(1,3), nbest=5,
                         treecrit="tau"){
  
  # Get varnames
  vn <- names(data0)
  
  # Find all unique D-Vine combinations 
  orders <- rbind(find_unique_orders(n=4, type="DVine"),
                  find_unique_orders(n=4, type="CVine"))
  
  # Kendall's Tau
  #s0 <- unique (unlist (lapply (u0, function (x) which (is.na (x)))))
  #s0 <- setdiff(1:nrow(u0),s0)
  #s1 <- unique (unlist (lapply (u1, function (x) which (is.na (x)))))
  #s1 <- setdiff(1:nrow(u1),s1)
  #Ktau0 <- TauMatrix(u0[s0,], weights=w0[s0])
  #Ktau1 <- TauMatrix(u1[s1,], weights=w1[s1])
  
  # Compute unjittered rank variables
  u0 <- apply(data0,2, function(x) rank.wt(x,w0))
  u1 <- apply(data1,2, function(x) rank.wt(x,w1))
  
  if(treecrit=="tau"){
    Ktau0 <- TauMatrix(u0, weights=w0)
    Ktau1 <- TauMatrix(u1, weights=w1)
  }else{
    Ktau0 <- cov.wt(u0,w=w0, cor=TRUE)$cor
    Ktau1 <- cov.wt(u1,w=w1, cor=TRUE)$cor  
  }
    
  # Compute sum of tau by tree
  #EdgesByKtau0 <- t(apply(orders,1,function(x) c(Ktau0[x[1],x[2]],Ktau0[x[2],x[3]],Ktau0[x[3],x[4]])))
  #EdgesByKtau1 <- t(apply(orders,1,function(x) c(Ktau1[x[1],x[2]],Ktau1[x[2],x[3]],Ktau1[x[3],x[4]])))
  
  EdgesByKtau0 <- t(apply(orders,1,function(x) gettaus(x,Ktau0)))
  EdgesByKtau1 <- t(apply(orders,1,function(x) gettaus(x,Ktau1)))
    
  TauSumOverEdges0 <- rowSums(abs(EdgesByKtau0))
  TauSumOverEdges1 <- rowSums(abs(EdgesByKtau1))
  
  os0 <- order(TauSumOverEdges0, decreasing =TRUE)
  os1 <- order(TauSumOverEdges1, decreasing =TRUE)
  
  # Reorder
  best0 <- orders[os0,]
  best1 <- orders[os1,]
  
  besttau0 <- EdgesByKtau0[os0,]
  besttau1 <- EdgesByKtau1[os1,] 
  
  ## Fit copulas and store logLik of best models
  # Add jittering to discrete data
  u0[,jit] <- apply(data0[,jit],2, function(x) add_jitter(x))
  u1[,jit] <- apply(data1[,jit],2, function(x) add_jitter(x))
  u0[,jit] <- apply(u0[,jit],2, function(x) rank.wt(x,w0))
  u1[,jit] <- apply(u1[,jit],2, function(x) rank.wt(x,w1))
  
  #logl0 <- NULL
  #logl1 <- NULL
  
  cAIC0 <- NULL
  cAIC1 <- NULL
  
  copulas0 <- list()
  copulas1 <- list()

  nbest <- min(nbest,nrow(best0))
  cat("Estimating best models:")
  for(i in 1:nbest){
    cat(paste0("\n",i,"...\n"))
    copulas0[[i]] <- kdeVineCop(u0, w=w0, order=best0[i,1:4], type=ifelse(best0[i,5]==0,"DVine","CVine"))
    copulas1[[i]] <- kdeVineCop(u1, w=w1, order=best1[i,1:4], type=ifelse(best0[i,5]==0,"DVine","CVine"))
    
    # Store sum of cAIC
    cAIC0 <- c(cAIC0, copulas0[[i]]$cAIC)
    cAIC1 <- c(cAIC1, copulas1[[i]]$cAIC)
    
  }

 return(list(copula0=copulas0[[1]],
             copula1=copulas1[[1]],
             order=best0[1,1:4],
             type=ifelse(best0[1,5]==0,"DVine","CVine"),
             data0=data0,data1=data1,
             u0=u0,u1=u1,
             varnames=vn, jit=jit,
             best0=best0,
             best1=best1,
             besttau0 = besttau0,
             besttau1 = besttau1,
             cAIC0=cAIC0,
             cAIC1=cAIC1,
             treecrit=treecrit
            ))
  
}

############################################################################
#### Function to get tau of first tree DVine-CVine orders
gettaus <- function(x,tau){ 
                  if(x[5]==0){
                    # DVine
                    c(tau[x[1],x[2]],tau[x[2],x[3]],tau[x[3],x[4]])
                  }else{
                    # CVine
                    c(tau[x[1],x[2]],tau[x[1],x[3]],tau[x[1],x[4]])
                  }
}

############################################################################
#### Function to find unique order

find_unique_orders <- function(n,type=c("DVine","CVine")){
  orders_list <- list()
  n <- 1:n
  m=1
  for(i in n){
    for(j in setdiff(n,i)){
      for(k in setdiff(n,c(i,j))){
        for(l in setdiff(n,c(i,j,k))){
          orders_list[[m]] <- c(i,j,k,l)
          m=m+1
        }
      }
    }
  }
  
  reverse_orders_list <- lapply(orders_list, function(x) rev(x))
  orders_string <- unlist(lapply(orders_list,function(x) paste0(x,collapse="")))
  reverse_orders_string <- unlist(lapply(reverse_orders_list,function(x) paste0(x,collapse="")))
  
  orders_matrix <- do.call("rbind", orders_list)
  reverse_orders_matrix <- do.call("rbind", reverse_orders_list)
  
  if(type[1]=="DVine"){
    res <- NULL
    for(i in 1:length(orders_string)){
      # check if new combination not already reverse order permutation 
      # of existing combination
      if(is.element(orders_string[i],reverse_orders_string[1:(i-1)])==FALSE){
      res <- rbind(res,orders_matrix[i,])
      
    }
    }
    res <- cbind(res,0)
  }else{
    sel <- seq(1,nrow(orders_matrix)-1,2)
    res <- orders_matrix[sel,]
    res <- cbind(res,1)
  }
  return(res)
}


############################################################################
#### Print best order

print_best_order <- function(VCopObj){
  
  VCopObj <- VCopObj$VCfit
  
  best0 <- VCopObj$best0
  best1 <- VCopObj$best1
  besttau0 <- VCopObj$besttau0
  besttau1 <- VCopObj$besttau1
  cAIC0 <- VCopObj$cAIC0
  cAIC1 <- VCopObj$cAIC1
  vn <- VCopObj$varnames
  nbest <- length(cAIC0)
  treecrit <- VCopObj$treecrit
  
  cat("\n# Edge weight:",ifelse(treecrit=="tau","Kendall's tau","Spearman's rho"),"\n")
  cat("\n# Best orders: Group 0\n")
  for(i in 1:nbest){
    
    if(best0[i,5]==0){
    # D-Vine  
    cat("# ",i,"  ",
          "D-vine: ",
          vn[best0[i,1]],
          " <-(",round(besttau0[i,1],2),")-> ",
          vn[best0[i,2]],
          " <-(",round(besttau0[i,2],2),")-> ",
          vn[best0[i,3]],
          " <-(",round(besttau0[i,3],2),")-> ",
          vn[best0[i,4]]," | cAIC: ", round(cAIC0[i],1),
          " | c(", paste(best0[i,1:4],collapse = ","),")",
          "\n", sep="")  
    }else{
    # C-Vine  
      cat("# ",i,"  ",
          "C-vine: ",
          vn[best0[i,1]], 
          " <-> ", vn[best0[i,2]], "(",round(besttau0[i,1],2),")|",
          " <-> ", vn[best0[i,3]], "(",round(besttau0[i,2],2),")|",
          " <-> ", vn[best0[i,4]], "(",round(besttau0[i,3],2),")",
          " | cAIC: ", round(cAIC0[i],1),
          " | c(", paste(best0[i,1:4],collapse = ","),")",
          "\n", sep="")    
    }
    
   
  }
  
  cat("\n# Best orders: Group 1\n")
  for(i in 1:nbest){
    if(best1[i,5]==0){
      cat("# ",i,"  ",
          "D-vine: ",
          vn[best1[i,1]],
          " <-(",round(besttau1[i,1],2),")-> ",
          vn[best1[i,2]],
          " <-(",round(besttau1[i,2],2),")-> ",
          vn[best1[i,3]],
          " <-(",round(besttau1[i,3],2),")-> ",
          vn[best1[i,4]]," | cAIC: ", round(cAIC1[i],1),
          " | c(", paste(best1[i,1:4],collapse = ","),")",
          "\n", sep="")  
    
  }else{
    # C-Vine  
    cat("# ",i,"  ",
        "C-vine: ",
        vn[best1[i,1]], 
        " <-> ", vn[best1[i,2]], "(",round(besttau1[i,1],2),")|",
        " <-> ", vn[best1[i,3]], "(",round(besttau1[i,2],2),")|",
        " <-> ", vn[best1[i,4]], "(",round(besttau1[i,3],2),")",
        " | cAIC: ", round(cAIC1[i],1),
        " | c(", paste(best1[i,1:4],collapse = ","),")",
        "\n", sep="")    
  }
    
  }
  
}
