##########################################################################
#### 11 - VCdeco -  estimation using kdvinecop
##########################################################################

require("kdevine")    # Non-parametric vine copula modeling
require("kdecopula")  # Copula density estimation

##########################################################################
#### Create copula matrix
create_copula_matrix <- function(x,type="Dvine"){
    d <- length(x)
    m <- matrix(0,d,d)
    rv <- NULL
    if(type=="Dvine"){
      diag(m)<- x
        for(i in 2:d){
          rv <- c(rv,x[d:i])
        }
    }else{
      diag(m)<- rev(x)
        for(i in (d-1):1){
        rv <- c(rv,x[i:1])
      }
    }
    m[lower.tri(m)] <- rv
    return(m)
}

##########################################################################
#### Jittering function
# Jitter dist. like in Example 1 of Nagler (2018: 5)
# unifom with theta=0, alternative: theta = 0.3
add_jitter <- function(x, theta=0, nu=5, bd=0.5) {
  n <- length(x)
  x <- round(x)+runif(n,-bd,bd)+theta*(rbeta(n,nu,nu)-0.5)
  return(x)
}

############################################################################
#### New wrapper function

VC_deco <- function(formula,data,group,weights,
                    ## Estimation
                    jit=c(1,3),            # which variables should be jittered?
                    matrix=NA,             # copula matrix; if NA, algorithm searches for best tree
                    type="DVine",          # D-Vine or C-Vine structure?
                    mult=1,                # bandwidth multiplier
                    print_info=TRUE,       # print progress info
                    group0=NULL,           # define group 0
                    treecrit="tau",        # Use Kendall's tau or Spearma's rho to select order?
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
                    cores=1                # number of cores
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
  
  # Add jittering to discrete data
  data0[,jit+1] <- apply(data0[,jit+1], 2, function(x) add_jitter(x))
  data1[,jit+1] <- apply(data1[,jit+1], 2, function(x) add_jitter(x))
  
  # dimension
  d <- ncol(data0)-3

  # Compute weighted percentile ranks
  u0 <- apply(data0[,2:(d+1)], 2, function(x) rank.wt(x,data0[,"(weights)"]))
  u1 <- apply(data1[,2:(d+1)], 2, function(x) rank.wt(x,data1[,"(weights)"]))
  
  
  #2 ---------------------------------------------------------------
  # Either estimate prespecified VCC or find optimal order
  if(print_info=TRUE){cat("VC density estimation group 0... \n")}
  VCfit0 <- wtd.kdevinecop(u0, matrix=matrix, weights=data0[,"(weights)"],
                                    treecrit=treecrit, info=TRUE, mult=mult,
                                    cores=cores)
    
  if(print_info=TRUE){cat("VC density estimation group 1... \n")}
  VCfit1 <- wtd.kdevinecop(u1, matrix=matrix, weights=data1[,"(weights)"],
                                     treecrit=treecrit, info=TRUE, mult=mult,
                                     cores=cores)
  
  #3 ---------------------------------------------------------------
  # Perform decomposition
  if(deco){
    VCdeco <- VC_deco_sim(VCfit0, VCfit1, formula=formula, jit=jit, 
                          data0=data0, data1=data1, reference=reference, 
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

