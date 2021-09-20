library("AER")
library(rvinecopulib)
setwd("C:/Users/David Gallusser/switchdrive/Uni Basel/Dissprojekte/6 Couples earnings decomposition/R code and data/functions")
source("4 - VCdeco - misc stats functions.R")
setwd("C:/Users/David Gallusser/switchdrive/Uni Basel/R Code/LSE/functions/rif-dfl-deco")
source("r-DFL-deco-2-2.R")

###############################################################
### Working directories

# Uni
setwd("C:/Users/gallusse/switchdrive/Uni Basel/Dissprojekte/6 Couples earnings decomposition/R code and data/CPS/")
# Laptop
setwd("C:/Users/David Gallusser/switchdrive/Uni Basel/Dissprojekte/6 Couples earnings decomposition/R code and data/CPS/")
sav_d <- getwd()

# Path to data directory
datafiles_d <- "/data/" 


##### Load CPS data
# Setwd
setwd(paste0(sav_d,datafiles_d))
years <- c(1976,2016)
cdataall <- NULL
for(i in years){
  load(paste0("CPScouples",i,".rda"))
  cdataall <- rbind(cdataall,cdata)
}
cdata <- cdataall
rm("cdataall")

names(cdata)

breaks <- c(0,1,5,10,15,20,30,60)
cdata$EXP2_f <- cut(cdata$EXP_f, breaks=breaks)
cdata$EXP2_m <- cut(cdata$EXP_m, breaks=breaks)

sort(unique(cdata$EXP2_m))
f1 <- EARNINGS ~ EDUC_m*EDUC_f + EDUC_m*EXP2_m + EDUC_f*EXP2_f + EXP2_m*EXP2_f +  
                 RACE_f*RACE_m + HISPAN_f*HISPAN_m + METRO_m*REGION_m 

cdata$year <- as.numeric(cdata$YEAR_m==2016)     
res <- dfl_deco(f1, cdata, group=YEAR_m, weights=ASECWTH_m, reference=0, tau=1:9/10)
ncol(model.matrix(f1,cdata))
dfl_deco_plot(res)

cdata$BLACK_m <- as.numeric(cdata$RACE_m=="200")
cdata$BLACK_f <- as.numeric(cdata$RACE_f=="200")
cdata$OTHER_m <- as.numeric(cdata$RACE_m=="300")
cdata$OTHER_f <- as.numeric(cdata$RACE_f=="300")
f2 <- EARNINGS ~  EDUC_m + EXP_m + BLACK_m + OTHER_m + HISPAN_m + 
                  EDUC_f + EXP_f + BLACK_f + OTHER_f + HISPAN_f 
res2 <- VC_deco_rw(f2, cdata, weights="ASECWTH_m", group="year")
dfplot <- tidyr::pivot_longer(res2, names(res)[2:ncol(res)], names_to="effect")
ggplot(dfplot, aes(tau, value, col=effect)) + geom_line() + geom_point()


##### Load CPS1985 data
data("CPS1985")
x <- CPS1985$education

p <- sampleProb(x)
predict.sampleProb(p,x)

res <- VC_deco_rw(log(wage) ~ experience + education + ethnicity, CPS1985, group="union")
dfplot <- tidyr::pivot_longer(res, names(res)[2:ncol(res)], names_to="effect")
ggplot(dfplot, aes(tau, value, col=effect)) + geom_line() + geom_point()

res <- dfl_deco(log(wage) ~ experience * education * ethnicity, CPS1985, group=union, reference=0, tau=1:9/10)
dfl_deco_plot(res)
dfl_diag(res)
sampleProb <- function(x, w=NULL){
          if(is.null(w)){w=rep(1,length(x))}
          u <- unique(x)
          u <- u[order(u)]
          p <- sapply(u, function(y) weighted.mean(as.numeric(x==y),w=w))
          q <- c(0,cumsum(p[1:(length(p)-1)]))
          p <- data.frame(u,p,q)
          return(p)
}


predict.sampleProb <- function(p, newdata){
         return(p[match(newdata,p$u),"p"])
}

wtd.quantile.ordered <- function(x, w=NULL, probs=c(0.25,0.5,0.75)){
  if(is.null(w)){w=rep(1,length(x))}
  u <- unique(x)
  u <- u[order(u)]
  p <- sapply(u, function(y) weighted.mean(as.numeric(x==y),w=w))
  p <- c(0,cumsum(p[1:(length(p)-1)]))
  q <- u[sapply(probs, function(y) max(which(y>p)))]
  names(q) <- paste0(probs*100,"%")
  return(q)
}



VC_deco_rw <- function(formula, data, weights=as.numeric(), group,
                         var_types = "d",
                         family_set="tll",      # specify parametric or nonparaemtric copula modell 
                         structure=NA,          # copula matrix; if NA, algorithm searches for best tree
                         nonpar_method="quadratic", #est. method nonparametric models: "quadratic"= std. trans estimator with local likelihood approx. of order two
                         mult=1,                # bandwidth multiplier
                         print_info=TRUE,       # print progress 
                         group0=NULL,           # define group 0
                         tree_crit="rho",       # Criterion for tree selection
                         cores=1                # number of cores
                         ){
  
  #0 ---------------------------------------------------------------
  # Prepare data
  # Define weights
  if(length(weights)==0){
    data$weights <- rep(1,nrow(data))
  }else{
    data$weights <- data[,weights]
  }
  # Get group
  data$group <- data[,group]

  # Get model.frame
  #mf <- get_all_vars(formula,data,weights=weights,group=group)
  mf <- model.frame(formula,data,weights=weights,group=group)
  # Rename weights and group variable
  names(mf)[is.element(names(mf),c("weights","group"))] <- c("(weights)","(group)")
  
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
  
  # Select groups
  data0 <- mf[which(mf[,"(group)"]==sel0),]
  data1 <- mf[which(mf[,"(group)"]!=sel0),] 
  
  # get dimension of copula and varnames variables
  d <- length(all.vars(formula))-1 #ncol(data0)-3
  varnames <- all.vars(formula)[2:(d+1)] 
  
  # Add data imputation here
  impute=NULL
  
  #1 ---------------------------------------------------------------
  # Estimate pseudo-copula obs (i.e. F(x))
  if(length(var_types)==1){
    var_types <- rep(var_types,d)
  }
  
  u0 <- apply(data0[,2:(d+1)], 2, function(x) rank.wt(x,data0[,"(weights)"], ties_method = "max"))
  u1 <- apply(data1[,2:(d+1)], 2, function(x) rank.wt(x,data1[,"(weights)"], ties_method = "max"))
  
  if(any(var_types=="d")){
    # Compute left-sided limit of discrete/integer variables (i.e. F(x^-)=F(x-1))
    sel <- which(var_types=="d")+1
    u0minus <- apply(as.data.frame(data0[,sel]), 2, function(x) rank.wt(x,data0[,"(weights)"], ties_method = "max", left_limit_CDF = TRUE))
    u1minus <- apply(as.data.frame(data1[,sel]), 2, function(x) rank.wt(x,data1[,"(weights)"], ties_method = "max", left_limit_CDF = TRUE))
    colnames(u0minus) <- colnames(u1minus) <- paste0(colnames(u0minus),"_minus")
    u0 <- cbind(u0,u0minus)
    u1 <- cbind(u1,u1minus)
  }
  
  #2 ---------------------------------------------------------------
  # Fit copula model
  structure0 <- structure
  structure1 <- structure
  show_trace <- print_info
  estcores <- cores
  
  if(print_info){cat("Density estimation group 0... \n")}
  VCfit0 <- vinecop(u0, var_types=var_types, family_set = family_set, 
                    nonpar_method = nonpar_method,
                    structure=structure0, mult=mult, 
                    weights=data0[,"(weights)"], tree_crit=tree_crit,
                    cores=estcores, show_trace=show_trace)
  
  if(print_info){cat("\nDensity estimation group 1... \n")}
  VCfit1 <- vinecop(u1, var_types=var_types, family_set = family_set, 
                    nonpar_method = nonpar_method,
                    structure=structure1, mult=mult, 
                    weights=data1[,"(weights)"], tree_crit=tree_crit,
                    cores=estcores, show_trace=show_trace)
  
  
  #3 ---------------------------------------------------------------
  # Fit marginal "densities"
  Mfit0 <- lapply(data0[,2:(d+1)], function(x) sampleProb(x))
  Mfit1 <- lapply(data1[,2:(d+1)], function(x) sampleProb(x))
  
  
  
  #4 ---------------------------------------------------------------
  # Perform decomposition
  VCfit <- list(VCfit0=VCfit0, VCfit1=VCfit1,
                Mfit0=Mfit0, Mfit1=Mfit1,
                formula=formula,
                data0=data0, data1=data1,
                u0=u0, u1=u1,
                d=d,
                impute=impute)
  
  deco_res <- VC_deco_rw_predict(VCfit, cores)
  return(deco_res)
}

VC_deco_rw_predict <- function(VCfit, cores=1, tau=1:9/10){
  
  d <- VCfit$d
  X0 <- do.call("list",VCfit$data0[,2:(d+1)])
  X1 <- do.call("list",VCfit$data1[,2:(d+1)])
  
  # Predict copula densities
  cd0 <- dvinecop(VCfit$u0, VCfit$VCfit0, cores = cores)
  #cd1 <- dvinecop(VCfit$u1, VCfit$VCfit1, cores = cores)
  cdC <- dvinecop(VCfit$u0, VCfit$VCfit1, cores = cores)
  
  # Predict marginal densities
  md0 <- mapply(function(x,y) predict.sampleProb(x,y),
                                 VCfit$Mfit0, X0)
  # md1 <- mapply(function(x,y) predict.sampleProb(x,y),
  #                                VCfit$Mfit1, X1)
  mdC <- mapply(function(x,y) predict.sampleProb(x,y),
                                  VCfit$Mfit1, X0)
   
  # Reweighting factor
  psiC <- (apply(mdC,1,prod)*cdC) / (apply(md0,1,prod)*cd0)
  
  # Sum stats 
  y0 <- wtd.quantile(VCfit$data0[,1], weights= VCfit$data0[,"(weights)"], probs=tau)
  y1 <- wtd.quantile(VCfit$data1[,1], weights= VCfit$data1[,"(weights)"], probs=tau)
  yC <- wtd.quantile(VCfit$data0[,1], weights= VCfit$data0[,"(weights)"]*psiC, probs=tau)
  
  # Return results 
  res <- data.frame(tau=tau, 
                    observed=y1-y0, 
                    wage_structure=y1-yC, 
                    composition=yC-y0)
  return(res)
}



