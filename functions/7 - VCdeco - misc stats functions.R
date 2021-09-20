###############################################################
###############################################################
### 7 - VCdeco - misc statistic functions

###############################################################
# Saving stargazer/xtable output in text files
mod_xtable <- function(output.file, ...) {
  output <- capture.output(print(...))
  cat(paste(output, collapse = "\n"), "\n", file=output.file, append=FALSE)
}

###############################################################
# Compute Gini
gini2   <- function(x,weights=NULL,trim=0){
  if(is.null(weights)){weights=rep(1,length(x))}
  if(sum(as.numeric(is.na(x)))>0){
    weights <- weights[-which(is.na(x))] 
    x <- x[-which(is.na(x))] 
  }
  if(trim!=0){
    sel <- which(x>=wtd.quantile(x,weights,probs=trim)&x<=wtd.quantile(x,weights,probs=1-trim))
    x <- x[sel]
    weights <- weights[sel]
  }
  mu <- weighted.mean(x,w=weights) 
  #ranks <- wtd.rank(x, weights=weights, normwt=TRUE)
  #ranks <- ranks/max(ranks, na.rm=TRUE)
  ranks <- rank.wt(x, weights)
  gini <- 2*cov(cbind(x,ranks), use = "pairwise.complete.obs")[1,2]/mu
  return(gini)
}

###############################################################
# Estimate "pseudo copula data" (i.e. weighted normalized ranks on (0,1)-interval)
rank.wt <- function(x, weights=NULL, left_limit_CDF=FALSE, ties_method="average") {
  if(is.null(weights)){weights=rep(1,length(x))}
  na <- which(is.na(x))
  sw <- sum(weights) - sum(weights[na]) + quantile(weights, probs=0.5)  #add median weightd instead n+1 in order to avoid problems at the margin
  weights[na] <- NA
  o <- order(x)
  wr <- cumsum(weights[o])/sw
  r <- rank(x, ties.method = ties_method) # max #Rank corresponds to P(X<=x): For ties highest rank is selected 
  x <- wr[r]  
  if(left_limit_CDF){ 
    #Left-sided limit of CDF of discrete variables, i.e. F_X(x^-)=F_X(x-1)
    #r <- rank(x, ties.method = "min")-1
    x0 <- x
    r <- rank(x0, ties.method = "min")-1
    x0 <- x0[order(x0)]
    x <- rep(0,length(x0))
    sel <- which(r==0)
    #x[-sel] <- wr[r[-sel]]
    x[-sel] <- x0[r[-sel]]
  }
  if(ties_method=="max"){x[x>0.998] <- 0.998} #trim
  return(x)  
} 

#rank.wt2 <- function(x, weights=NULL) {
#  if(is.null(weights)){weights=rep(1,length(x))}
#  sw <- sum(weights) - sum(weights[which(is.na(x))]) + quantile(weights, probs=0.5)  #add median weightd instead n+1 in order to avoid problems at the margin
#  x <- sapply(x, function(y)  if(is.na(y)){NA}else{sum(weights[which(x<=y)])/(sw)})
#  return(x)  
#} 


###############################################################
# Estimate share of observations with y=0
share0 <- function(x,weights) (sum(weights[which(x==0)])/sum(weights))

###############################################################
# Wrapper function for weigthed summary statistics
wtd.sum <- function(x, weight, 
                    na.rm=TRUE, 
                    trim=0, log_trans=TRUE) { 
  x_t <- if(log_trans){log(x)}else{x}
  mu <- weighted.mean(x,weights = weight, na.rm=na.rm)
  sd <- sqrt(wtd.var(x,weights = weight, na.rm=na.rm))
  res <- c(if(log_trans){log(mu)}else{mu},
           if(log_trans){log(sd)}else{sd},
           min(x,na.rm=na.rm),
           p10 <- wtd.quantile(x_t,weights = weight,probs=0.1, na.rm=na.rm),
           p25 <- wtd.quantile(x_t,weights = weight,probs=0.25, na.rm=na.rm),
           p50 <- wtd.quantile(x_t,weights = weight,probs=0.5, na.rm=na.rm),
           p75 <- wtd.quantile(x_t,weights = weight,probs=0.75, na.rm=na.rm),
           p90 <- wtd.quantile(x_t,weights = weight,probs=0.9, na.rm=na.rm),
           p99 <- wtd.quantile(x_t,weights = weight,probs=0.99, na.rm=na.rm),
           max(x,na.rm=na.rm),
           share0(x,weight), 
           if(log_trans){p90-p10}else{p90/p10},
           if(log_trans){p50-p10}else{p50/p10},
           if(log_trans){p90-p50}else{p90/p50},
           gini2(x, weight,trim=trim),
            #theil(x, weight),
           sum(as.numeric((is.na(x)==FALSE)))
)
names(res) <- c("Mean","Std. Dev.","Min.",
                "P10","P25","P50","P75","P90","P99",
                "Max.","%obs.=0", 
                "P90/P10","P50/P10", "P90/P50", 
                "Gini",#"Theil L", "Theil T",
                "N")
return(res)
}

###############################################################
# Weighted correlation function
wtd.cor <- function(data,weights,trim=0){
  if(dim(data)[2]<2){stop("Two variables min!")}
  if(is.null(weights)){weights=rep(1,nrow(data))}
  #Trimming for Pearson Correlation
  if(length(trim)==1&trim[1]==0){
   sel <- 1:nrow(data)
  }else{
   if(length(trim)!=ncol(data)){trim <- rep(trim[1],ncol(data))}
   trim_thresholds <- mapply(function(x,y) c(wtd.quantile(x,weights,probs=y), wtd.quantile(x,weights,probs=1-y)),
                             x=as.list(data), y=trim)
   sel <- mapply(function(x,y) which(x>=y[1]&x<=y[2]),
                    x=as.list(data), y=as.list(as.data.frame(trim_thresholds)))
   sel <- Reduce(intersect,sel)
   if(length(sel)==0){stop("Too much trimming!")}
  }
  
  #Coefficient estimation
  u <- apply(data,2, function(x) rank.wt(x,weights))
  pearson <- cov.wt(data[sel,], wt=weights[sel], cor=TRUE)$cor
  spearman <- cov.wt(u, wt=weights, cor=TRUE)$cor
  pearson <- pearson[lower.tri(pearson)]
  spearman <- spearman[lower.tri(spearman)]
  
  # Get variable names
  n <- names(data)
  k <- length(n)
  vcomb <- matrix(paste0(matrix(n,k,k,byrow=TRUE),matrix("-",k,k),matrix(n,k,k,byrow=FALSE)),k,k,byrow=TRUE)
  vcomb <- vcomb[lower.tri(vcomb)]
  
  # Return results
  res <-  data.frame(vars=rep(vcomb,2),
                        rho=c(pearson,spearman),
                        type=rep(c("Pearson","Spearman"),each=length(pearson)))
  return(res)
}


###############################################################
# Theil index

theil <- function(x,weights=NULL){ 
         if(is.null(weights)){weights=rep(1,length(x))}
         mu <- mean(x, na.rm=TRUE)
         Tl <- weighted.mean(x/mu*ifelse(is.infinite(log(x/mu)),0,log(x/mu)), weights=weights, na.rm=TRUE)
         Tt <- weighted.mean(ifelse(is.infinite(log(mu/x)),0,log(mu/x)), weights=weights, na.rm=TRUE)
         return(c(Tl,Tt))
        }


###############################################################
# Compute conditional ecdf of a copula
CondCop <- function(u1,u2,u=0.1,w=NULL,q=10){
  if(is.null(w)){w <- rep(1,length(u2))}
  cat <- cut(u1,(0:q)/q, right=FALSE)
  unlist(lapply(split(data.frame(u2,w),cat), 
                function(x) weighted.mean(as.numeric(x[,1]<=u), w=x[,2])))
}


############################################################################
#### Function to select trimmed samples
trims <- function(data,weights,trim){
  if(length(trim)!=ncol(data)){trim <- rep(trim[1],ncol(data))}
  trim_thresholds <- mapply(function(x,y) c(wtd.quantile(x,weights,probs=y), wtd.quantile(x,weights,probs=1-y)),
                            x=as.list(data), y=trim)
  sel <- mapply(function(x,y) which(x>=y[1]&x<=y[2]),
                x=as.list(data), y=as.list(as.data.frame(trim_thresholds)))
  sel <- Reduce(intersect,sel)
  if(length(sel)==0){stop("Too much trimming!")}
  return(sel)
}
