########################################################
## 10 - VCdeco - topcoded income imputation functions
## Functions to impute CPS top incomes
########################################################

# Function to perfom random imputation of top coded individuals
# i.e. assign random draws from Pareto distribution
impute_top_inc <- function(df,year,cell_means){
  pos <- which(cell_means$YEAR==year)
  var <- cell_means[pos,"VARIABLE"]
  xmin <- cell_means[pos,"xmin"]
  m <- cell_means[pos,"MEAN_WAGE"]
  alpha <- cpareto(m,xmin)[1]
  sel <- which(df[,var]>=xmin)
  df[sel,var] <- rpareto(length(sel),xmin,alpha)
  if(var=="INCLONGJ"){
    df[sel,"INCWAGE"] <- ifelse(df[sel,"SRCEARN"]==1,df[sel,"INCLONGJ"],0) + df[sel,"OINCWAGE"]
  }
  return(df)
}

# Pareto distribution function
ppareto <- function(x,xmin,alpha){
  1-(xmin/x)^alpha
}

# Random draws from pareto (inverse transform sampling)
rpareto <- function(n,xmin,alpha){
  u <- runif(n)
  x <- xmin*(1-u)^(-1/alpha)
  return(x)
}

# Compute pareto coefficient
cpareto <- function(m,xmin){
  if(is.na(m)|is.na(xmin)){return(NA)}
  alpha <- m/(m-xmin)
  beta <- m/xmin
  return(c(alpha=alpha, beta=beta))
}

# MLE Hill estimator for pareto coefficient with top censored data
# see Armour, Burkhauser and Larrimore (2016) 
hillest <- function(y,w,taumin=0.99,xmin="max"){
  if(xmin=="max"){
    xt <- max(y,na.rm=TRUE)
  }else{
    xt <- xmin
  }
  xc <- wtd.quantile(y,weights=w,probs=taumin)
  Tt <- sum(w[which(y>=xt)])
  sel <- which(y<xt&y>=xc)
  M <- sum(w[sel])
  Sct <- sum(w[sel]*log(y[sel]))
  alpha <- M/(Tt*log(xt) + Sct - (M+Tt)*log(xc))
  return(alpha=alpha)
}
