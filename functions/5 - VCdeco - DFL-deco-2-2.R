############################################
###### DFL decomposition function  #########
############################################

## Version 2.2, 26. August 2020

## Note: 
## Group 1 is the reference group. Its wage structure 
# is used for the counterfactuals.
# In order to sequentially decompose the composition effect, 
# variable have to be entered seperated by |
# If sequence="marginal" the marginal of the last variable entered
# is reweighted first.
# If firstrw="conditional" the conditional distribution of the first
# variable entered in the formula is reweigted first.

# library(AER)
# data("CPS1985")
# f <- wage ~ gender | experience | education | region

# result <- dfl_deco(f, data=CPS1985, group=union, tau=c(0.5,0.75,0.8))
# result$quantile
# result$other.stats

###########################################################
## Package to be loaded
require(Formula) #extended formulas
require(Hmisc) #weighted distributional statistics (i.e. quantiles, etc.)
require(survey) #for glm models with survey data
require(ggplot2) #plots
require(reshape2) #melting data.frames

###########################################################
# Actual decomposition function  

dfl_deco <- function(formula, data, weights, group, 
                     ws_formula=NULL,
                     reference=1, 
                     na.action = na.exclude, 
                     tau=c(10,50,90,99)/100, 
                     firstrw=c("marginal","conditional"),
                     stats=TRUE){
  
  # group: indicator variable defining the distribution to be compared
  # reference: defining the group which will be reweighted, i.e. which wage structure is to be used. 
  #            Composition effect: Difference between observed reference group and reweighted reference group
  #            Wage structure effect: Difference between reweighted reference group and comparison group
  # tau: quantiles to be evaluated
  # firstrw=conditional: sequential decompositon by reweighting conditional distribution
  #            of first entered covariate first, otherwise marginal distribution of
  #            last covariate is reweigted fist. 

  ##########################################################
  ##########################################################
  ## 1) Set up data
  
  #Use match.call function to call data.vectors
  mf = match.call()
  m = match(c("formula", "data", "weights", "na.action","group"), names(mf), 0)
  mf = mf[c(1, m)]
  mf$drop.unused.levels = TRUE
  
  # Retrieve formula as.Formula and model.frame
  f <- as.Formula(formula)
  
  # Make sure that all variables (incl. for detailed WS effect are selected)
  if(!is.null(ws_formula)){
    fws <- Formula(ws_formula)
    f <- update(f,as.Formula(paste(". ~ . + ", as.character(fws, lhs=1)[2], "+", as.character(fws, rhs=1)[3])))
  }
  
  mf[[1]] = as.name("model.frame")
  mf$formula <- f
  mf = eval.parent(mf)
  mt = attr(mf, "terms")
  
  # Store "orignal" decomposition model again
  f <- as.Formula(formula)
  
  # Extract variables, weights and group identifier
  #reg = get_all_vars(mt, mf)
  dep = model.response(mf, "numeric")
  weight = model.weights(mf)
  if (!is.null(weight) && !is.numeric(weight)) {
    stop("'weights' must be a numeric vector")
  }
  if (!is.null(weight)) {
    weight = weight
  }
  else {
    weight = rep(1, length(dep))
  }
  
  groupN = mf[, ncol(mf)]
  mf$groupN <- groupN
  
  # Number of covariates for decomposition of composition effect
  nvar = length(f)[2]
  
  # Check if group variable is a binary variable
  if(is.numeric(groupN)&length(unique(groupN))>2|
     is.factor(groupN)&length(unique(groupN))>2){
    stop("Group variable must either be a binary numeric or factor variable.")
  }
  
  # If binary numeric variable with value different than c(0,1) create factor
  if(is.numeric(groupN)&identical(sort(unique(groupN)),c(0,1))==FALSE){
       rvalue <- range(groupN, na.rm=TRUE)
       groupN <- factor(as.numeric(groupN==rvalue[2]))
       levels(groupN) <- rvalue
       mf$groupN <- groupN
       if(reference %in% c(0,1) == FALSE){
         reference <- which(rvalue==reference)-1
       }
  }
  
  # Set reference group to 1 if not determined
  if(is.character(reference)&is.factor(groupN)){
    reference <- which(levels(groupN)==reference)-1
  }
  if(is.factor(groupN)){
    reference_print <- levels(groupN)[reference+1]
  }else{
    reference_print <- reference
  }
  cat("Reweighted reference group: ", reference_print,"\n \n") 
  
  # If first reweighting not specified choose marginal first as default
  firstrw = ifelse(length(firstrw)==2,"marginal",firstrw)
  
  # Set Progress Bar 
  cat("Logit estimation...\n")
  pb <- txtProgressBar(min=0,max=1, style=3)
  nmods <- 1+nvar+ifelse(!is.null(ws_formula),4,0)
  
  ##########################################################
  ##########################################################
  ### 2) Fit logit models & estimate reweighting factors
  
  ##########################################################
  ### Unconditional/Sample probabilities
  mod <- groupN ~ 1
  p1 <- mean(pfit(mod,mf,weight))
  p0 <- 1-p1
  probs <- rep(p0/p1,nrow(mf))
  setTxtProgressBar(pb,1/nmods) 
  
  ##########################################################
  ### Conditional probabilities
  ### Beginning with the LAST variable entered in Formula,
  ### estimating then with prob
  
  m <- 0 # only for ProgressBar
  
  for(i in nvar:1){
    
  mod <- update(formula(f, rhs=nvar:i, collapse=TRUE), groupN ~ .)
  p1 <- pfit(mod,mf,weight)
  p0 <- 1-p1
  probs <- cbind(probs,p0/p1) 

  m <- m + 1 #Progressbar
  setTxtProgressBar(pb,(1+m)/nmods) 
  }
  
  ##########################################################
  ### Reweigthing factors
  ###
  ### psi_s contains relative probabilites: 
  ### first [P(t=1)/P(P=0)]*[P(t=1|X1)/P(P=0|X1)],
  ### second [P(t=0|X1)/P(P=1|X1)]*[P(t=1|X2,X1)/P(P=0|X2,X1)]
  ### etc.
  ### X1 is the variable entered last, X2 the variable entered second last etc.
  
  psi_s <- (probs[,1]^-1)*probs[,2]
  if(nvar==1){
  #if there is only one variable stick to psi_s and use it as weight
  psi <- as.matrix(psi_s)
  }
  else{
  #procedure to compute the weights if nvar>1
  for(i in 2:nvar){
    psi_s <- cbind(psi_s,(probs[,i]^-1)*probs[,i+1])
  }
  
  # Marginal of last entered variable or conditional of first entered 
  # variable to be reweighted first?
  first <- ifelse(firstrw=="marginal",1,nvar)
  last <- ifelse(firstrw=="marginal",nvar,1)
  correct <- ifelse(firstrw=="marginal",1,-1)
  loopvals <- seq(first+correct,last)
  
  # psi contains the weights actually used to reweighted the distribution
  # psi contains on the first position the weight used first
  # (i.e. [P(t=1)/P(P=0)]*[P(t=1|X1)/P(P=0|X1)] if the "marginal" and
  # [P(t=1|X1,...,X(M-1))/P(P=0||X1,...,X(M-1))]*[P(t=1|X1,...,X(M))/P(P=0|X1,...,X(M))]
  # if "conditional")
  psi <- as.matrix(psi_s[,first])
  colnames(psi) <- paste("psi_X",first,sep="")
  
  for(i in loopvals){
    psi <- cbind(psi,apply(psi_s[,first:i],1,prod))
    colnames(psi)[length(colnames(psi))] <- paste("psi_X",first,"to",i,sep="")
  }
  
  #end weight competion procedure if nvar>1
  }
  
  ###########################################################
  ### Weights for decomposition of wage structure effect
  
  if(!is.null(ws_formula)){
    
  # What's the group variable value of reference group?
  if(is.factor(groupN)){
    group_val <- levels(groupN)[1+reference]
  } else {
    group_val = reference
  }
    
  # Set up model for unconditional prob
  mod <- as.Formula(paste(as.character(fws, lhs=1)[2],"~ 1"))
  p_S_1 <- pfit(mod,mf[which(groupN==group_val),],weight[which(groupN==group_val)])  
  p_S_0 <- pfit(mod,mf[which(groupN!=group_val),],weight[which(groupN!=group_val)]) 
  
  setTxtProgressBar(pb,(nmods-2)/nmods)
  # Model for conditional prob
  mod <- fws 
  p_S_X_1 <- pfit(mod,mf[which(groupN==group_val),],weight[which(groupN==group_val)]) 
  p_S_X_0 <- pfit(mod,mf[which(groupN!=group_val),],weight[which(groupN!=group_val)]) 
  
  psi_S_1 <- (1-p_S_1)/(1-p_S_X_1)
  psi_S_0 <- (1-p_S_0)/(1-p_S_X_0)
  
  # Add wage structure weight to main data.frame
  psi_S <- rep(1,nrow(mf))
  select <- which(groupN==group_val)
  psi_S[select] <- psi_S_1
  select <- which(groupN!=group_val)
  psi_S[select] <- psi_S_0
  

  # set text bar
  setTxtProgressBar(pb,1)
  }
  
  cat("\n\n")
  
  ###########################################################
  ###########################################################
  ### Compute decomposition terms
  
  # compute decomposition if stats==TRUE
  if(stats==TRUE){

  ##########################################################
  ### Observed distributions
  F1 <- stat(dep,weight,groupN,group=1,rwfactor=rep(1,length(weight)),
             tau=tau)
  F0 <- stat(dep,weight,groupN,group=0,rwfactor=rep(1,length(weight)), 
             tau=tau)

  
  ##########################################################
  # Counterfactual distribution(s)
  #if reference group==0 all rw factors are the inverse of the computed
  pow <- ifelse(reference==1,1,-1) 
  first <- ifelse(reference==1,1,nvar)
  last <- ifelse(reference==1,nvar,1)
  loopvals <- seq(first,last)
  
  FC <- NULL
  for(i in 1:nvar){
    FC <- cbind(FC,stat(dep,weight,groupN,group=reference,rwfactor=psi[,i]^pow,
                        tau=tau))
  }
  
  if(nvar==1){
    FC <- as.matrix(FC)
  }

  ##########################################################
  # Decomposition of aggregate effects and detailed composition effects
  Delta <- cbind(F1 - F0, F1 - FC[,nvar], FC[,nvar] - F0)
  if(reference==1){
  colnames(Delta) <- c("Delta_O","Delta_X","Delta_S")
  } else {
  colnames(Delta) <- c("Delta_O","Delta_S","Delta_X")
  }
  
  if(nvar>1){
  if(reference==1){
      Delta <- cbind(Delta,F1-FC[,1])
      colnames(Delta)[length(colnames(Delta))] <- paste("Delta_X",1,sep="")
      
      for(i in 2:nvar){
      Delta <- cbind(Delta,FC[,i-1]-FC[,i])
      colnames(Delta)[length(colnames(Delta))] <- paste("Delta_X",i,sep="")

      }
  } else {
      for(i in nvar:2){
      Delta <- cbind(Delta,FC[,i]-FC[,i-1])
      colnames(Delta)[length(colnames(Delta))] <- paste("Delta_X",i,sep="")
      
      }
      Delta <- cbind(Delta,FC[,1]-F0) 
      colnames(Delta)[length(colnames(Delta))] <- paste("Delta_X",1,sep="")
    }
  }
 
  ###########################################################
  ### Decomposition of wage structure effect
  if(!is.null(ws_formula)){
    
  ws_var <- as.character(fws, lhs=1)[2]
  ws_var_val <- 1
  if(is.factor(mf[,ws_var])){
    ws_var_val <- levels(mf[,ws_var])[2]
  }
  
  # Select only observations with ws group 0
  select <- which(mf[,ws_var]!=ws_var_val)
  
  # Compute counterfactual values  
  FCW1 <- stat(dep[select],weight[select],groupN[select],group=1,rwfactor=psi_S[select],
               tau=tau)
  FCW0 <- stat(dep[select],weight[select],groupN[select],group=0,rwfactor=psi_S[select], 
               tau=tau) 
  
  Delta_WS_X1 <- (F1-FCW1) - (F0-FCW0)
  Delta_WS_other <- Delta[,ifelse(reference==1,3,2)] - Delta_WS_X1
  
  Delta <- cbind(Delta,Delta_WS_X1,Delta_WS_other)
  
  }
  
  ##########################################################
  # Prepare results of decomposition for export
  quantile=cbind(tau,Delta[1:length(tau),]) 
  other.stats=Delta[(length(tau)+1):nrow(Delta),]
  
  } else {
  ##########################################################
  #if no stats return empty objects
  quantile=NULL
  other.stats=NULL
  }
 
  ##########################################################
  ### Export results
  
  res <- list(quantile=quantile,
              other.stats=other.stats,
              formula=formula, 
              mf=mf, 
              weight=weight,
              psi=psi, 
              reference=reference, 
              tau=tau,
              firstrw=firstrw)
  
  if(!is.null(ws_formula)){
    # Add WS weights to the weight matrix for export
    psi <- cbind(psi, psi_S)
    res <- list(quantile=quantile,
                other.stats=other.stats,
                formula=formula, 
                mf=mf, 
                weight=weight,
                psi=psi, 
                reference=reference,
                reference_print=reference_print, 
                tau=tau,
                firstrw=firstrw,
                formula.rw=fws)
  }
  
  return(res)   
}


#############################################################
## Function for fitting and predicting conditional Probabilities

pfit <- function(mod,df,w, newdata=NULL){
  
  # Without survey package
  #dep <- model.frame(mod,df)[,1]
  #reg <- model.matrix(mod,df)
  
  #logit <- glm(dep~reg, weights=w, family = binomial(link = "logit"), na.action=na.exclude, y=FALSE, model=FALSE)
  
  # With survey package
  design <- survey::svydesign(~0, data=df, weights=~w)
  logit <- survey::svyglm(mod, data=df, design=design,family=quasibinomial(link="logit"))
  
  # Predict
  p_X_1  <- predict.glm(logit, newdata=newdata, type="response", na.action = na.exclude)
  
  # Truncate weights 
  #p_X_1[which(p_X_1 < 0.01)] <- 0.01
  #p_X_1[which(p_X_1 > 0.99)] <- 0.99
  
  return(p_X_1)
}

#############################################################
### Gini function (code by Rothe(2015))
Gini <- function (x, w) {
  n <- length(x)
  w <- w/sum(w)
  G <- sum(x[order(x)] * 1:n * w[order(x)])
  G <- 2 * G/(n*sum(x[order(x)]  * w[order(x)]))
  G - 1 - (1/n)
}

#############################################################
### Function for distributional statistics
stat <- function(dep,weight,groupN,group=c(0,1),rwfactor,
                 tau=c(10,50,90,99)/100){
  
  # Factor variables for group selection allowed
  if(is.factor(groupN)){
    group <- levels(groupN)[1+group]
  }
  
  # Select variables
  dep <- dep[which(groupN==group)]
  weight <- weight[which(groupN==group)]
  rwfactor <- rwfactor[which(groupN==group)]
  
  w <- weight*rwfactor
  
  ### Are the necessary quantiles estimated?
  tau1 <- union(c(0.1,0.5,0.9,0.95,0.99),tau)
  tau1 <- tau1[order(tau1)]
  
  # get quantiles statistics
  quantile <- Hmisc::wtd.quantile(dep,weight=w,probs=tau1)
  
  #Get mean and var
  mu <- Hmisc::wtd.mean(dep, weight=w)
  va <-  Hmisc::wtd.var(log(dep), weight=w)
  
  # Estimate additional stats if all stats required
  #Overall gini and income share of top 5%
    gini <- Gini(dep, w) 
    select <- which(dep>=quantile[match(0.95,tau1)])
    s_top05 <- (wtd.mean(dep[select], weight=w[select])/mu)*0.05
    
    #Decile ratios
    quantile1 <- log(quantile)
    p90p10 <- quantile1[match(0.9,tau1)]-quantile1[match(0.1,tau1)]
    p90p50 <- quantile1[match(0.9,tau1)]-quantile1[match(0.5,tau1)]
    p50p10 <- quantile1[match(0.5,tau1)]-quantile1[match(0.1,tau1)]
    p99p90 <- quantile1[match(0.99,tau1)]-quantile1[match(0.9,tau1)]
    
    res <- c(quantile[match(tau,tau1)],mu,va,gini,
             p90p10,p90p50,p50p10,p99p90,
             s_top05)
    names(res)[(length(tau)+1):length(res)] <- c("Mean","Var(log y)","Gini",
                                                 "p90-p10","p90-p50","p50-p10","p99-p90",
                                                 "Top 5% share")
    return(res)
}

#############################################################
### Function for kernel density estimates

kden <- function(dep,weight=NULL,
                 groupN=NULL,group=c(0,1),
                 rwfactor=NULL,
                 px=NULL,
                 bw = "nrd0",
                 kernel="gaussian",
                 n=512,
                 na.rm = TRUE){
  
  # Factor variables for group selection allowed
  if(is.factor(groupN)){
    group <- levels(groupN)[1+group]
  }
  if(is.null(groupN)){group <- "all"}
  
  # Prepare weights
  if(is.null(weight)){weight <- rep(1,length(dep))} 
  if(is.null(rwfactor)){rwfactor <- rep(1,length(dep))}      
  
  # Select variables
  if(is.null(groupN)==FALSE){
    select <- which(groupN==group)
    dep <- dep[select]
    weight <- weight[select]
    rwfactor <- rwfactor[select]
  }
  
  #Adjust weights
  w <- weight*rwfactor
  wsum <- sum(w)
  w <- w/wsum
  if(is.null(px)){px <- 1} else {px <- wsum/px}
  
  if(sum(is.na(dep))!=0&na.rm==TRUE){
    rm <- which(is.na(dep))
    dep <- dep[-rm]
    w <- w[-rm]
  }
  
  #Estimate density
  d <- density(dep,weights=w, 
               kernel=kernel, bw=bw,n=n)
  
  #Return results
  d <- data.frame(group=rep(group,n),x=d$x,density=d$y*px)
  return(d)
}

#############################################################
### Plot function for composition effect results

dfl_deco_plot <- function(result,type=c(1,2,3)){
  
  result <- result[["quantile"]]
  
  if(type==1|ncol(result)==4){
    ## type 1: Observed difference and main decomposition terms (S,X)
    diff <- as.data.frame(result[,c(1:4)])
    
  } else if(type==2){
    ## type 2: All individual terms besides observed difference
    diff <- as.data.frame(result[,-2])
    
  } else{
    ## type 3: Only detailed terms
    diff <- as.data.frame(result[,-c(2:4)])  
  } 
  diff <- melt(diff, id.vars="tau", measure.vars = names(diff)[-1], variable.name= "effect", value.name="delta")
  plot <- ggplot(diff, aes(tau,delta, colour = effect)) + geom_hline(yintercept = 0, colour="grey")  + geom_line()+geom_point(aes(shape=effect, color=effect)) + scale_shape_manual(values=c(15:20,0:14,15:20,0:14))
  
  return(plot)
}



#############################################################
### dfl_diag():
### Diagnosis tool to compare covariates distribution
### of actual and reweighted distribution

dfl_diag <- function(result, compareRef=FALSE, psi.2=FALSE){

  #model and reference group
  f <- as.Formula(result$formula)
  reference <- result$reference
  
  #data
  mf <- result$mf
  weight <- result$weight

  #weights
  if(psi.2==FALSE){
  psi <- as.matrix(result$psi)  
  } else {
  psi <- as.matrix(result$psi.2) 
  }
  # Select psi 
  if(ncol(psi)==1){
  psi <- psi[,1]
  } else if(colnames(psi)[ncol(psi)]=="psi_S"){
  psi <- psi[,ncol(psi)-1]  
  } else {
  psi <- psi[,ncol(psi)]
  }
  
  # Select observations of reference group
  if(is.factor(mf$groupN)){
    reference <- levels(mf$groupN)[reference + 1]
  }
  
  selectRef <- which(mf$groupN==reference)
  # If cond==FALSE use comparison group
  # for comparison to actual values;
  # else use reference group. 
  if(compareRef==FALSE){
  selectCom <- which(mf$groupN!=reference)
  } else {
  selectCom <- selectRef  
  }
  
  #Prepare df
  mod <- formula(f, collapse=TRUE)
  #mRef <- model.matrix(mod,mf)[selectRef,-1]
  #mCom <- model.matrix(mod,mf)[selectCom,-1]
  mRef <- get_all_vars(mod,mf)[selectRef,]
  mCom <- get_all_vars(mod,mf)[selectCom,]
  f <- as.formula(paste0(names(mRef)[1],"~", paste0(names(mRef)[2:ncol(mRef)], collapse="+")))
  mRef <- model.matrix(f,mRef)[,-1]
  mCom <- model.matrix(f,mCom)[,-1]
  wRef <- weight[selectRef]
  wCom <- weight[selectCom]
  psi <- psi[selectRef]
  
  #Find means, diff in means, var/sd 
  mean_obs <- apply(mCom,2,function(x) wtd.mean(x, weights=wCom))
  mean_rw <- apply(mRef,2,function(x) wtd.mean(x, weights=psi*wRef))
  
  sd_obs <- apply(mCom,2,function(x) wtd.var(x, weights=wCom))
  sd_rw <- apply(mRef,2,function(x) wtd.var(x, weights=psi*wRef))
 
  mean_diff <- mean_obs - mean_rw
  
  sd_diff <- sqrt(sd_ob + sd_rw)
  sd_obs <- sqrt(sd_ob)
  sd_rw <- sqrt(sd_rw)
  
  #Export table
  res <- t(rbind(mean_obs,mean_rw,mean_diff,sd_obs, sd_rw,sd_diff))
  return(res)
}

#############################################################
### dfl_stat():
### Returns descriptive statistics of covariates

dfl_stat <- function(formula,
                     data,
                     weights,
                     group,
                     na.action = na.exclude,
                     reference=1,
                     constant=FALSE){
  
  mf = match.call()
  m = match(c("formula", "data", "weights", "na.action","group"), names(mf), 0)
  mf = mf[c(1, m)]
  mf$drop.unused.levels = TRUE
  
  # Retrieve formula as.Formula and model.frame
  f <- as.Formula(formula)
  
  mf[[1]] = as.name("model.frame")
  mf$formula <- f
  mf = eval.parent(mf)
  mt = attr(mf, "terms")
  
  # Store "orignal" decomposition model again
  f <- as.Formula(formula)
  
  # Extract variables, weights and group identifier
  #reg = get_all_vars(mt, mf)
  weight = model.weights(mf)
  if (!is.null(weight) && !is.numeric(weight)) {
    stop("'weights' must be a numeric vector")
  }
  if (!is.null(weight)) {
    weight = weight
  }
  else {
    weight = rep(1, nrow(mf))
  }
  
  #weight = weight/sum(weight)
  groupN = mf[, ncol(mf)]
  #reg$groupN <- groupN
  mf$groupN <- groupN
  
  # Select observations of reference group
  comparison <- ifelse(reference==1,0,1) 
  if(is.factor(mf$groupN)){
    reference <- levels(mf$groupN)[reference + 1]
    comparison <- levels(mf$groupN)[which(levels(mf$groupN)!=reference)]
  }
  
  selectRef <- which(mf$groupN==reference)
  selectCom <- which(mf$groupN!=reference)
  
  #Prepare df
  mod <- formula(f, collapse=TRUE) #include reference group of cat. variables by +0
  if(constant==FALSE){
    mRef <- as.matrix(model.matrix(mod,mf)[selectRef,-1])
    mCom <- as.matrix(model.matrix(mod,mf)[selectCom,-1])
  }else{
    mRef <- as.matrix(model.matrix(mod,mf)[selectRef,])
    mCom <- as.matrix(model.matrix(mod,mf)[selectCom,])
  }
  wRef <- weight[selectRef]
  wCom <- weight[selectCom]
  
  #Find means, diff in means, var/sd 
  mean_Ref <- apply(mRef,2,function(x) wtd.mean(x, weights=wRef))
  mean_Com <- apply(mCom,2,function(x) wtd.mean(x, weights=wCom))
  
  sd_Ref <- apply(mRef,2,function(x) wtd.var(x, weights=wRef))
  sd_Com <- apply(mCom,2,function(x) wtd.var(x, weights=wCom))
  
  mean_diff <- mean_Ref - mean_Com
  
  sd_diff <- sqrt(sd_Ref + sd_Com)
  sd_Ref <- sqrt(sd_Ref)
  sd_Com <- sqrt(sd_Com)
  
  # Sum of weights
  N <- matrix(c(length(wRef),length(wCom),sum(wRef),sum(wCom)),ncol=2,byrow=TRUE)
  colnames(N) <- c(reference,comparison)
  rownames(N) <- c("Obs.","Sum of weights")
  
  #Export table
  res <- t(rbind(mean_Ref,mean_Com,mean_diff,sd_Ref, sd_Com,sd_diff))
  colnames(res) <-  c(paste(rep("mean",3),c(reference,comparison,"diff"),sep="_"),paste(rep("sd",3),c(reference,comparison,"diff"),sep="_"))

  res <- list(means=res, N=N)
  return(res)
}



