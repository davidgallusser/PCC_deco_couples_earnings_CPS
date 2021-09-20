##########################################################################
#### 6 - VCdeco - wage imputation with a Heckman sample selection model
##########################################################################

## Wrapper function
impute_wages <- function(selection, outcome, data, weights=NULL, method="ml",
                         group=NULL, dropNA=TRUE, print_info=TRUE, short_info=FALSE){
  
  if(is.null(group)){
    data$`(group)` <- 1
    group <- "(group)"
  }
  
  res <- lapply(split(data,data[,group]), function(x)  fit_and_predict(x,selection, outcome, weights=weights,method=method,dropNA=dropNA))
  data <- do.call("rbind", lapply(res, function(x) x[[1]]))
  #res <- 1
  #df <- fit_and_predict(df,selection, outcome, weights=weights,method=method,dropNA=dropNA)
  
  if(print_info){
    cat("\n")
    for(i in 1:length(res)){
      cat("Group", names(res)[i],"\n")
      if(short_info){
      print(coef(summary(res[[i]][[2]]))[c("sigma","rho"), ]) 
      }else{
      print(summary(res[[i]][[2]]))
      }
      cat("\n")
    }
    
    cat("Number of missing wages: \n")
    for(i in 1:length(res)){
      cat("Group", names(res)[i], res[[i]][[3]],"\n")
    }
    
    cat("\nNumber of observations dropped (still missing wages): \n")
    for(i in 1:length(res)){
      cat("Group", names(res)[i], res[[i]][[4]],"\n")
    }
  }
  return(data)
}

### Fit and predict function
fit_and_predict <-  function(df,
                             selection, outcome, weights=NULL, 
                             method="2step",
                             dropNA=TRUE){
  
  # Prepare variables
  selection_var <- as.character(selection)[2]
  wage_var <- as.character(outcome)[2]
  if(length(grep("log\\(",wage_var))>0){
    log_trans = TRUE
    wage_var <- substr(wage_var, 5, nchar(wage_var)-1)
  }
  
  #Add "weights" if missing
  wdrop <- FALSE
  if(is.null(weights)){
    df$weights <- rep(1,nrow(df))
    weights <- "weights"
    wdrop <- TRUE
  }
  
  ## Fit model
  est <- selection(selection, 
                    outcome, 
                    data=df, weights=df[,weights], method = method)

  # Get list of observation with missing wage variable 
  # and complete observations of variables of outcome model
  selection_Xvars <- names(get_all_vars(selection,df))[-1]
  outcome_Xvars <- names(get_all_vars(outcome,df))[-1]
  miX <- unique (unlist (lapply (df[,outcome_Xvars], function(df) which (is.na (df)))))
  s <- setdiff(which(df[,selection_var]==0),miX)
  
  ##Get and sample residuals
  #if(method=="ml"){
  # By hand if "ml": 
  #allvars <- as.formula(paste0(as.character(selection)[2]," ~ ", as.character(selection)[3],"+",as.character(outcome)[2],"+",as.character(outcome)[3]))
  #dfoutcome <- get_all_vars(allvars,df)
  #sel <- unique (unlist (lapply (dfoutcome, function(x) which (is.na (x)))))
  #if(length(sel)>0){
  #dfoutcome <- dfoutcome[-sel,]
  #}
  
  #nselcoef <- length(est$twoStep$probit$estimate)
  #ntotcoef <- length(est$estimate)
  #outcomecoefs <- est$estimate[-c((1:nselcoef),ntotcoef-1,ntotcoef)]
  #lambda <- est$estimate["sigma"] * est$estimate["rho"]
  
  #sel <- which(dfoutcome[,selection_var]==1|dfoutcome[,selection_var]==TRUE)
  #predlink <- model.matrix(selection,dfoutcome[sel,])%*%est$estimate[1:nselcoef]
  #predoutcome <- model.matrix(outcome,dfoutcome[sel,])%*%outcomecoefs + lambda* dnorm(predlink) / pnorm(predlink)  
  #res <- as.numeric(model.frame(outcome, dfoutcome[sel,])[,1]-predoutcome)
  #}else{
  # Using resid() if 2step estimation
  #res <- resid(est)
  #sel <- which(is.na(res)==FALSE)
  #res <- res[sel]
  #}
 
  # Sample second stage residuals (trim sample to avoid extreme residuals)
  #res <- res[which(res>quantile(res,0.01)&res<quantile(res,1-0.01))]
  #res <- sample(res,length(s),replace=TRUE) 
  
  ## Impute wages
  # Select independent variables of outcome variable for censored observations
  outcome2 <- as.formula(paste0(as.character(selection)[2]," ~ ", as.character(outcome)[3]))
  mm <- model.matrix(outcome2, df[s,])
  
  # Make predictions for censored observationss
  # i.e. conditional predicition if E(Y|X,Z,w=1)
  if(method=="ml"){
     selmm <- model.matrix(selection,df[s,])
     nselcoef <- length(est$twoStep$probit$estimate)
     selcoefs <- est$estimate[1:nselcoef]
     matchedselcoefs <- match(names(selcoefs),colnames(selmm))
     predlink <- selmm[,matchedselcoefs]%*%est$estimate[1:nselcoef]
     ntotcoef <- length(est$estimate)
     outcomecoefs <- est$estimate[-c((1:nselcoef),ntotcoef-1,ntotcoef)]
     matchedcoefs <- match(names(outcomecoefs),colnames(mm))
     lambda <- est$estimate["sigma"] * est$estimate["rho"]
  }else{
    predlink <- predict(est, type="link", part="selection")[s]
    outcomecoefs <- est$lm$coefficients[-length(est$lm$coefficients)] 
    matchedcoefs <- match(names(outcomecoefs),paste0("XO",colnames(mm)))
    lambda <- est$rho * est$sigma  
  }  
  
  # Select dependent variable and impute
  dep <- df[,wage_var]
  dep[s] <- mm[,matchedcoefs]%*%outcomecoefs + lambda * dnorm(predlink) / pnorm(predlink)  
  
  # Add error term (normal distribution with variance \hat{sigma}^2)
  res <- rnorm(length(s),0,est$estimate["sigma"]^2)
  dep[s] <- dep[s] + res
  
  # Transform variable
  if(log_trans==TRUE){
    dep[s] <- exp(dep[s])
  }
  

  # Add to original data.frame
  df[,wage_var] <- dep 
 
  #Number of originally missing values
  nimputed <- length(s)
  nomissing <- length(which(df[,selection_var]==0))
  
  # Drop observation with still missing values
  ndropped <- NULL
  if(dropNA==TRUE){
    s <- which(is.na(df[,wage_var])|df[,wage_var]==0)
    ndropped <- length(s)
    if(ndropped>0){
    df <- df[-s,]
    }
  }
  
  if(wdrop==TRUE){
  df$weights <- NULL  
  }  

  return(list(data=df, est=est, nomissing=nomissing, ndropped=ndropped, nimputed=nimputed))
}



#### The following functions imputes wages for individuals who do not 
#### have paid work

# impute_wages <- function(selection, outcome, data, weights=NULL, group,
#                         dropNA=TRUE, print_info=TRUE){
#  
#  selection_var <- as.character(selection)[2]
#  wage_var <- as.character(outcome)[2]
#  if(length(grep("log\\(",wage_var))>0){
#  log_trans = TRUE
#  wage_var <- substr(wage_var, 5, nchar(wage_var)-1)
#  }
  
  #Add "weights" if missing
#  if(is.null(weights)){
#    weights <- rep(1,nrow(data))
#  }
  
#  # Set reference group 
#  if(is.factor(data[,group])){
#    sel0 <- levels(data[,group])[1]
#  }else{
#    sel0 <- sort(unique(data[,group]))[1]
#  }
  
  # Select 
#  data0 <- data[which(data[,group]==sel0),]
#  data1 <- data[which(data[,group]!=sel0),]
  
  #Add "weights" if missing
#  weights0 <- weights[which(data[,group]==sel0)]
#  weights1 <- weights[which(data[,group]!=sel0)]
  
  #Retrieve outcome variable
#  dep0 = data0[,wage_var]
#  dep1 = data1[,wage_var]
  
  ### Run model 
#  est0 <- selection(selection, 
#                     outcome, 
#                     data=data0, weights =weights0, method = "2step")
  
#  est1  <- selection(selection, 
#                            outcome, 
#                            data=data1, weights =weights1, method = "2step")

#  if(print_info==TRUE){
#    print(summary(est0)) 
#    print(summary(est1)) 
#  }
  
  # Select censored observations where we observe
  # all exploratory variables of outcome equation
#  selection_Xvars <- names(get_all_vars(selection,data0))[-1]
#  outcome_Xvars <- names(get_all_vars(outcome,data0))[-1]
#  miX0 <- unique (unlist (lapply (data0[,outcome_Xvars], function (x) which (is.na (x)))))
#  miX1 <- unique (unlist (lapply (data1[,outcome_Xvars], function (x) which (is.na (x)))))
  
#  s0 <- setdiff(which(data0[,selection_var]==0),miX0)
#  s1 <- setdiff(which(data1[,selection_var]==0),miX1)
  
  # Get second-stage residuals 
#  res0 <- resid(est0)
#  res1 <- resid(est1)
#  res0 <- res0[which(is.na(res0)==FALSE)]
#  res1 <- res1[which(is.na(res1)==FALSE)]

  # Number of censored observations
  #N00 <- summary(est0)$param$N0
  #N01 <- summary(est1)$param$N0
  
  # Sample second-stage residuals
#  res0 <- sample(res0,length(s0),replace=TRUE)
#  res1 <- sample(res1,length(s1),replace=TRUE)
  
  ### Impute missing wages
  # Define outcome equation with selection variable as dependent variable (for predictions)
#  outcome2 <- as.formula(paste0(as.character(selection)[2]," ~ ", as.character(outcome)[3]))

  # Conditional predicition if E(Y|X,Z,w=1)
#  predlink0 <- predict(est0, type="link", part="selection")[s0]
#  dep0[s0] <-  model.matrix(outcome2, data0[s0,])%*%est0$lm$coefficients[-length(est0$lm$coefficients)] + est0$rho* est0$sigma * dnorm(predlink0) / pnorm(predlink0)  
  #predict(est0, newdata=data0[s0,], type="unconditional", part="outcome")
  
  # Conditional predicition if E(Y|X,Z,w=1)
#  predlink1 <- predict(est1, type="link", part="selection")[s1]
#  dep1[s1] <- model.matrix(outcome2, data1[s1,])%*%est1$lm$coefficients[-length(est1$lm$coefficients)] + est1$rho* est1$sigma * dnorm(predlink1) / pnorm(predlink1)
  # predict(est0, newdata=data0[s0,], type="unconditional", part="outcome")
  # Conditional predicition if E(Y|X,Z,w=0)
  # predict(est1, newdata=data1[s1,], type="unconditional", part="outcome") - est1$rho* est1$sigma * dnorm(-predlink1) / pnorm(-predlink1)
  
#  dep0[s0] <- res0 + dep0[s0]
#  dep1[s1] <- res1 + dep1[s1]

#  if(log_trans==TRUE){
#    dep0[s0] <- exp(dep0[s0])
#    dep1[s1] <- exp(dep1[s1])
#  }

  # Add to original data.frame
#  data0[,wage_var] <- dep0 
#  data1[,wage_var] <- dep1 
  

  # Drop 
#  if(dropNA==TRUE){
    
    # Show characteristics of observations where imputation failed
    #print(data0[which(is.na(data0[,wage_var])),c(selection_var,wage_var,outcome_Xvars, selection_Xvars)])
    #print(data1[which(is.na(data1[,wage_var])),c(selection_var,wage_var,outcome_Xvars, selection_Xvars)])
    
#    s0 <- which(is.na(data0[,wage_var]))
#    data0 <- data0[-s0,]
    
#    s1 <- which(is.na(data1[,wage_var]))
#    data1 <- data1[-s1,]
    
#    if(print_info==TRUE){cat("Due to missing independent variables:",length(s0), "observations dropped in group 0,", length(s1), "observations dropped in group 1.")}
    
    
#  }
  
#  return(data=rbind(data0,data1))
  
#}

