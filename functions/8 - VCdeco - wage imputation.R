##########################################################################
#### 8 - VCdeco - wage imputation with a Heckman sample selection model
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
    log_t = TRUE
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
  # and complete observations in all included variables (i.e. in outcome & selection model)
  selection_Xvars <- names(get_all_vars(selection,df))[-1]
  outcome_Xvars <- names(get_all_vars(outcome,df))[-1]
  miX <- unique (unlist (lapply (df[,outcome_Xvars], function(df) which(is.na(df)))))
  miS <- unique (unlist (lapply (df[,selection_Xvars], function(df) which(is.na(df)))))
  s <- setdiff(which(df[,selection_var]==0),union(miS,miX))
  ## Impute wages
  # Select independent variables of outcome variable for censored observations
  outcome2 <- as.formula(paste0(as.character(selection)[2]," ~ ", as.character(outcome)[3]))
  mm <- model.matrix(outcome2, df[s,])
  
  # Make predictions for censored observationss
  # i.e. conditional predicition if E(Y|X,Z,w=1)
  if(method=="ml"){
     # Get covariates and coefficients
     selmm <- model.matrix(selection,df[s,])
     nselcoef <- length(est$twoStep$probit$estimate)
     selcoefs <- est$estimate[1:nselcoef]
     matchedselcoefs <- match(names(selcoefs),colnames(selmm))
     ntotcoef <- length(est$estimate)
     outcomecoefs <- est$estimate[-c((1:nselcoef),ntotcoef-1,ntotcoef)]
     matchedcoefs <- match(names(outcomecoefs),colnames(mm))
     
     # Predict probit link and get estimated elements error covariance matrix
     predlink <- selmm[,matchedselcoefs]%*%est$estimate[1:nselcoef]
     sigma <-  est$estimate["sigma"]
     rho <- est$estimate["rho"]
  }else{
    predlink <- predict(est, type="link", part="selection")[s]
    outcomecoefs <- est$lm$coefficients[-length(est$lm$coefficients)] 
    matchedcoefs <- match(names(outcomecoefs),paste0("XO",colnames(mm)))
    sigma <- est$sigma  
    rho <- est$rho 
  }  
  
  # Get parameters of conditional distribution of Y given 'truncated' variable X
  lambda <- dnorm(predlink) / (1 - pnorm(predlink)) 
  delta <- lambda*(lambda-predlink)
  
  # Sqrt of conditional variance of observation i 
  sigma_i <- sqrt((1-rho^2*delta)*sigma^2)
  
  # Select dependent variable 
  dep <- df[,wage_var]
  
  # Impute:
  # y = x'b + rho*sigma*lambda(z'c) + nu with nu ~ N(0,(1-rho^2*delta)*sigma^2) and delta = lambda(z'c)*(lambda(z'c) - z'c)
  dep[s] <- mm[,matchedcoefs]%*%outcomecoefs + rho*sigma*lambda + rnorm(length(s),0,sigma_i) 
  
  # Transform variable
  if(log_t==TRUE){
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


