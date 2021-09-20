##########################################################################
#### 9 - VCdeco - data prep functions
##########################################################################

# This function merges couples
merge_couples <- function(pdata, hdata=NULL,
                          hid=NULL, pid, parid, coupleid=NULL, 
                          gendervar, level_male=1,
                          pvars=NULL, hhvars=NULL){
  
  # Check if pvars/hhvars null
  if(is.null(pvars)){pvars <- 1:ncol(pdata)}
  if(is.null(hhvars)&is.null(hdata)==FALSE){hhvars <- 1:ncol(hdata)}
  
  # Select variables
  pdata <- pdata[,pvars]
  if(is.null(hdata)==FALSE){hdata <- hdata[,hhvars]}
  
  # Select male/female
  pdata[,gendervar] <- as.factor(pdata[,gendervar])
  xm <- pdata[which(pdata[,gendervar]==levels(pdata[,gendervar])[level_male]),]
  xf <- pdata[which(pdata[,gendervar]!=levels(pdata[,gendervar])[level_male]),]
  
  # Rename variable names with gender suffix
  names(xm) <- paste0(names(xm),"_m")
  names(xf) <- paste0(names(xf),"_f")
  
  # ID var names with gender suffic
  pid_m <-  paste0(pid,"_m")
  parid_m <-  paste0(parid,"_m")
  pid_f <-  paste0(pid,"_f")
  parid_f <-  paste0(parid,"_f")
  
  if(is.null(hdata)==FALSE){
    hid_m <-  paste0(hid,"_m")
    hid_f <-  paste0(hid,"_f")
  }
  
  # Match females to male spouses
  m_f_ds <- match(xf[,pid_f],xm[,parid_m])
  na_f <- which(is.na(m_f_ds))
  
  # Match females to female spouse
  m_f_ss <- match(xf[,pid_f],xf[,parid_f])
  na_f_ss <- which(is.na(m_f_ss))
  n_f_ss <- sum(is.na(m_f_ss)==FALSE)
  
  # Match males to male spouses 
  m_m_ss <- match(xm[,pid_m],xm[,parid_m])
  na_m_ss <- which(is.na(m_m_ss))
  n_m_ss <- sum(is.na(m_m_ss)==FALSE)
  
  ## Create new data.frame 
  df <- cbind(xm[m_f_ds,],xf)
  df <- df[-na_f, ]
  df$samesex <- FALSE
  
  ## Same sex couples: Women
  if(n_f_ss>0){
    df_f_ss <- cbind(xf[m_f_ss,],xf)
    df_f_ss <- df_f_ss[-na_f_ss, ]
    df_f_ss <- df_f_ss[which(duplicated(df_f_ss[,parid_f])==FALSE),]
    df_f_ss$samesex <- TRUE
  }else{
    df_f_ss <- NULL
  }
  
  # Same sex couples: Men
  if(n_m_ss>0){
    df_m_ss <- cbind(xm[m_m_ss,],xm)
    df_m_ss <- df_m_ss[-na_m_ss, ]
    df_m_ss <- df_m_ss[which(duplicated(df_m_ss[,parid_m])==FALSE),]
    df_m_ss$samesex <- TRUE
  }else{
    df_m_ss <- NULL
  }
  
  # Add household variables 
  if(is.null(hdata)==FALSE){
    m_hh <- match(df[,hid_m],hdata[,hid])
    df <- cbind(df,hdata[m_hh,])
    
    if(n_f_ss>0){
    m_hh <- match(df_f_ss[,hid_f],hdata[,hid])
    df_f_ss <- cbind(df_f_ss,hdata[m_hh,])
    }
    if(n_m_ss>0){
    m_hh <- match(df_m_ss[,hid_m],hdata[,hid])
    df_m_ss <- cbind(df_m_ss,hdata[m_hh,])
    }
  }
  
  # Export results
  dfall <- list(different_sex=df,
                same_sex_f=df_f_ss,
                same_sex_m=df_m_ss)
  return(dfall)
}

# Recode NA
recode_na  <- function(x, 
                       f = c(".S (-3)", ".M (-2)", ".C (-1)"),
                       n = -3:-1
                       ){
  if(is.factor(x)){
    s <- which(is.element(x,f)) 
  }else{
    s <- which(is.element(x,n))  
  }
  x[s] <- NA
  return(x)
}

