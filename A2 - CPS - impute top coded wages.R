##############################################################
## A2 Impute top coded wages
##############################################################

# See: https://cps.ipums.org/cps/topcodes_tables.shtml
# See: https://onlinelibrary.wiley.com/doi/pdf/10.1111/ecin.12299 
# See: https://www.bls.gov/opub/mlr/2009/08/art1full.pdf 

## There are 3 earnings variables

# INCWAGE  (Census: 1962-1987: I51A or WSAL_VAL, since 1988: WSAL-VAL): "Wage and salary income" or "Wages and salaries"
# INCLONGJ (Census: ERN-VAL): "Earnings from longest job" or "Primary earnings"
# OINCWAGE (Census: WS-VAL): "Earnings from other work included wage and salary earnings"

# Since 1988: INCWAGE = ifelse(SRCEARN==1,INCLONGJ,0) + OINCWAGE
# Topcoding: INCWAGE (1962-1987), INCLONGJ (since 1988) and OINCWAGE (since 1988)

# Different types of topcoding
  # Topcoding (i.e. max value) 1962-1995
  # Mean replacement above thresholds 1996-2010
  # Rank proximity swapping 2010-2020 (does not affect distribution!)
# However, even internal CPS data is "hard" topcoded 


# We use Armour, Burkhauser and Larrimore (2016) 
# estimated cell means to impute top coded
# wages with a Pareto distribution. For years from 2011
# we use their method to calculate cell means. 

#########################################################
### Import cell means above thresholds by Armout et al. (2016)

cell_means <- read.table("data_public/Armour_et_al-2016-Table-A1.txt", dec=".", header=TRUE)
cell_means1 <- cell_means[,c(3,4)]
names(cell_means1) <- names(cell_means)[1:2]
cell_means <- rbind(cell_means[,c(1,2)],cell_means1)
cell_means <- cell_means[-which(is.na(cell_means$YEAR)),]
cell_means$VARIABLE <- c(rep("INCWAGE", length(1967:1986)),rep("INCLONGJ", length(1987:2013)))
cell_means$YEAR <- cell_means$YEAR+1 

#########################################################
### Topcoding thersholds 1962-1995
tc1 <- data.frame(YEAR=c(1962:1967, 1968:1975, 1976:1981, 1982:1984, 1985:1987,1988:1995),
           INCWAGE=c(rep(50000,6),rep(50000,8), rep(50000,6), rep(75000,3),rep(99999,3),rep(199998,8)),
           INCLONGJ=c(rep(NA,6+8+6+3+3),rep(99999,8)),
           OINCWAGE=c(rep(NA,6+8+6+3+3),rep(99999,8)))

#########################################################
### Thresholds for replacement values 1996-2010
tc2 <- data.frame(YEAR=c(1996:2002, 2003:2010),
           INCLONGJ=c(rep(150000,7),rep(200000,8)),
           OINCWAGE=c(rep(25000,7),rep(35000,8)))

#########################################################
# Add threshold as minimu values for pareto extrapoation
cell_means$xmin <- NA
cell_means[1:20,"xmin"] <- tc1[7:26,"INCWAGE"]
cell_means[21:28,"xmin"] <- tc1[27:nrow(tc1),"INCLONGJ"]
cell_means[29:43,"xmin"] <- tc2[,"INCLONGJ"]

#########################################################
# Add pareto coefficients
cell_means$alpha <-  apply(cbind(cell_means$MEAN_WAGE,cell_means$xmin),1, function(x) cpareto(x[1],x[2])[1])

#########################################################
### Compute Pareto coeffiecient with the method
### of Armour et al (2016) for years after 2011

# Load data
setwd(paste0(wd,datafiles_d))
years <- c(2011:2020)
CPSl <- NULL 
for(i in years){
  print(i)
  load(file=paste0("CPS",i,".rda"))
  CPSl <- rbind(CPSl,CPS)
}
CPS <- CPSl
rm("CPSl")

# Estimate alpha 
alpha <- unlist(lapply(split(CPS,CPS$YEAR), function(x) hillest(x$INCLONGJ,x$ASECWT,taumin=0.99,xmin=1099999)))

# Add alpha to table
cell_means2 <- data.frame(YEAR=years,
                          MEAN_WAGE=(alpha/(alpha-1))*1099999,
                          VARIABLE="INCLONGJ",
                          xmin=1099999,
                          alpha=alpha)
cell_means <- rbind(cell_means[-(44:47),],cell_means2)


#########################################################
# Perform imputation for all datasets 

# Top income are ranomdly imputed
# INCWAGE is imputed directly until 1987
# Since 1988 INCLONGJ is imputed. 
# INCWAGE is constructed as the sum of INCLONGJ and OINCWAGE if SRCEARN==1,
# Note: OINCWAGE is also

years <- cell_means$YEAR
for(i in years){
  cat("\n",i)
  load(file=paste0("CPS",i,".rda"))
  
  #Recode missings
  CPS$INCWAGE  <- recode_na(CPS$INCWAGE, n=c(99999999,99999998))
  CPS$INCLONGJ <- recode_na(CPS$INCLONGJ, n=c(99999999,99999998))
  CPS$OINCWAGE <- recode_na(CPS$OINCWAGE, n=c(99999999,99999998))
  
  #Create "non-imputed" variable
  CPS$INCWAGE_NI <- CPS$INCWAGE
  CPS$INCLONGJ_NI <- CPS$INCLONGJ
  
  #Impute top coded values
  if(i>1987){
  sel1 <- which(CPS$OINCWAGE>0)
  sel2 <- which(CPS$OINCWAGE==max(CPS$OINCWAGE, na.rm=TRUE))
  cat(" - Max OINCWAGE:",max(CPS$OINCWAGE, na.rm=TRUE), "- share topcoded OINCWAGE",
      round(sum(CPS[sel2,"ASECWT"])/sum(CPS[sel1,"ASECWT"],6)))
  }
  
  #Impute top coded values
  set.seed(123)
  CPS <- impute_top_inc(CPS,i,cell_means)
  
  # Save imputed data 
  save(CPS, file=paste0("CPS",i,".rda"))
}

