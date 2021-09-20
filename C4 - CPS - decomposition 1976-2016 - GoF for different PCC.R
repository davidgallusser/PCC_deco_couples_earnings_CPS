##############################################################
## C4 - CPS - decomposition 1976-2016 - GoF for different PCC
##############################################################

##### Load data
# Setwd
setwd(paste0(wd,datafiles_d))
years <- c(1976,2016)
cdataall <- NULL
for(i in years){
  load(paste0("CPScouples",i,".rda"))
  cdataall <- rbind(cdataall,cdata)
}
cdata <- cdataall
rm("cdataall")


##### Set formulas
f <- EARNINGS ~ WAGE_m:ANNUALHOURS_m + WAGE_f:ANNUALHOURS_f
impute=NULL

structure0 <- NA
structure1 <- NA

#########################################################
#### Find best models
cat("Find best orders\n")
tic("Find best orders")
set.seed(123)
VC_best_models <- VC_deco(f,cdata, 
                  weights="ASECWTH_m", group="YEAR_m", 
                  var_types = c("c","d","c","d"),
                  family_set = "tll", #c("onepar","twopar")
                  nonpar_method = "quadratic",
                  structure=list(structure0,
                                 structure1),
                  tree_crit = "AIC",
                  find_best_tree = TRUE,
                  mult=1,
                  multiply=1,
                  tau=c(seq(5,95,5),99)/100,
                  deco=FALSE, n=50000, quasi = TRUE, 
                  detailed_deco = FALSE,
                  impute=impute,          
                  impute_bs_only=TRUE,
                  bootstrap=FALSE, it=200, cores=4,
                  print_info=TRUE)
toc()

#########################################################
#### 
setwd(paste0(sav_d,datafiles_d))
save(VC_best_models, file=paste0("copula_best_models_main_",years[1],"-",years[2],".rda"))
setwd(sav_d)

#########################################################
#### Test single models
##### Set formulas
f <- EARNINGS ~ WAGE_m:ANNUALHOURS_m + WAGE_f:ANNUALHOURS_f
impute=NULL

setwd(paste0(sav_d,datafiles_d))
load(paste0("copula_best_models_main_",years[1],"-",years[2],".rda"))
setwd(sav_d)

VC_best_models$vine_IC
VC_best_models$cvine[[4]]
structure0 <- VC_best_models$cvine[[4]]
structure1 <- VC_best_models$cvine[[4]]

#########################################################
#### Find best models
cat("Fit single model\n")
tic("Fit single model")
set.seed(123)
VC_single_model <- VC_deco(f,cdata, 
                          weights="ASECWTH_m", group="YEAR_m", 
                          var_types = c("c","d","c","d"),
                          family_set = "tll", #c("onepar","twopar")
                          nonpar_method = "quadratic",
                          structure=list(structure0,
                                         structure1),
                          tree_crit = "rho",
                          find_best_tree = FALSE,
                          mult=1,
                          multiply=1,
                          tau=c(seq(5,95,5),99)/100,
                          deco=TRUE, n=50000, quasi = TRUE, 
                          detailed_deco = FALSE,
                          impute=impute,          
                          impute_bs_only=TRUE,
                          bootstrap=FALSE, it=200, cores=4,
                          print_info=TRUE)
toc()


# Check overall fit
VC_deco_qqplot(VC_single_model,formula=NULL, type=1, logscale=TRUE, group=years)
VC_deco_qqplot(VC_single_model,formula=NULL, type=2, logscale=FALSE)
VC_deco_qqplot(VC_single_model,formula=NULL, type=4, logscale=TRUE, tau=5:99/100)
VC_deco_qqplot(VC_single_model,formula=NULL, type=4, logscale=FALSE, tau=5:99/100)

