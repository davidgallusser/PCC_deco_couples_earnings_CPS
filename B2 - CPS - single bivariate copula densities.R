##############################################################
## B2 - CPS - single bivariate copula densities 1976-2016
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
setwd(sav_d)

##### Define models for wage imputation
impute <- list(list(sf <- LFP_f ~ NCHLT5_f + NCHILD518_f + asinh(abs(INCIDR_hh)) + asinh(INCWAGE_m) + EDUC_f + AGE_f + I(AGE_f^2) + METRO_f  + REGION_f,
                    wf <- log(WAGE_f) ~ EDUC_f*(EXP_f + I(EXP_f^2)) + METRO_f + REGION_f + RACE_f + HISPAN_f,
                    weights="ASECWT_f"),
               list(sm <- LFP_m ~ NCHLT5_m + NCHILD518_m + asinh(abs(INCIDR_hh)) + asinh(INCWAGE_f) + EDUC_m + AGE_m + I(AGE_m^2) + METRO_m  + REGION_m,
                    wm <- log(WAGE_m) ~ EDUC_m*(EXP_m + I(EXP_m^2)) + METRO_m + REGION_m + RACE_m + HISPAN_m ,
                    weights="ASECWT_m"))


#########################################################
# Copulas between male wages and female hours 
f_c <-  I(WAGE_m+ANNUALHOURS_f) ~ WAGE_m + ANNUALHOURS_f 

tic("Decomposition")
set.seed(123)
VCopObj_MWage_FHours <- VC_deco(f_c,cdata, 
                                weights="ASECWTH_m", group="YEAR_m", 
                                var_types = c("c","d"),
                                family_set = "tll", #c("onepar","twopar")
                                nonpar_method = "quadratic",
                                structure=NA,
                                tree_crit = "rho",
                                mult=1,
                                multiply=1,
                                tau=c(seq(5,95,5),99)/100,
                                deco=FALSE, n=50000, quasi = TRUE, detailed_deco = FALSE,
                                impute=impute,          
                                impute_bs_only=TRUE,
                                bootstrap=FALSE, it=200, cores=4,
                                print_info=TRUE)
toc()

setwd(paste0(sav_d,datafiles_d))
save(VCopObj_MWage_FHours , file=paste0("copula_female_hours_male_wages_",years[1],"-",years[2],"_new.rda"))
setwd(sav_d)
