##############################################################
## C1 CPS decomposition 1976-2016
##############################################################

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

##############################################################
### Perform decomposition

##### Set formulas
f <- EARNINGS ~ WAGE_m:ANNUALHOURS_m + WAGE_f:ANNUALHOURS_f

##### Define models for wage imputation
impute <- list(list(sf <- LFP_f ~ NCHLT5_f + NCHILD518_f + asinh(abs(INCIDR_hh)) + asinh(INCWAGE_m) + EDUC_f + AGE_f + I(AGE_f^2) + METRO_f  + REGION_f,
                    wf <- log(WAGE_f) ~ EDUC_f*(EXP_f + I(EXP_f^2)) + METRO_f + REGION_f + RACE_f + HISPAN_f,
                    weights="ASECWT_f"),
               list(sm <- LFP_m ~ NCHLT5_m + NCHILD518_m + asinh(abs(INCIDR_hh)) + asinh(INCWAGE_f) + EDUC_m + AGE_m + I(AGE_m^2) + METRO_m  + REGION_m,
                    wm <- log(WAGE_m) ~ EDUC_m*(EXP_m + I(EXP_m^2)) + METRO_m + REGION_m + RACE_m + HISPAN_m ,
                    weights="ASECWT_m"))

# Prefered model
structure0=as_rvine_structure(matrix(c(3,2,4,1,2,4,3,0,2,4,0,0,2,0,0,0),4,4))
structure1=as_rvine_structure(matrix(c(2,3,1,4,3,2,1,0,3,2,0,0,3,0,0,0),4,4))

structure0 <- dvine_structure(c(2,1,3,4))# get_structure(VCopObj$VCfit$VCfit1)#rvine_structure(c(2,1,3,4))
structure1 <- dvine_structure(c(2,1,3,4))# get_structure(VCopObj$VCfit$VCfit1)#rvine_structure(c(2,1,3,4))

########################################################################
### Perfom main decomposition
cat("Main decomposition \n")
tic("Decomposition")
set.seed(123)
VCopObj<- VC_deco(f,cdata, 
                  weights="ASECWTH_m", group="YEAR_m", 
                  var_types = c("c","d","c","d"),
                  family_set = "tll", #c("onepar","twopar")
                  nonpar_method = "quadratic",
                  structure=list(structure0,
                                 structure1),
                  tree_crit = "rho",
                  mult=1,
                  multiply=1,
                  tau=c(seq(5,95,5),99)/100,
                  deco=TRUE, n=50000, quasi = TRUE, 
                  detailed_deco = TRUE,
                  impute=impute,          
                  impute_bs_only=TRUE,
                  bootstrap=TRUE, it=200, cores=4,
                  print_info=TRUE)
toc()

setwd(paste0(sav_d,datafiles_d))
save(VCopObj, file=paste0("copula_decomposition_main_",years[1],"-",years[2],".rda"))
setwd(sav_d)

########################################################################
### Perform decomposition with alternative reference group
cat("Decomposition: Alternative reference group \n")
tic("Decomposition: Alternative reference group")
set.seed(123)
VCopObj_ref0 <- VC_deco(f,cdata, 
                  weights="ASECWTH_m", group="YEAR_m", 
                  var_types = c("c","d","c","d"), 
                  family_set = "tll", #c("onepar","twopar")
                  nonpar_method = "quadratic",
                  structure=list(structure0,
                                 structure1),
                  tree_crit = "rho",
                  reference=0,
                  mult=1,
                  multiply=1,
                  tau=c(seq(5,95,5),99)/100,
                  deco=TRUE, n=50000, quasi = TRUE,
                  detailed_deco = TRUE,
                  impute=impute,          
                  impute_bs_only=TRUE,
                  bootstrap=TRUE, it=200, cores=4,
                  print_info=TRUE)
toc()

setwd(paste0(sav_d,datafiles_d))
save(VCopObj_ref0, file=paste0("copula_decomposition_ref0_",years[1],"-",years[2],".rda"))
setwd(sav_d)


########################################################################
### Perform decomposition with alternative copula structure/order

structure0alt <- dvine_structure(c(1,2,3,4))
structure1alt <- dvine_structure(c(1,2,3,4))

cat("Decomposition: Alternative order \n")
set.seed(123)
VCopObj_alt_order <- VC_deco(f,cdata, 
                        weights="ASECWTH_m", group="YEAR_m", 
                        var_types = c("c","d","c","d"), 
                        family_set = "tll", #c("onepar","twopar")
                        nonpar_method = "quadratic",
                        structure=list(structure0alt,
                                       structure1alt),
                        tree_crit = "rho",
                        reference=1,
                        mult=1,
                        multiply=1,
                        tau=c(seq(5,95,5),99)/100,
                        deco=TRUE, n=50000, quasi = TRUE,
                        detailed_deco = TRUE,
                        impute=impute,          
                        impute_bs_only=TRUE,
                        bootstrap=TRUE, it=200, cores=4,
                        print_info=TRUE)
toc()

setwd(paste0(sav_d,datafiles_d))
save(VCopObj_alt_order, file=paste0("copula_decomposition_alt_order_",years[1],"-",years[2],".rda"))
setwd(sav_d)

########################################################################
### Perform decomposition with parametric copulas

structure0 <- dvine_structure(c(2,1,3,4))# get_structure(VCopObj$VCfit$VCfit1)#rvine_structure(c(2,1,3,4))
structure1 <- dvine_structure(c(2,1,3,4))# get_structure(VCopObj$VCfit$VCfit1)#rvine_structure(c(2,1,3,4))

cat("Decomposition parametric PPC \n")
tic("Decomposition")
set.seed(123)
VCopObj_parest <- VC_deco(f,cdata, 
                  weights="ASECWTH_m", group="YEAR_m", 
                  var_types = c("c","d","c","d"),
                  family_set = "parametric", #"tll"c("onepar","twopar")
                  nonpar_method = "quadratic",
                  par_method = "mle",
                  structure=list(structure0,
                                 structure1),
                  tree_crit = "rho",
                  mult=1,
                  multiply=1,
                  tau=c(seq(5,95,5),99)/100,
                  deco=TRUE, n=50000, quasi = TRUE,
                  detailed_deco = TRUE,
                  impute=impute,          
                  impute_bs_only=TRUE,
                  bootstrap=TRUE, it=200, cores=4,
                  print_info=TRUE)
toc()

setwd(paste0(sav_d,datafiles_d))
save(VCopObj_parest, file=paste0("copula_decomposition_parametric_",years[1],"-",years[2],".rda"))
setwd(sav_d)

########################################################################
### Perform decomposition conditioning on educational composition 

structure0 <- dvine_structure(c(2,1,3,4))# get_structure(VCopObj$VCfit$VCfit1)#rvine_structure(c(2,1,3,4))
structure1 <- dvine_structure(c(2,1,3,4))# get_structure(VCopObj$VCfit$VCfit1)#rvine_structure(c(2,1,3,4))


# DFL decomposition
# deco <- dfl_deco(log(EARNINGS) ~ EDUC_f * EDUC_m + EDUC_m*EXP2_m + EDUC_f*EXP2_f + EXP2_f * EXP2_m ,cdata,weights=ASECWTH_m,group = YEAR_m,
#                  tau=c(1:9/10,0.95,0.99))
# dfl_deco_plot(deco)

# Formula reweighting on educational 
f_c <-  EARNINGS ~ WAGE_m:ANNUALHOURS_m + WAGE_f:ANNUALHOURS_f | EDUC2_f * EDUC2_m + EDUC2_m*EXP2_m + EDUC2_f*EXP2_f + EXP2_f * EXP2_m

# Add professional school and PhD graduates to Master's degree graduates
# in 2016 since they are not distinguished in 1976:
cdata$EDUC2_f <- cdata$EDUC_f
cdata$EDUC2_m <- cdata$EDUC_m
cdata[which(cdata$YEAR_m==2016&cdata$EDUC2_m%in%c("7","8")),"EDUC2_m"] <- "6"
cdata[which(cdata$YEAR_f==2016&cdata$EDUC2_f%in%c("7","8")),"EDUC2_f"] <- "6"

cat("Decomposition conditional on education \n")
tic("Decomposition conditional on education")
set.seed(123)
VCopObj_cond <- VC_deco(f_c,cdata, 
                  weights="ASECWTH_m", group="YEAR_m", 
                  var_types = c("c","d","c","d"),
                  family_set = "tll", #c("onepar","twopar")
                  nonpar_method = "quadratic",
                  structure=list(structure0,
                                 structure1),
                  tree_crit = "rho",
                  mult=1,
                  multiply=1,
                  tau=c(seq(5,95,5),99)/100,
                  deco=TRUE, n=50000, quasi = TRUE, 
                  detailed_deco = TRUE,
                  impute=impute,          
                  impute_bs_only=TRUE,
                  bootstrap=TRUE, it=200, cores=4,
                  print_info=TRUE)
toc()

setwd(paste0(sav_d,datafiles_d))
save(VCopObj_cond, file=paste0("copula_decomposition_cond_on_edu_",years[1],"-",years[2],".rda"))
setwd(sav_d)

#################################################################
#### Inspect fits
# setwd(paste0(sav_d,datafiles_d))
# #load(paste0("copula_decomposition_",years[1],"-",years[2],".rda"))
# load(paste0("copula_decomposition_main_",years[1],"-",years[2],".rda"))
# setwd(sav_d)
# 
# VC_deco_plot_obs_vs_fit(VCopObj)
# # Check overall fit
# VC_deco_qqplot(VCopObj,formula=NULL, type=1, logscale=TRUE, group=c("1975","2015"))
# VC_deco_qqplot(VCopObj,formula=NULL, type=2, logscale=FALSE)
# VC_deco_qqplot(VCopObj,formula=NULL, type=3, logscale=TRUE, tau=3:99/100)
# VC_deco_qqplot(VCopObj,formula=NULL, type=3, logscale=FALSE, tau=3:99/100)
# 
# # Check bivariate fits
# f <- EARNINGS ~ ANNUALHOURS_f
# VC_deco_qqplot(VCopObj,formula=f, type=2, logscale=FALSE, tau=1:99/100)
# f <- EARNINGS ~ ANNUALHOURS_m
# VC_deco_qqplot(VCopObj,formula=f, type=2, logscale=FALSE, tau=1:99/100)
# f <- EARNINGS ~ I(ANNUALHOURS_f>0)*I(ANNUALHOURS_m>0)
# VC_deco_qqplot(VCopObj,formula=f, type=3, logscale=FALSE, tau=1:99/100)
# f <- EARNINGS ~ I(ANNUALHOURS_f+ANNUALHOURS_m)
# VC_deco_qqplot(VCopObj,formula=f, type=3, logscale=FALSE, tau=1:99/100)
# f <- EARNINGS ~ I((WAGE_f+WAGE_m)/2)
# VC_deco_qqplot(VCopObj,formula=f, type=3, logscale=TRUE, tau=2:95/100)
# f <- EARNINGS ~ I(ANNUALHOURS_m*WAGE_m)
# VC_deco_qqplot(VCopObj,formula=f, type=3, logscale=TRUE, tau=2:95/100)
# f <- EARNINGS ~ I(ANNUALHOURS_f*WAGE_f)
# VC_deco_qqplot(VCopObj,formula=f, type=3, logscale=TRUE, tau=2:95/100)
# 
# # Plot structures
# 
# plot(VCopObj$VCfit$VCfit0, tree=c(1,2,3), edge_labels = "tau")
# plot(VCopObj$VCfit$VCfit1, tree=c(1,2,3), edge_labels = "tau")
# 
# 
# # Check densities
# plotlist <- get_density_plots(VCopObj,varnames=c("Wm","Hm","Wf","Hf"))
# compare_plots(plotlist,edge=1, groups=c(1976,2016))
# compare_plots(plotlist,edge=2, groups=c(1976,2016))
# compare_plots(plotlist,edge=3, groups=c(1976,2016))
# compare_plots(plotlist,edge=4, groups=c(1976,2016))
# compare_plots(plotlist,edge=5, groups=c(1976,2016))
# compare_plots(plotlist,edge=6, groups=c(1976,2016))
# 
# ### Perform partial correlation test to test for simplifying assumption
# # vars <- c(4,3,2)
# # e1 <- bicop(VCopObj$VCfit$u1[,c(vars[1],vars[3])],family_set="tll", weights=VCopObj$VCfit$data1[,"(weights)"])
# # e2 <- bicop(VCopObj$VCfit$u1[,c(vars[2],vars[3])],family_set="tll", weights=VCopObj$VCfit$data1[,"(weights)"])
# # zr1 <- hbicop(VCopObj$VCfit$u1[,c(vars[1],vars[3])],e1,cond_var=2)
# # zr2 <- hbicop(VCopObj$VCfit$u1[,c(vars[2],vars[3])],e1,cond_var=2)
# # e21 <- bicop(cbind(zr1,zr2),family_set="tll", weights=VCopObj$VCfit$data1[,"(weights)"])
# # plot(e1)
# # plot(e2)
# # plot(e21)
# # pacotest(cbind(zr1,zr2),VCopObj$VCfit$u1[,vars[3]],'CCC')
# # 
# # vars <- c(2,3,1)#c(2,1,3)
# # e1 <- bicop(VCopObj$VCfit$u1[,c(vars[1],vars[3])],family_set="tll", weights=VCopObj$VCfit$data1[,"(weights)"])
# # e2 <- bicop(VCopObj$VCfit$u1[,c(vars[2],vars[3])],family_set="tll", weights=VCopObj$VCfit$data1[,"(weights)"])
# # zr3 <- hbicop(VCopObj$VCfit$u1[,c(vars[1],vars[3])],e1,cond_var=2)
# # zr4 <- hbicop(VCopObj$VCfit$u1[,c(vars[2],vars[3])],e2,cond_var=2)
# # e22 <- bicop(cbind(zr3,zr4),family_set="tll", weights=VCopObj$VCfit$data1[,"(weights)"])
# # plot(e1)
# # plot(e2)
# # plot(e22)
# # pacotest(cbind(zr3,zr4),VCopObj$VCfit$u1[,vars[3]],'CCC')
# # 
# # zr5 <- hbicop(cbind(zr1,zr2),e21,cond_var=2)
# # zr6 <- hbicop(cbind(zr3,zr4),e22,cond_var=1)
# # e31 <- bicop(cbind(zr5,zr6),family_set="tll", weights=VCopObj$VCfit$data1[,"(weights)"])
# # plot(e31)
# # pacotest(cbind(zr5,zr6),zr2,'CCC')
# 
# # Inspect density estimates
# density_plots  <- VC_deco_plot_density(VCopObj, varnames=c("Wm","Hm","Wf","Hf"))
# VC_deco_plot_compare(density_plots, conditional=0)
# VC_deco_plot_compare(density_plots, conditional=1)
# 
# VC_deco_plot_compare(density_plots, edge=1, diff=TRUE)
# VC_deco_plot_compare(density_plots, edge=2, diff=TRUE)
# VC_deco_plot_compare(density_plots, edge=3, diff=TRUE)
# VC_deco_plot_compare(density_plots, edge=4, diff=TRUE)
# VC_deco_plot_compare(density_plots, edge=5, diff=TRUE)
# VC_deco_plot_compare(density_plots, edge=6, diff=TRUE)
# 
# # Evaluate fit:
# # Plot differences between observed and simulated quantiles
# VC_deco_plot_obs_vs_fit(VCopObj)
# VCopObj$VCdeco$aggregate$quants[c(1,7),]
# 

# ##############################################################
# # Evaluate fit:
# # Plot differences between observed and simulated quantiles
# 
# VC_deco_plot_obs_vs_fit(VCopObj)
# VCopObj$VCdeco$aggregate$quants[c(1,7),]
# 
# # Kolmogorov-Smirnov and Cramèr-von-Mises-Statistics
# cat("p-values: observed dist. = fitted dist.",
#     paste0("\n",names(VCopObj$VCdeco_se$aggregate$quants_sim_vs_obs_KS_CMS_test),": ",VCopObj$VCdeco_se$quants_sim_vs_obs_KS_CMS_test))

