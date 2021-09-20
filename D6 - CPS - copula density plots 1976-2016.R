##############################################################
## D6 - CPS - copula density plots 1976-2016
##############################################################

years <- c(1976,2016)
setwd(paste0(wd,datafiles_d))
load(paste0("copula_decomposition_main_",years[1],"-",years[2],".rda"))
setwd(wd)

##############################################################
## Plot checkerboard copula
# Base on gaussian copula with \rho=0.3
#        0.5  0.5
# 0.7    0.4, 0.3
# 0.3    0.1, 0.2

u <- expand.grid(X=1:100,Y=1:100)/100
u$d <- NA
u[which(u$X<0.7&u$Y<0.5),"d"] <- 0.4/(0.5*0.7)
u[which(u$X>=0.7&u$Y<0.5),"d"] <- 0.1/(0.5*0.3)
u[which(u$X<0.7&u$Y>=0.5),"d"] <- 0.3/(0.5*0.7)
u[which(u$X>=0.7&u$Y>=0.5),"d"] <- 0.2/(0.5*0.3)
#u <- matrix(u$d, nrow=100, ncol=100, byrow=TRUE)
chequerboard <- ggplot(u, aes(X,Y,fill=d)) + geom_tile() +  scale_fill_gradientn(colours = divcolors[length(divcolors):1]) + labs(fill="density")


# Save
fn <- paste0(pfolder,"chequerboard_copula",ptype)
ggsave(fn, plot=chequerboard, units=punits, scale=pscale, width=pheight*1.3, height=pheight)


# ##############################################################
# ## Plot vine graph
# cop1 <- VCopObj$VCfit$VCfit1
# cop1$names[1:4] <- c("Wage Men","Hours Men", "Wage Women", "Hours Women")
# #plot(cop1, tree=c(1,2,3), var_names = "use", edge_labels = "pair")
# vg <- vinegraph(cop1, tree=c(1,2,3), var_names = "use", edge_labels = NULL)
# plot(vg)
# 
# # Save
# fn <- paste0(pfolder,"vine_graph",ptype)
# ggsave(fn, plot=vg, units=punits, scale=0.7, width=pwidth, height=pheight)
# 


##############################################################
## Plot densities

#varnames=c("Wm","Hm","Wf","Hf")
#
varnames=c("Wage men", "Hours men", "Wage women", "Hours women")
plotlist <- get_density_plots(VCopObj, varnames=varnames, colorpalette=NULL)

# Density 1
density1 <- compare_plots(plotlist,
                          edge=1, groups=c(1975,2015),
                          colorpalette=rev(divcolors),
                          title=TRUE,
                          zlim=5,
                          bias=2.5)
density1
fn <- paste0(pfolder,"density_hours_wages_men",ptype)
ggsave(fn, plot=density1, units=punits, scale=pscale, width=pwidth, height=pheight)

# Assortative mating
density2 <- compare_plots(plotlist,
                          edge=2, groups=c(1975,2015),
                          colorpalette=rev(divcolors),
                          title=FALSE,
                          bias=2)
density2
fn <- paste0(pfolder,"density_mating",ptype)
ggsave(fn, plot=density2, units=punits, scale=pscale, width=pwidth, height=pheight)


# Density 3
density3 <- compare_plots(plotlist,
                          edge=3, groups=c(1975,2015),
                          colorpalette=rev(divcolors),
                          title=FALSE,
                          zlim=2.5,
                          bias=1.25)
density3
fn <- paste0(pfolder,"density_3",ptype)
ggsave(fn, plot=density3, units=punits, scale=pscale, width=pwidth, height=pheight)


varnames=c("Wage M", "Hours M", "Wage W", "Hours W")
plotlist <- get_density_plots(VCopObj, varnames=varnames, colorpalette=NULL)

# Density 4
density4 <- compare_plots(plotlist,
                          edge=4, groups=c(1975,2015),
                          colorpalette=rev(divcolors),
                          title=TRUE,
                          zlim=2,
                          bias=2)
density4
fn <- paste0(pfolder,"density_4",ptype)
ggsave(fn, plot=density4, units=punits, scale=pscale, width=pwidth, height=pheight)


# Hours women, wages men
density5 <- compare_plots(plotlist,
                          edge=5, groups=c(1975,2015),
                          colorpalette=rev(divcolors),
                          title=FALSE,
                          zlim=2,
                          bias=1.25)
density5
fn <- paste0(pfolder,"density_5",ptype)
ggsave(fn, plot=density5, units=punits, scale=pscale, width=pwidth, height=pheight)


# Hours women, hours men
density6 <- compare_plots(plotlist,
                          edge=6, groups=c(1975,2015),
                          colorpalette=rev(divcolors),
                          title=FALSE,
                          zlim=2.5,
                          bias=1.5)
density6
# Save
fn <- paste0(pfolder,"density_6",ptype)
ggsave(fn, plot=density6, units=punits, scale=pscale, width=pwidth, height=pheight)

### Combine all plots 
# Unconditional copulas
cdplots <- grid.arrange(
  density1,
  density2,
  density3, 
  ncol=1, nrow=3,
  padding=0,
  heights=c(1,1,1),#unit(c(0.5,5,5,5),rep("cm",4)),
  widths=c(2)#unit(c(5,5),rep("cm",2))
)

fn <- paste0(pfolder,"densities_main_deco_1975-2015_unconditional",ptype)
ggsave(fn, plot=cdplots, units=punits, scale=pscale, width=pwidth, height=pheight*2.3)

cdplots <- grid.arrange(
  density4,
  density5,
  density6, 
  ncol=1, nrow=3,
  padding=0,
  heights=c(1,1,1),#unit(c(0.5,5,5,5),rep("cm",4)),
  widths=c(2)#unit(c(5,5),rep("cm",2))
)

fn <- paste0(pfolder,"densities_main_deco_1975-2015_conditional",ptype)
ggsave(fn, plot=cdplots, units=punits, scale=pscale, width=pwidth, height=pheight*2.3)


# Retrieve only plots of assortative mating and male hours/wages
varnames=c("Wages men", "Hours men", "Wages women", "Hours women")
plotlist <- get_density_plots(VCopObj, varnames=varnames, colorpalette=NULL)
plotlist_sel <- plotlist[c(2,1)] 


########################################################################
### Plot density between male wages and female hours

setwd(paste0(wd,datafiles_d))
load(paste0("copula_female_hours_male_wages_",years[1],"-",years[2],"_new.rda"))
setwd(wd)

varnames=c("Wages men", "Hours women")
plotlist <- get_density_plots(VCopObj_MWage_FHours, varnames=varnames, colorpalette=NULL)

# Female hours and male wages
density1 <- compare_plots(plotlist,
                          edge=1, groups=c(1975,2015),
                          colorpalette=rev(divcolors),
                          bias=1, 
                          zlim=2)
density1
fn <- paste0(pfolder,"density_female_hours_male_wages",ptype)
ggsave(fn, plot=density1, units=punits, scale=pscale, width=pwidth, height=pheight)

########################################################################
### Plot density between male/female wages, male hours/wage, 
### and male wages/female hours, respectively

plotlist_sel <- c(plotlist_sel, plotlist)


density1 <- compare_plots(plotlist_sel,
                          edge=1, groups=c(1975,2015),
                          colorpalette=rev(divcolors),
                          bias=2, 
                          title=TRUE)

density2 <- compare_plots(plotlist_sel,
                          edge=2, groups=c(1975,2015),
                          colorpalette=rev(divcolors),
                          zlim=5,
                          bias=2.5, 
                          title=FALSE)

density3 <- compare_plots(plotlist_sel,
                          edge=3, groups=c(1975,2015),
                          colorpalette=rev(divcolors),
                          bias=1, 
                          zlim=2,
                          title=FALSE)

# Combine plots
cdplots <- grid.arrange(
  density1,
  density2,
  density3, 
  ncol=1, nrow=3,
  padding=0,
  heights=c(1,1,1),#unit(c(0.5,5,5,5),rep("cm",4)),
  widths=c(2)#unit(c(5,5),rep("cm",2))
)

fn <- paste0(pfolder,"density_mating_male_wage_labor_division",ptype)
ggsave(fn, plot=cdplots, units=punits, scale=pscale, width=pwidth, height=pheight*2.3)
  
