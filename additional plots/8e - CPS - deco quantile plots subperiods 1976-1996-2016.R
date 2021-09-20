##############################################################
## 8 - CPS - deco quantile plots subperiods  1976 - 1996 - 2016
##############################################################

years <- c(1976,1996)
setwd(paste0(wd,datafiles_d))
load(paste0("copula_decomposition_",years[1],"-",years[2],".rda"))
setwd(sav_d)

VCopObj_7696 <- VCopObj

years <- c(1996,2016)
setwd(paste0(wd,datafiles_d))
load(paste0("copula_decomposition_",years[1],"-",years[2],".rda"))
setwd(sav_d)

VCopObj_9616 <- VCopObj
rm("VCopObj")

##############################################################
# Decomposition quantile plots
plotdf1 <- VCopObj_7696$VCdeco$aggregate$quants_deco
plotdf2 <- VCopObj_9616$VCdeco$aggregate$quants_deco
plotdf1$effect <- "1976-1996"
plotdf2$effect <- "1996-2016"
plotdf <- rbind(plotdf1,plotdf2)

plotdf <- plotdf[, c(1:4,ncol(plotdf),6:(ncol(plotdf)-1))]
plotdf <- melt(plotdf, id.vars=names(plotdf)[1:5],
               measure.vars=names(plotdf)[6:ncol(plotdf)],
               variable.name = "tau",
               value.name = "delta")
plotdf$tau <- as.numeric(as.character(plotdf$tau))

# Add s.e. (if available)
if(is.null(VCopObj_cond$VCdeco_se$aggregate)==FALSE){
  plotdf2a <- VCopObj_7696$VCdeco_se$composition_effect$quants_deco_se
  plotdf2b <- VCopObj_9616$VCdeco_se$wage_structue_effect$quants_deco_se
  plotdf2a$year <- "1976-1996"
  plotdf2b$year <- "1996-2016"
  plotdf2 <- rbind(plotdf2a,plotdf2b)
  plotdf2 <- plotdf2[, c(1:4,ncol(plotdf2),5:(ncol(plotdf2)-1))]
  plotdf2 <- melt(plotdf2, id.vars=names(plotdf2)[1:5],
                  measure.vars=names(plotdf2)[6:ncol(plotdf2)],
                  variable.name = "tau",
                  value.name = "se")
  plotdfa <- subset(plotdf,effect =="1976-1996")
  plotdf2a <- subset(plotdf2,effect =="1976-1996")
  plotdfb <- subset(plotdf,effect !="1976-1996")
  plotdf2b <- subset(plotdf2,effect !="1976-1996")
  plotdfa$se <- plotdf2a$se[match(plotdfa$names1,plotdf2a$names1)]
  plotdfb$se <- plotdf2b$se[match(plotdfb$names1,plotdf2b$names1)]
  plotdfa$t <- VCopObj_7696$VCdeco_se$aggregate$quants_deco_critval[match(plotdfa[,1],VCopObj_7696$VCdeco_se$composition_effect$quants_deco_critval$names1),"t"]
  plotdfb$t <- VCopObj_9616$VCdeco_se$aggregate$quants_deco_critval[match(plotdfb[,1],VCopObj_9616$VCdeco_se$wage_structure_effect$quants_deco_critval$names1),"t"]
  plotdf <- rbind(plotdfa, plotdfb)
}else{
  plotdf$se <- 0 
  plotdf$t <- 0
}

plotdf$names2 <- dplyr::recode(plotdf$names2,
                               "1 WAGE_m" = "Wage Men",
                               "2 ANNUALHOURS_m" = "Hours Men",
                               "3 WAGE_f" = "Wage Women",
                               "4 ANNUALHOURS_f" = "Hours Women",
                               "1 WAGE_m,ANNUALHOURS_m" = "Wage Men,Hours Men",
                               "2 WAGE_m,WAGE_f" = "Wage Men, Wage Women",
                               "3 WAGE_f,ANNUALHOURS_f" = "Wage Women,Hours Women",
                               "4 ANNUALHOURS_m,WAGE_f;WAGE_m" = "Hours Men,Wage Women;Wage Men",
                               "5 WAGE_m,ANNUALHOURS_f;WAGE_f" = "Wage Men, Hours Women; Wage Women",
                               "6 ANNUALHOURS_m,ANNUALHOURS_f;WAGE_m,WAGE_f" = "Hours Men,Hours Women;Wage Men,Wage Women")
                               
##############################################################
## Aggregate decomposition
plot_agg <- ggplot(plotdf[which(plotdf$cat2=="aggregate"),],
                   aes(tau,delta,
                       col=effect,
                       fill=effect)) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept=0, col="darkgrey") + 
  facet_wrap(~names2, ncol=3)+  #+ theme(legend.position = "none")+ 
  labs(x="Quantile",y="Log difference",
       col="", fill="") 
plot_agg

# Add pointwise- and uniform-confidence bands
plot_agg <-  plot_agg +  
  #geom_ribbon(aes(ymin = delta - 1.96*se, ymax = delta + 1.96*se), alpha=0.3)  +
  geom_ribbon(aes(ymin = delta - t*se, ymax = delta + t*se),  alpha=0.3, color=NA)
plot_agg 

fn <- paste0(pfolder,"deco_quantiles_aggregate_subperiods",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

##############################################################
## Detailed decomposition: Marginals direct effects
plotdf1 <- subset(plotdf, cat1=="Marginal"&
                          cat2=="direct"&
                          names2!="5 M Interactions > 0"&
                          tau>0.04)
plot_marg <- ggplot(plotdf1 , 
                    aes(tau,delta,col=effect,fill=effect)) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept=0, col="darkgrey") + 
  facet_wrap(~names2, ncol=2) +# theme(legend.position = "none")+ 
  labs(x="Quantile",y="Log difference", color="", fill="")#+
  #coord_cartesian(ylim = c(-0.3, 0.3))
plot_marg 

# Add confidence bands
plot_marg  <- plot_marg + 
  #geom_ribbon(aes(ymin = delta - 1.96*se, ymax = delta + 1.96*se),alpha=0.3) +
  geom_ribbon(aes(ymin = delta - t*se, ymax = delta + t*se),  alpha=0.4, color=NA)
plot_marg 

fn <- paste0(pfolder,"deco_quantiles_detailed_marginals1_subperiods",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

## Detailed decomposition: Marginals Interaction
plotdf1 <- subset(plotdf, cat1=="Marginal"&
                    cat2=="interaction"&
                    names2!="5 M Interactions > 0"&
                    tau>0.04)

plot_marg <- ggplot(plotdf1,
                    aes(tau,delta,
                        col=effect,fill=effect)) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept=0, col="lightgrey") + 
  facet_wrap(~names2, ncol=2) + #theme(legend.position = "none")+ 
  labs(x="Quantile",y="Log difference", color="Contribution", fill="Contribution") +
  coord_cartesian(ylim = c(-0.1, 0.1))
plot_marg 

fn <- paste0(pfolder,"deco_quantiles_detailed_marginals2_subperiods",ptype)
ggsave(fn, units=punits, width=pwidth, height=pheight)

##############################################################
## Detailed decomposition: Copula 1
plotdf1 <- subset(plotdf, 
                   names1 %in% c("C 1","C 2","C 3","C 4","C 5","C 6") &
                    tau>0.04)

plot_cop <- ggplot(plotdf1,
                   aes(tau,delta,
                   col=effect,fill=effect)) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept=0, col="lightgrey") + 
  facet_wrap(~names2, ncol=2, as.table=FALSE) + #theme(legend.position = "none")+ 
  labs(x="Quantile",y="Log difference", color="", fill="") +
  coord_cartesian(ylim = c(-0.1, 0.1))
plot_cop 

# Add confidence bands
plot_cop   <- plot_cop +  
  #geom_ribbon(aes(ymin = delta - 1.96*se, ymax = delta + 1.96*se), alpha=0.3)+
  geom_ribbon(aes(ymin = delta - t*se, ymax = delta + t*se), alpha=0.4, color=NA)
plot_cop

fn <- paste0(pfolder,"deco_quantiles_detailed_copula_subperiods",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

##############################################################
## Detailed decomposition: Copula Interaction
plotdf1 <- subset(plotdf, 
                  cat1=="Copula"&
                  cat2=="interaction"&
                  tau>0.04)

plot_cop <- ggplot(plotdf1,
                   aes(tau,delta,
                       col=effect,fill=effect)) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept=0, col="lightgrey") + 
  facet_wrap(~names2, ncol=2, as.table=FALSE) + #theme(legend.position = "none")+ 
  labs(x="Quantile",y="Log difference", color="Contribution", fill="Contribution") +
  coord_cartesian(ylim = c(-0.05, 0.05))
plot_cop 

sel <- unique(plotdf$names2)[13]
ggplot(plotdf[which(plotdf$names2==sel),], aes(tau,delta,col=effect,fill=effect)) + 
  geom_line() + geom_point() + 
  geom_hline(yintercept=0, col="lightgrey") + 
  #facet_wrap(~names2, ncol=2) + theme(legend.position = "none")+ 
  labs(x="Quantile",y="Log difference", color="", fill="") #+ 
#coord_cartesian(ylim = c(-0.1, 0.1))


