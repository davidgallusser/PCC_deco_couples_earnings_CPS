##############################################################
## D3 - CPS - deco quantile plots 1976 - 2016 - baseline decomposition
##############################################################

years <- c(1976,2016)
setwd(paste0(wd,datafiles_d))
load(paste0("copula_decomposition_main_",years[1],"-",years[2],".rda"))
setwd(wd)

##############################################################
# Prepare data, rename variables
plotdf <- VCopObj$VCdeco$aggregate$quants_deco
plotdf <- melt(plotdf, id.vars=names(plotdf)[1:4],
               measure.vars=names(plotdf)[5:ncol(plotdf)],
               variable.name = "tau",
               value.name = "delta")
plotdf$tau <- as.numeric(as.character(plotdf$tau))

# Add s.e. (if available)
if(is.null(VCopObj$VCdeco_se$aggregate)==FALSE){
  plotdf2 <- VCopObj$VCdeco_se$aggregate$quants_deco_se
  plotdf2 <- plotdf2[match(VCopObj$VCdeco$aggregate$quants_deco$names1,plotdf2$names1),]
  plotdf2 <- melt(plotdf2, id.vars=names(plotdf2)[1:4],
                  measure.vars=names(plotdf2)[5:ncol(plotdf2)],
                  variable.name = "tau",
                  value.name = "se")
  plotdf$se <- plotdf2$se
  plotdf$t <- VCopObj$VCdeco_se$aggregate$quants_deco_critval[match(plotdf[,1],VCopObj$VCdeco_se$aggregate$quants_deco_critval$names1),"t"]
}else{
  plotdf$se <- 0 
  plotdf$t <- 0
}

plotdf$names2 <- dplyr::recode_factor(plotdf$names2,
                                      "1 Total difference"="Observed difference",
                                      "2 Marginals"="Aggregate marginal effect",
                                      "3 Copula"="Aggregate dependence effect",
                                      "1 WAGE_m"="Wage men",
                                      "2 ANNUALHOURS_m"="Hours men",
                                      "3 WAGE_f"="Wage women",
                                      "4 ANNUALHOURS_f"="Hours women",
                                      "1 WAGE_m,ANNUALHOURS_m"="Wage men, hours men",
                                      "2 WAGE_m,WAGE_f"="Wage men, wage women",
                                      "3 WAGE_f,ANNUALHOURS_f"="Wage women, hours women",
                                      "4 ANNUALHOURS_m,WAGE_f;WAGE_m"="Hours men, wage women\n|wage men", 
                                      "5 WAGE_m,ANNUALHOURS_f;WAGE_f"="Wage men, hours women\n|wage women",
                                      "6 ANNUALHOURS_m,ANNUALHOURS_f;WAGE_m,WAGE_f"="Hours men, hours women\n|wage men, wage women")
                                      


       
##############################################################
## Aggregate decomposition (Panel A)
plot_agg <- ggplot(plotdf[which(plotdf$cat2=="aggregate"),],
                   aes(tau,delta,
                       col=names2,
                       fill=names2,
                       shape=names2)) + 
  geom_hline(yintercept=0, col="darkgrey") + 
  #facet_wrap(~names2, ncol=3) + theme(legend.position = "none")+ 
  labs(x="Quantile",y="Log difference", color="", fill="", shape="") + 
  geom_ribbon(aes(ymin = delta - t*se, ymax = delta + t*se),  alpha=0.3, color=NA) + 
  geom_line() + geom_point() 
plot_agg

fn <- paste0(pfolder,"deco_1976_2016_main_quantiles_aggregate_new",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)


##############################################################
## Detailed decomposition: Direct marginal effects (Panel B)
plot_marg1 <- ggplot(plotdf[which(plotdf$cat1=="Marginal"&
                                   plotdf$cat2=="direct"&
                                   plotdf$names2!="5 M Interactions > 0"&
                                   plotdf$tau>0.04),], 
                                   aes(tau,delta,
                                       col=names2,
                                       fill=names2, 
                                       shape=names2
                                       )) + 
  geom_hline(yintercept=0, col="darkgrey") + 
  #geom_ribbon(aes(ymin = delta - 1.96*se, ymax = delta + 1.96*se),alpha=0.3) +
  geom_ribbon(aes(ymin = delta - t*se, ymax = delta + t*se),  alpha=0.4, color=NA) + 
  geom_line() + geom_point() + 
  #facet_wrap(~names2, ncol=2) + theme(legend.position = "none") + 
  labs(x="Quantile",y="Log difference", color="", fill="", shape="") +
  coord_cartesian(ylim = c(-0.3, 1.1))
plot_marg1 

fn <- paste0(pfolder,"deco_1976_2016_main_quantiles_detailed_marginals1",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)


##############################################################
## Detailed decomposition: Direct dependence effects 1 (Panel C)
plot_cop1 <- ggplot(plotdf[which(plotdf$cat1=="Copula"&
                                  plotdf$cat2=="direct"&
                                  is.element(plotdf$names1,c("C 1","C 2","C 3","C 5"))&
                                  plotdf$tau>0.04),], 
                   aes(tau,delta,col=names2,
                       fill=names2,
                       shape=names2)) + 
  geom_hline(yintercept=0, col="darkgrey") + 
  #geom_ribbon(aes(ymin = delta - 1.96*se, ymax = delta + 1.96*se), alpha=0.3)+
  geom_ribbon(aes(ymin = delta - t*se, ymax = delta + t*se), alpha=0.1, color=NA)+
  geom_line() + geom_point() + 
  #facet_wrap(~names2, ncol=2, scales = "free") + theme(legend.position = "none")+ 
  labs(x="Quantile",y="Log difference", color="", fill="", shape="") + 
  coord_cartesian(ylim = c(-0.15, 0.15))
plot_cop1 

fn <- paste0(pfolder,"deco_1976_2016_main_quantiles_detailed_copula1",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)


##############################################################
## Detailed decomposition: Direct dependence effects (Panel D)
plot_cop2 <- ggplot(plotdf[which(plotdf$cat1=="Copula"&
                                  plotdf$cat2=="direct"&
                                  is.element(plotdf$names1,c("C 4","C 6"))&
                                  plotdf$tau>0.04),],
                   aes(tau,delta,col=names2,
                       fill=names2,
                       shape=names2)) + 
  geom_hline(yintercept=0, col="darkgrey") + 
  #geom_ribbon(aes(ymin = delta - 1.96*se, ymax = delta + 1.96*se), alpha=0.3)+
  geom_ribbon(aes(ymin = delta - t*se, ymax = delta + t*se), alpha=0.1, color=NA)+
  geom_line() + geom_point() + 
  #facet_wrap(~names2, ncol=2, scales = "free") + theme(legend.position = "none")+ 
  labs(x="Quantile",y="Log difference", color="", fill="", shape="") + 
  coord_cartesian(ylim = c(-0.1, 0.1))
plot_cop2

fn <- paste0(pfolder,"deco_1976_2016_main_quantiles_detailed_copula2",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

##############################################################
## Detailed decomposition: Marginal and Dependence Interaction Effect (Panel E)
dfplot <- plotdf[which(plotdf$names2 %in% c("5 M Interactions > 0", "7 C Interactions > 0")&
                         plotdf$tau>0.04),]
dfplot$names2 <- dplyr::recode(dfplot$names2,
                               "5 M Interactions > 0"="Marginal interaction",
                               "7 C Interactions > 0"="Dependence interaction")
plot_inter <- ggplot(dfplot, 
                     aes(tau,delta,
                         col=names2,
                         fill=names2, 
                         shape=names2
                     )) + 
  geom_hline(yintercept=0, col="darkgrey") + 
  #geom_ribbon(aes(ymin = delta - 1.96*se, ymax = delta + 1.96*se),alpha=0.3) +
  geom_ribbon(aes(ymin = delta - t*se, ymax = delta + t*se),  alpha=0.4, color=NA) + 
  geom_line() + geom_point() + 
  labs(x="Quantile",y="Log difference", color="", fill="", shape="") #+
#coord_cartesian(ylim = c(-0.3, 0.3))
plot_inter 

fn <- paste0(pfolder,"deco_1976_2016_main_quantiles_detailed_interactions",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)


##############################################################
## Detailed decomposition: Combine most important results
labels <- c("A Aggregate decomposition",
            "B Direct marginal effects",
            "C Most important direct dependence effects")

plot_agg <- plot_agg +  
            theme(plot.margin= unit(c(0.9,0,0.1,0), "cm"))  + 
            scale_y_continuous(labels=function(x) scales::comma(x, accuracy = 0.01)) + 
            labs(y="Log difference", x="", shape="", col="", fill="")
plot_marg1 <- plot_marg1 +  
              theme(plot.margin= unit(c(0.9,0,0.1,0), "cm")) +
              scale_y_continuous(labels=function(x) scales::comma(x, accuracy = 0.01)) + 
              labs(y="Log difference", x="", shape="", col="", fill="")
plot_cop1 <- plot_cop1 + 
             theme(plot.margin= unit(c(0.9,0,0,0), "cm"))

plot_summary <- egg::ggarrange(plot_agg, 
               plot_marg1, 
               plot_cop1,
               labels = labels, 
               label.args = list(gp=gpar(font=2), x=unit(0.5,"line"), hjust=0, vjust=2))
plot_summary

fn <- paste0(pfolder,"deco_1976_2016_main_quantiles_summary",ptype)
ggsave(fn, plot=plot_summary, units=punits, scale=pscale, width=pwidth, height=pheight*2.45)


##############################################################
## Detailed decomposition: Combine remaining results
labels <- c("D Remaining direct dependence effects",
            "E Interaction effects")

plot_cop2 <- plot_cop2 +  
  theme(plot.margin= unit(c(0.9,0,0.1,0), "cm"))  + 
  scale_y_continuous(labels=function(x) scales::comma(x, accuracy = 0.01)) + 
  labs(y="Log difference", x="", shape="", col="", fill="")
plot_inter  <- plot_inter  + 
  theme(plot.margin= unit(c(0.9,0,0,0), "cm"))

plot_summary2 <- egg::ggarrange(plot_cop2, 
                               plot_inter,
                               labels = labels, 
                               label.args = list(gp=gpar(font=2), x=unit(0.5,"line"), hjust=0, vjust=2))
plot_summary2

fn <- paste0(pfolder,"deco_1976_2016_main_quantiles_summary_remaining",ptype)
ggsave(fn, plot=plot_summary2, units=punits, scale=pscale, width=pwidth, height=pheight*1.75)

###############################################################
## Additional plots
# 
# 
# ##############################################################
# ## Detailed decomposition: Marginals Interaction
# plot_marg2 <- ggplot(plotdf[which(plotdf$cat1=="Marginal"&plotdf$cat2=="interaction"&plotdf$tau>0.04),],
#                      aes(tau,
#                          delta,
#                          col=names2,
#                          fill=names2,
#                          #shape=names2
#                      )) + 
#   geom_hline(yintercept=0, col="lightgrey") + 
#   #geom_ribbon(aes(ymin = delta - t*se, ymax = delta + t*se),  alpha=0.4, color=NA) + 
#   #facet_wrap(~names2, ncol=2) + theme(legend.position = "none")+ 
#   geom_line() + geom_point() + 
#   labs(x="Quantile",y="Log difference", color="", fill="", shape="") +
#   coord_cartesian(ylim = c(-0.1, 0.1))
# plot_marg2 
# 
# fn <- paste0(pfolder,"deco_1976_2016_main_quantiles_detailed_marginals2",ptype)
# ggsave(fn, units=punits, width=pwidth, height=pheight)
# 
# 
# 
# ##############################################################
# ## Detailed decomposition: Copula 3
# plot_cop3 <- ggplot(plotdf[which(plotdf$cat1=="Copula"&
#                                    plotdf$cat2=="direct"&
#                                    is.element(plotdf$names1,c("C1","C2","C3","C5"))&
#                                    plotdf$tau>0.04),], aes(tau,delta,col=names2,fill=names2)) + 
#   geom_line() + geom_point() + 
#   geom_hline(yintercept=0, col="darkgrey") + 
#   facet_wrap(~names2, ncol=2, scales = "free") + theme(legend.position = "none")+ 
#   labs(x="Quantile",y="Log difference", color="Effect", fill="Effect") + 
#   coord_cartesian(ylim = c(-0.1, 0.1))
# plot_cop3
# 
# # Add confidence bands
# plot_cop   <- plot_cop +  
#   #geom_ribbon(aes(ymin = delta - 1.96*se, ymax = delta + 1.96*se), alpha=0.3)+
#   geom_ribbon(aes(ymin = delta - t*se, ymax = delta + t*se), alpha=0.4, color=NA)
# plot_cop
# 
# fn <- paste0(pfolder,"deco_quantiles_detailed_copula3",ptype)
# ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)
# 
# ##############################################################
# ## Detailed decomposition: Copula Interaction
# plot_cop <- ggplot(plotdf[which(plotdf$cat1=="Copula"&plotdf$cat2=="interaction"&plotdf$tau>0.04),], aes(tau,delta,col=names2,fill=names2)) + 
#   geom_line() + geom_point() + 
#   geom_hline(yintercept=0, col="lightgrey") + 
#   facet_wrap(~names2, ncol=4, scales = "free") + theme(legend.position = "none")+ 
#   labs(x="Quantile",y="Log difference", color="Effect", fill="Effect") + 
#   coord_cartesian(ylim = c(-0.05, 0.05))
# plot_cop 
# 
# sel <- levels(plotdf$names2)[13]
# ggplot(plotdf[which(plotdf$names2==sel),], aes(tau,delta,col=names2,fill=names2)) + 
#   geom_line() + geom_point() + 
#   geom_hline(yintercept=0, col="lightgrey") + 
#   #facet_wrap(~names2, ncol=2) + theme(legend.position = "none")+ 
#   labs(x="Quantile",y="Log difference", color="Effect", fill="Effect") #+ 
# #coord_cartesian(ylim = c(-0.1, 0.1))

