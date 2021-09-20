##############################################################
## D1 - CPS - Descriptive statistics plots
##############################################################

##############################################################
# Marginal statistics all years

# Load data
setwd(paste0(wd,datafiles_d))
load("couple_stats_CPS_1976-2019.rda")
recessions <- read.csv2("nber_recessions.csv",header=TRUE)
setwd(wd)

# Recession years
recessions <- data.frame(peak=recessions[,3]/12+1800-1/12,
                         trough=recessions[,4]/12+1800-1/12)
recessions <- subset(recessions, recessions$peak>1975|recessions$trough>1975)

##############################################################
# Prep data
names(stats_all_years)[1] <- "EARNINGS"
stats_all_years2 <- melt(stats_all_years, 
                         id.vars=c("year","stat"),
                         measure.vars=names(stats_all_years)[1:5],
                         variable.name="var",
                         value.name="value") 
stats_all_years2$cat <- "Employed & non-active"
stats_all_years3  <- melt(stats_pos_all_years, 
                          id.vars=c("year","stat"),
                          measure.vars=names(stats_pos_all_years)[1:4],
                          variable.name="var",
                          value.name="value")
stats_all_years3$cat <- "Employed"
stats_all_years2 <- rbind(stats_all_years2,stats_all_years3)
rm("stats_all_years3")

stats_all_years2$gender <- memisc::recode(stats_all_years2$var, 
                                          "Women" <-c("WAGE_f","HOURS_f"),    
                                          "Men" <- c("WAGE_m","HOURS_m"),
                                          NA <-c("EARNINGS"))
stats_all_years2$var2 <- memisc::recode(stats_all_years2$var, 
                                        "Weekly hours" <-c("HOURS_m","HOURS_f"),    
                                        "Hourly wages" <- c("WAGE_m","WAGE_f"),
                                        "Earnings" <- c("EARNINGS"))                                           
stats_all_years2$year <- stats_all_years2$year-1 


##############################################################
## Figure 1: Earnings quantiles
sel <- which(is.element(stats_all_years2$stat,c("P10","P50","P90"))&stats_all_years2$var=="EARNINGS")
dfplot <- stats_all_years2[sel,]
dfplot$stat <- dplyr::recode(dfplot$stat,
                             "P10"="0.1",
                             "P50"="0.5",
                             "P90"="0.9")
ggplot() + geom_line(data=dfplot, aes(year,value, col=stat)) + 
  geom_point(data=dfplot, aes(year,value, col=stat)) + 
  geom_hline(yintercept=0, color="darkgrey") +
  labs(y="2019 US dollars", x="", color="Quantile")+
  geom_rect(data=recessions, aes(xmin=peak, xmax=trough, ymin=-Inf, ymax=+Inf), alpha=0.2, fill="grey") + 
  scale_y_continuous(labels = scales::number_format(big.mark=","))
ggsave(paste0(pfolder,"quantiles_earnings",ptype),scale=pscale,units=punits,width=pwidth,height=pheight)
ggsave(paste0(pfolder,"quantiles_earnings_paper",ptype),scale=pscale,units=punits,width=pwidth,height=pheight)


##############################################################
## Figure C.1: Wage quantiles
sel <- which(is.element(stats_all_years2$stat,c("P10","P50","P90"))&
               is.element(stats_all_years2$var,c("WAGE_m","WAGE_f")))
dfplot <- stats_all_years2[sel,]
dfplot$stat <- dplyr::recode(dfplot$stat,
                             "P10"="0.1",
                             "P50"="0.5",
                             "P90"="0.9")
dfplot$cat <- dplyr::recode(dfplot$cat,
                             "Employed"="Employed",
                             "Employed & non-active"="Employed &\nnon-active")
ggplot() + 
  geom_line(data=dfplot, aes(year,value, color=stat, linetype=cat)) + 
  geom_hline(yintercept=0, color="darkgrey") +
  facet_grid( ~ gender) + 
  labs(y="2019 US dollars", x="", color="Quantile", linetype="Population", shape="Population") + 
  geom_rect(data=recessions, aes(xmin=peak, xmax=trough, ymin=-Inf, ymax=+Inf), alpha=0.2, fill="grey")
ggsave(paste0(pfolder,"quantiles_wages",ptype), scale=pscale, units=punits,width=pwidth,height=pheight)
ggsave(paste0(pfolder,"quantiles_wages_paper",ptype), scale=pscale, units=punits, width=pwidth, height=pheight)

##############################################################
## Figure C.2: Employment rate and average hours of labor force participants
sel1 <- which(is.element(stats_all_years2$stat,c("%obs.=0"))&
               is.element(stats_all_years2$var,c("HOURS_m","HOURS_f"))&
               is.element(stats_all_years2$cat,c("Employed & non-active")))
sel2 <- which(is.element(stats_all_years2$stat,c("Mean"))&
                is.element(stats_all_years2$var,c("HOURS_m","HOURS_f"))&
                is.element(stats_all_years2$cat,c("Employed")))
dfplot <- stats_all_years2[c(sel1,sel2),]
dfplot[which(dfplot$stat=="%obs.=0"),"value"] <- 1-dfplot[which(dfplot$stat=="%obs.=0"),"value"]
dfplot$stat <- as.factor(dfplot$stat)  
levels(dfplot$stat) <- c("Employment Rate","Avg. Weekly Hours of Employees")
ggplot() + 
  geom_line(data=dfplot, aes(year,value, color=gender)) + 
  geom_point(data=dfplot, aes(year,value, color=gender)) +
  facet_wrap(~ stat, scales="free") + 
  labs(y=" ", x="", color="") + 
  geom_rect(data=recessions, aes(xmin=peak, xmax=trough, ymin=-Inf, ymax=+Inf), alpha=0.2, fill="grey")
ggsave(paste0(pfolder,"employment_rate_and_average_hours",ptype), scale=pscale, units=punits,width=pwidth,height=pheight)
ggsave(paste0(pfolder,"employment_rate_and_average_hours_paper",ptype), scale=pscale, units=punits,width=pwidth,height=pheight*0.9)


##############################################################
## Figure C.3: Corrrelations
dfplot <- corrs_all_years
dfplot$year <- dfplot$year-1 
dfplot$var1 <- substr(corrs_all_years$vars,1,regexpr("-",corrs_all_years$vars)-1) 
dfplot$var2 <- substr(corrs_all_years$vars,regexpr("-",corrs_all_years$vars)+1,nchar(as.character(corrs_all_years$vars))) 
dfplot$vars <- memisc::recode(corrs_all_years$vars, 
                              "Female Wage - Male Wage" <- "WAGE_f-WAGE_m",
                              "Female Hours - Male Hours" <- "HOURS_f-HOURS_m",
                              "Female Wage - Female Hours" <- "HOURS_f-WAGE_f",
                              "Male Wage - Male Hours" <- "HOURS_m-WAGE_m",
                              "Female Wage - Male Hours" <- "WAGE_f-HOURS_m",
                              "Male Wage - Female Hours" <- "HOURS_f-WAGE_m"
)

ggplot(dfplot, aes(year,rho, color=type)) + 
  geom_line() + geom_point() + geom_hline(yintercept=0, color="darkgrey") +
  facet_wrap( ~ vars) + 
  labs(y="", x="", color="Correlation")
ggsave(paste0(pfolder,"correlations",ptype),units=punits,width=pwidth,height=pheight*1.2)




##############################################################
# Additional plots (not used for paper)

# ## Gini plots
# # Earnings
# sel <- which(stats_all_years2$stat=="Gini"&stats_all_years2$var=="EARNINGS")
# dfplot <- stats_all_years2[sel,]
# ggplot() + geom_line(data=dfplot, aes(year,value, color=stat)) +
#   geom_point(data=dfplot, aes(year,value, color=stat)) +
#   labs(y="Gini", x="Year", color="") + 
#   theme(legend.position = "none")+
#   geom_rect(data=recessions, aes(xmin=peak, xmax=trough, ymin=-Inf, ymax=+Inf), alpha=0.2, fill="grey")
# ggsave(paste0(pfolder,"gini_earnings",ptype),units=punits,width=pwidth,height=pheight, bg=pbg)
# 
# # Wages and hours
# sel <- which(stats_all_years2$stat=="Gini"&stats_all_years2$var!="EARNINGS")
# dfplot <- stats_all_years2[sel,]
# ggplot() + geom_line(data=dfplot, aes(year,value, color=cat)) +
#   geom_point(data=dfplot, aes(year,value, color=cat)) + 
#   facet_grid(var2 ~ gender) +
#   labs(y="Gini", x="Year", color="Empl. Status")  + 
#   geom_rect(data=recessions, aes(xmin=peak, xmax=trough, ymin=-Inf, ymax=+Inf), alpha=0.2, fill="grey")
# ggsave(paste0(pfolder,"gini_wage_hours",ptype),scale=pscale,units=punits,width=pwidth,height=pheight)


# # Quantile hours all
# sel <- which(is.element(stats_all_years2$stat,c("P10","P25","P50","P90"))&
#                is.element(stats_all_years2$cat,c("Employed & non-active"))&
#                is.element(stats_all_years2$var,c("HOURS_m","HOURS_f"))
# )
# dfplot <- stats_all_years2[sel,]
# ggplot() + 
#   geom_line(data=dfplot, aes(year,value, color=stat)) + 
#   geom_point(data=dfplot, aes(year,value, color=stat)) + 
#   geom_hline(yintercept=0, color="lightgrey") +
#   facet_grid(~ gender) + 
#   labs(y="Avg. Weekly Hours", x="Year", color="Quantile") + 
#   geom_rect(data=recessions, aes(xmin=peak, xmax=trough, ymin=-Inf, ymax=+Inf), alpha=0.2, fill="grey")
# ggsave(paste0(pfolder,"quantiles_hours",ptype),scale=pscale, units=punits,width=pwidth,height=pheight)
# 
# # Quantile hours by employment status
# sel <- which(is.element(stats_all_years2$stat,c("P10","P50"))&
#                #is.element(stats_all_years2$cat,c("Employed & non-active"))&
#                is.element(stats_all_years2$var,c("HOURS_m","HOURS_f"))
# )
# dfplot <- stats_all_years2[sel,]
# ggplot() + 
#   geom_line(data=dfplot, aes(year,value, color=stat, linetype=cat)) + 
#   geom_point(data=dfplot, aes(year,value, color=stat, shape=cat)) + 
#   geom_hline(yintercept=0, color="lightgrey") +
#   facet_grid(~ gender) + 
#   #labs(y="Avg. Weekly Hours", x="Year", color="Quantile") + 
#   labs(y="Avg. Weekly Hours", x="Year", color="Quantile", linetype="Empl. Status", shape="Empl. Status") + 
#   geom_rect(data=recessions, aes(xmin=peak, xmax=trough, ymin=-Inf, ymax=+Inf), alpha=0.2, fill="grey")
# ggsave(paste0(pfolder,"quantiles_hours_cat",ptype),scale=pscale,units=punits,width=pwidth,height=pheight)
# 
# # Zero hours
# sel <- which(is.element(stats_all_years2$stat,c("%obs.=0"))&
#                is.element(stats_all_years2$var,c("HOURS_m","HOURS_f"))&
#                is.element(stats_all_years2$cat,c("Employed & non-active")))
# dfplot <- stats_all_years2[sel,]
# ggplot() + 
#   geom_line(data=dfplot, aes(year,value*100, color=cat)) + 
#   geom_point(data=dfplot, aes(year,value*100, color=cat)) +
#   geom_hline(yintercept=0, color="lightgrey") +
#   facet_grid(~ gender) + 
#   theme(legend.position = "none")+
#   labs(y="% Not Employed", x="Year", color="share") + 
#   geom_rect(data=recessions, aes(xmin=peak, xmax=trough, ymin=-Inf, ymax=+Inf), alpha=0.2, fill="grey")
# ggsave(paste0(pfolder,"zero_hours_lf",ptype),units=punits,width=pwidth,height=pheight)


# ##############################################################
# # Compare two years
# 
# # Load data
# years <- c(1976,2016)
# setwd(paste0(wd,datafiles_d))
# load(paste0("marginal_stats_",years[1],"-",years[2],".rda"))
# setwd(wd)
# 
# # Quantile functions
# plotdf <- marginals_stats$estimates$quants
# names(plotdf)[3:7] <- c("Couples' Earnings",
#                         "Male Wages", "Male Hours", 
#                         "Female Wages", "Female Hours")
# plotdf <- melt(plotdf, id.vars=names(plotdf)[1:2],
#                measure.vars=names(plotdf)[3:7],
#                variable.name = "var",
#                value.name = "y")
# 
# plotdf$group <- factor(plotdf$group)
# levels(plotdf$group) <- years
# lev <- levels(plotdf$var)
# lev
# 
# if(is.null(marginals_stats$se)){
#   plotdf$se <- 0
# }else{
#   plotdfse <- marginals_stats$se$quants_se
#   plotdfse <- melt(plotdfse, id.vars=names(plotdfse)[1:2],
#                    measure.vars=names(plotdfse)[3:7],
#                    variable.name = "var",
#                    value.name = "y")
#   plotdf$se <- plotdfse$y
# }
# 
# 
# plot_hours <- ggplot(plotdf[which(is.element(plotdf$var,lev[c(3,5)])),], aes(tau,y, col=group, fill=group)) +
#   geom_line() + geom_point() + 
#   geom_ribbon(aes(ymin = ifelse(y - 1.96*se>0,y - 1.96*se,0), ymax = y + 1.96*se),alpha=0.3) + 
#   facet_grid( ~ var) +labs(y="Weekly hours", x="Quantile", color="Year", fill="Year")
# 
# plot_wages <- ggplot(plotdf[which(is.element(plotdf$var,lev[c(2,4)])),], aes(tau,log(y), col=group, fill=group)) +
#   geom_line() + geom_point() + 
#   geom_ribbon(aes(ymin = log(y - 1.96*se), ymax = log(y + 1.96*se)),alpha=0.3) + 
#   facet_grid( ~ var) +labs(y="Log hourly wage", x="Quantile", color="Year", fill="Year")
# 
# plot_earnings <- ggplot(plotdf[which(is.element(plotdf$var,lev[c(1)])),], aes(tau,log(y), col=group, fill=group)) +
#   geom_line() + geom_point() + 
#   geom_ribbon(aes(ymin = log(y - 1.96*se), ymax = log(y + 1.96*se)),alpha=0.3) + 
#   #facet_wrap( ~ var, ncol=1) + 
#   labs(y="Labor earnings", x="Quantile", color="Year", fill="Year")
# 
# plot_earnings2 <- ggplot(plotdf[which(is.element(plotdf$var,lev[c(1)])),], aes(tau,y, col=group, fill=group)) +
#   geom_line() + geom_point() + 
#   geom_ribbon(aes(ymin = y - 1.96*se, ymax = y + 1.96*se),alpha=0.3) + 
#   #facet_wrap( ~ var, ncol=1) + 
#   labs(y="Labor earnings", x="Quantile", color="Year", fill="Year")
# 
# 
# plot_hours_wages  <- ggplot(plotdf[which(is.element(plotdf$var,lev[2:5])),], aes(tau,y, col=group, fill=group)) +
#   geom_line() + geom_point() +
#   facet_wrap( ~ var,scales = "free", dir="v", labeller = ) + 
#   labs(y="Hourly wages/Weekly hours", x="Quantile", color="Year", fill="Year")
# 
# #plot_hours_wages <- grid.arrange(plot_hours, plot_wages, plot_earnings)
# plot_hours_wages <- grid.arrange(plot_hours, plot_wages)
# 
# plot_hours_wages
# 
# # Save
# fn <- paste0(pfolder,"quantiles_marginals",ptype)
# ggsave(fn, units=punits, width=pwidth, height=pheight)
# 
# ##############################################################
# # Correlations plots
# 
# # Correlation plots
# plotdf <- marginals_stats$estimates$corrs[order(marginals_stats$estimates$corrs$group,
#                                                 marginals_stats$estimates$corrs$vars,
#                                                 marginals_stats$estimates$corrs$type),]
# if(is.null(marginals_stats$se)){
#   plotdf$se <- 0
# }else{
#   plotdf$se <- marginals_stats$se$corrs[order(marginals_stats$se$corrs$group,
#                                               marginals_stats$se$corrs$vars,
#                                               marginals_stats$se$corrs$type),"rho"]
# }
# 
# plotdf$group <- factor(plotdf$group)
# levels(plotdf$group) <- c("1979","2019")
# 
# # All association measures
# ggplot(plotdf, aes(type,rho, fill=group)) + 
#   geom_bar(stat="identity", position=position_dodge())+ 
#   #facet_grid(vars ~ cat )
#   geom_errorbar(aes(ymin=rho-1.96*se, ymax=rho+1.96*se), width=.2,
#                 position=position_dodge(0.9))+ 
#   facet_wrap(~vars)
# 
# ## Spearman's rho only
# plotdf$vars2 <- plotdf$vars
# 
# stat_sel <- "rank cor."
# 
# vs <- c("Wages Men","Hours Men","Wages Women","Hours Women")
# vs <- c(paste0(vs[1],"-",vs[2:4]),paste0(vs[2],"-",vs[3:4]),paste0(vs[3],"-",vs[4]))
# levels(plotdf$vars2) <- vs[c(5,4,6,3,1,2)]
# plotdf$order1 <- rank(abs(plotdf[which(plotdf$type==stat_sel&plotdf$group==levels(plotdf$group)[2]),"rho"])) # 6:1  
# ggplot(plotdf[which(plotdf$type==stat_sel),], aes(reorder(vars2,order1),rho,fill=group,order)) + 
#   #geom_bar(stat="identity",  width=0.9, position=position_dodge(0.9))+ 
#   geom_bar(stat="identity",  width=0.8, position=position_dodge2(0.9,reverse = TRUE)) +
#   geom_hline(yintercept = 0, col="darkgrey") + 
#   geom_errorbar(aes(ymin=rho-1.96*se, ymax=rho+1.96*se),color="black", width=0.7,
#                 #               position=position_dodge(0.9)) + 
#                 position=position_dodge2(1.5, padding=0.5, reverse = TRUE))+ 
#   coord_flip() +  labs(x="",y="Spearman's rho", fill="Year") 
# 
# # Save
# fn <- paste0(pfolder,"spearmansrho",ptype)
# ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

