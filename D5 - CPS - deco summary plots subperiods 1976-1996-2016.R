##############################################################
## D5 - CPS - deco summary plots subperiods  1976 - 1996 - 2016
##############################################################

years <- c(1976,1996)
setwd(paste0(wd,datafiles_d))
load(paste0("copula_decomposition_main_",years[1],"-",years[2],".rda"))
#load(paste0("copula_decomposition_",years[1],"-",years[2],".rda"))
setwd(wd)

VCopObj_7696 <- VCopObj

years <- c(1996,2016)
setwd(paste0(wd,datafiles_d))
load(paste0("copula_decomposition_main_",years[1],"-",years[2],".rda"))
#load(paste0("copula_decomposition_",years[1],"-",years[2],".rda"))
setwd(wd)

VCopObj_9616 <- VCopObj
rm("VCopObj")


#######################################################################
## Summary bar plot for decomposition 
#######################################################################

plotdf3 <- VCopObj_7696$VCdeco$aggregate$stats_deco
plotdf4 <- VCopObj_9616$VCdeco$aggregate$stats_deco
plotdf3$years <- "1975-1995"
plotdf4$years  <- "1995-2015"

plotdf3 <- rbind(plotdf3,plotdf4)
plotdf3 <-  subset(plotdf3,cat2=="direct"|cat1=="Observed")

# Rename variables
plotdf3$names2 <- as.character(plotdf3$names2)
plotdf3$names2 <- gsub("WAGE_m","Wage men", plotdf3$names2)
plotdf3$names2 <- gsub("WAGE_f","Wage women", plotdf3$names2)
plotdf3$names2 <- gsub("ANNUALHOURS_f","Hours women", plotdf3$names2)
plotdf3$names2 <- gsub("ANNUALHOURS_m","Hours men", plotdf3$names2)
plotdf3$names2 <- substr(plotdf3$names2,3,nchar(plotdf3$names2))
plotdf3$names2 <- paste0(1:length(plotdf3$names2)," ",plotdf3$names2)
plotdf3$names2[c(9,10,11,22,23,24)]<- rep(c("9 Hours M,Wage W; Wage M", 
                                            "10 Wage M,Hours W; Wage W",
                                            "11 Hours M,Hours W; Wage M,Wage W"),2) 
plotdf3$names2 <- substr(plotdf3$names2,3,nchar(plotdf3$names2))
plotdf3$names2 <- stringr::str_trim(plotdf3$names2, side='both')
plotdf3$names2 <- gsub(";"," |", plotdf3$names2)

plotdf3 <- plotdf3[, c(1:4,ncol(plotdf3),5:(ncol(plotdf3)-1))]
plotdf3 <- melt(plotdf3, id.vars=names(plotdf3)[1:5],
                measure.vars=names(plotdf3)[6:ncol(plotdf3)],
                variable.name = "stat",
                value.name = "delta")

if(is.null(VCopObj_7696$VCdeco_se$aggregate)==FALSE){
  plotdf2a <- VCopObj_7696$VCdeco_se$aggregate$stats_deco
  plotdf2b <- VCopObj_9616$VCdeco_se$aggregate$stats_deco
  plotdf2a$years <- "1975-1995"
  plotdf2b$years  <- "1995-2015"
  plotdf2 <- rbind(plotdf2a,plotdf2b)
  plotdf2 <-  subset(plotdf2,cat2=="direct"|cat1=="Observed")
  plotdf2 <- plotdf2[, c(1:4,ncol(plotdf2),5:(ncol(plotdf2)-1))]
  plotdf2 <- melt(plotdf2, id.vars=names(plotdf2)[1:5],
                  measure.vars=names(plotdf2)[6:ncol(plotdf2)],
                  variable.name = "stat",
                  value.name = "se")
  plotdf3a <- subset(plotdf3,years == "1975-1995")
  plotdf2a <- subset(plotdf2,years == "1975-1995")
  plotdf3b <- subset(plotdf3,years != "1975-1995")
  plotdf2b <- subset(plotdf2,years != "1975-1995")
  plotdf3a$se <- plotdf2a$se[match(plotdf3a$names1,plotdf2a$names1)]
  plotdf3b$se <- plotdf2b$se[match(plotdf3b$names1,plotdf2b$names1)]
  plotdf3 <- rbind(plotdf3a, plotdf3b)
}else{
  plotdf3$se <- 0 
}

plotdf3$cat1 <- dplyr::recode_factor(plotdf3$cat1, 
                                     "Observed"="Observed",
                                     "Marginal"="Marginal",
                                     "Copula"="Dependence")

# plotdf3 <- plotdf3 %>% dplyr::group_by(names1,cat1,cat2,stat) %>% 
#     dplyr::mutate(total=sum(delta),
#                   `1975-1995`=delta[which(years=="1975-1995")])  %>% as.data.frame()
plotdf3$years <- factor(plotdf3$years, levels=c("1975-1995","1995-2015"))

## Decile Ratio
plotdf4 <- plotdf3[which(plotdf3$stat=="P90/P10"),]
plotdf4 <- plotdf4 %>% dplyr::group_by(years) %>% dplyr::mutate(order1=abs(delta)) %>% as.data.frame()
plotdf4$names2 <- dplyr::recode(plotdf4$names2, 
                               "Total difference"="Observed difference",
                               "M Interactions > 0"="Marginal interactions",
                               "C Interactions > 0"="Dependence interactions")


# Plot
ggplot(plotdf4, 
       aes(reorder(names2,order1),
           delta,
           fill=cat1, 
           alpha=years)) + 
  geom_bar(stat="identity", position = "dodge", width=0.8) + 
  geom_hline(yintercept = 0, col="darkgrey") + 
  geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),
                color="black", width=0.75,
                position=position_dodge2(1.5, padding=0.1, reverse = FALSE) )+  
  coord_flip() +  
  labs(x="",y="Difference in log(D9/D1)", fill="Effect", alpha="Period") +
  scale_alpha_manual(values=c("1995-2015"=1,"1975-1995"=0.4)) 

# Save
fn <- paste0(pfolder,"summary_deco_subperiods_decile_ratio",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

# ## Gini coefficient
# plotdf4 <- plotdf3[which(plotdf3$stat=="Gini"),]
# plotdf4 <- plotdf4 %>% dplyr::group_by(years) %>% dplyr::mutate(order1=abs(delta)) %>% as.data.frame()
# plotdf4$names2 <- dplyr::recode(plotdf4$names2, 
#                                 "Total difference"="Observed difference",
#                                 "M Interactions > 0"="Marginal interactions",
#                                 "C Interactions > 0"="Dependence interactions")
# 
# 
# ggplot(plotdf4, 
#        aes(reorder(names2,order1),
#            delta,
#            fill=cat1, 
#            alpha=years)) + 
#   geom_bar(stat="identity", position = "dodge", width=0.8) + 
#   geom_hline(yintercept = 0, col="darkgrey") + 
#   geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),color="black", width=0.7,
#   #               #               position=position_dodge(0.9)) + 
#                  position=position_dodge2(1.5, padding=0.5, reverse = TRUE))+ 
#   # )+
#   coord_flip() +  
#   labs(x="",y="Difference in Gini", fill="Marginal/Copula", alpha="Years") +
#   scale_alpha_manual(values=c(0.4,1)) 
# 
# # Save
# fn <- paste0(pfolder,"summary_deco_education_gini",ptype)
# ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)


