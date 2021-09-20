##############################################################
## D4 - CPS - Deco summary plots 1976-2016 - summary plots
##############################################################

years <- c(1976,2016)
setwd(paste0(wd,datafiles_d))
load(paste0("copula_decomposition_main_",years[1],"-",years[2],".rda"))
load(paste0("copula_decomposition_ref0_",years[1],"-",years[2],".rda"))
load(paste0("copula_decomposition_alt_order_",years[1],"-",years[2],".rda"))
load(paste0("copula_decomposition_parametric_",years[1],"-",years[2],".rda"))
load(paste0("copula_decomposition_cond_on_edu_",years[1],"-",years[2],".rda"))
setwd(wd)

#######################################################################
## Plot summary bar plot for main decomposition (log(D90/D10) and Gini)
#######################################################################

plotdf <- VCopObj$VCdeco$aggregate$stats_deco 
plotdf <- plotdf[which(plotdf$cat2=="direct"|plotdf$cat1=="Observed"),]

# Rename variables
plotdf$names2 <- as.character(plotdf$names2)
plotdf$names2 <- gsub("WAGE_m","Wage men", plotdf$names2)
plotdf$names2 <- gsub("WAGE_f","Wage women", plotdf$names2)
plotdf$names2 <- gsub("ANNUALHOURS_f","Hours women", plotdf$names2)
plotdf$names2 <- gsub("ANNUALHOURS_m","Hours men", plotdf$names2)
plotdf$names2 <- substr(plotdf$names2,3,nchar(plotdf$names2))
plotdf$names2 <- paste0(1:length(plotdf$names2)," ",plotdf$names2)
plotdf$names2[c(9,10,11)]<- c("9 Hours M,Wage W; Wage M", "10 Wage M,Hours W; Wage W", "11 Hours M,Hours W; Wage M,Wage W")  
plotdf$names2 <- substr(plotdf$names2,3,nchar(plotdf$names2))
plotdf$names2 <- gsub(";"," |", plotdf$names2)


plotdf <- reshape2::melt(plotdf, id.vars=names(plotdf)[1:4],
               measure.vars=names(plotdf)[5:ncol(plotdf)],
               variable.name = "stat",
               value.name = "delta")

if(is.null(VCopObj$VCdeco_se$aggregate)==FALSE){
  plotdf2 <- VCopObj$VCdeco_se$aggregate$stats_deco_se
  plotdf2 <-  plotdf2[which(plotdf2$cat2=="direct"|plotdf2$cat1=="Observed"),]
  plotdf2 <- reshape2::melt(plotdf2, id.vars=names(plotdf2)[1:4],
                 measure.vars=names(plotdf2)[5:ncol(plotdf2)],
                 variable.name = "stat",
                 value.name = "se")
  plotdf$se <- plotdf2$se[match(plotdf$names1,plotdf2$names1)]
}else{
  plotdf$se <- 0 
}


plotdf$cat1 <- dplyr::recode_factor(plotdf$cat1, 
                             "Observed"="Observed",
                             "Marginal"="Marginal",
                             "Copula"="Dependence")

plotdf$order1 <- rank(abs(plotdf[,"delta"])) # 6:1  
plotdf$names2 <- dplyr::recode(plotdf$names2, 
                                       "Total difference"="Observed difference",
                                       " M Interactions > 0"="Marginal interactions",
                                       " C Interactions > 0"="Dependence interactions")

## Decile Ratio
ggplot(plotdf[which(plotdf$stat=="P90/P10"),], aes(reorder(names2,order1),delta,fill=cat1)) + 
  geom_bar(stat="identity",  width=0.8) + #, position=position_dodge(0.9))+ 
  #geom_bar(stat="identity",  width=0.8, position=position_dodge2(0.9,reverse = TRUE)) +
  geom_hline(yintercept = 0, col="darkgrey") + 
  geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),color="black", width=0.7,
                #               position=position_dodge(0.9)) + 
                # position=position_dodge2(1.5, padding=0.5, reverse = TRUE))+ 
                )+
  coord_flip() +  labs(x="",y="Difference in log(D9/D1)", fill="") 

# Save
fn <- paste0(pfolder,"summary_deco_main",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

# ## Gini
# ggplot(plotdf[which(plotdf$stat=="Gini"),], aes(reorder(names2,order1),delta,fill=cat1)) + 
#   geom_bar(stat="identity",  width=0.8) + #, position=position_dodge(0.9))+ 
#   #geom_bar(stat="identity",  width=0.8, position=position_dodge2(0.9,reverse = TRUE)) +
#   geom_hline(yintercept = 0, col="darkgrey") + 
#   geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),color="black", width=0.7,
#                 #               position=position_dodge(0.9)) + 
#                 # position=position_dodge2(1.5, padding=0.5, reverse = TRUE))+ 
#   )+
#   coord_flip() +  labs(x="",y="Difference in Gini coefficient", fill="") 
# 
# # Save
# fn <- paste0(pfolder,"summary_deco_main_gini",ptype)
# ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)
# 

#######################################################################
## Plot summary bar plot for decomposition with switched reference group (log(D90/D10))
#######################################################################

plotdf3 <-  VCopObj_ref0$VCdeco$aggregate$stats_deco
plotdf3 <-  plotdf3[which(plotdf3$cat2=="direct"|plotdf3$cat1=="Observed"),]

# Rename variables
plotdf3$names2 <- as.character(plotdf3$names2)
plotdf3$names2 <- gsub("WAGE_m","Wage men", plotdf3$names2)
plotdf3$names2 <- gsub("WAGE_f","Wage women", plotdf3$names2)
plotdf3$names2 <- gsub("ANNUALHOURS_f","Hours women", plotdf3$names2)
plotdf3$names2 <- gsub("ANNUALHOURS_m","Hours men", plotdf3$names2)
plotdf3$names2 <- substr(plotdf3$names2,3,nchar(plotdf3$names2))
plotdf3$names2 <- paste0(1:length(plotdf3$names2)," ",plotdf3$names2)
plotdf3$names2[c(9,10,11)]<- c("9 Hours M,Wage W; Wage M", "10 Wage M,Hours W; Wage W", "11 Hours M,Hours W; Wage M,Wage W")  
plotdf3$names2 <- substr(plotdf3$names2,3,nchar(plotdf3$names2))
plotdf3$names2 <- gsub(";"," |", plotdf3$names2)

plotdf3 <- melt(plotdf3, id.vars=names(plotdf3)[1:4],
               measure.vars=names(plotdf3)[5:ncol(plotdf3)],
               variable.name = "stat",
               value.name = "delta")

if(is.null(VCopObj_ref0$VCdeco_se$aggregate)==FALSE){
  plotdf2 <- VCopObj_ref0$VCdeco_se$aggregate$stats_deco_se
  plotdf2 <-  plotdf2[which(plotdf2$cat2=="direct"|plotdf2$cat1=="Observed"),]
  plotdf2 <- melt(plotdf2, id.vars=names(plotdf2)[1:4],
                  measure.vars=names(plotdf2)[5:ncol(plotdf2)],
                  variable.name = "stat",
                  value.name = "se")
  plotdf3$se <- plotdf2$se[match(plotdf3$names1,plotdf2$names1)]
}else{
  plotdf3$se <- 0 
}

plotdf3$cat1 <- dplyr::recode_factor(plotdf3$cat1, 
                                    "Observed"="Observed",
                                    "Marginal"="Marginal",
                                    "Copula"="Dependence")

plotdf3$order1 <- rank(abs(plotdf3[,"delta"])) # 6:1  
plotdf3$names2 <- dplyr::recode(plotdf3$names2, 
                                "Total difference"="Observed difference",
                               " M Interactions > 0"="Marginal interactions",
                               " C Interactions > 0"="Dependence interactions")

## Decile Ratio
ggplot(plotdf3[which(plotdf3$stat=="P90/P10"),], aes(reorder(names2,order1),delta,fill=cat1)) + 
  geom_bar(stat="identity",  width=0.8) + #, position=position_dodge(0.9))+ 
  #geom_bar(stat="identity",  width=0.8, position=position_dodge2(0.9,reverse = TRUE)) +
  geom_hline(yintercept = 0, col="darkgrey") + 
  geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),color="black", width=0.7,
                #               position=position_dodge(0.9)) + 
                # position=position_dodge2(1.5, padding=0.5, reverse = TRUE))+ 
  )+
  coord_flip() +  labs(x="",y="Difference in log(D9/D1)", fill="") 

# Save
fn <- paste0(pfolder,"summary_deco_switched_ref",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

### Compare decomposition with different reference groups
plotdf$reference <- "2015 (baseline)"
plotdf3$reference <- "1975"
plotdf4 <- rbind(plotdf,plotdf3)

## Decile Ratio
ggplot(plotdf4[which(plotdf4$stat=="P90/P10"&plotdf4$cat1!="Total"),], aes(reorder(names2,order1),delta,fill=reference)) + 
  #geom_bar(stat="identity",  width=0.8) + #, position=position_dodge(0.9))+ 
  geom_bar(stat="identity",  width=0.8, position=position_dodge2(0.9,reverse = TRUE)) +
  geom_hline(yintercept = 0, col="darkgrey") + 
  geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),color="black", width=0.7,
                        #     position=position_dodge(0.9)
                 position=position_dodge2(1.5, padding=0.2, reverse = TRUE)
  )+
  coord_flip() +  labs(x="",y="Difference in log(D9/D1)", fill="Reference copula") 

# Save
fn <- paste0(pfolder,"summary_deco_compare_reference_groups",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)


#######################################################################
## Summary bar plot for decomposition based on alternative DVine order (1,2,3,4)
#######################################################################

plotdf3 <- VCopObj_alt_order$VCdeco$aggregate$stats_deco
plotdf3 <-  plotdf3[which(plotdf3$cat2=="direct"|plotdf3$cat1=="Observed"),]

# Rename variables
plotdf3$names2 <- as.character(plotdf3$names2)
plotdf3$names2 <- gsub("WAGE_m","Wage men", plotdf3$names2)
plotdf3$names2 <- gsub("WAGE_f","Wage women", plotdf3$names2)
plotdf3$names2 <- gsub("ANNUALHOURS_f","Hours women", plotdf3$names2)
plotdf3$names2 <- gsub("ANNUALHOURS_m","Hours men", plotdf3$names2)
plotdf3$names2 <- substr(plotdf3$names2,3,nchar(plotdf3$names2))
plotdf3$names2 <- paste0(1:length(plotdf3$names2)," ",plotdf3$names2)
plotdf3$names2[c(9,10,11)]<-  c("9 Wage M,Wage W; Hours M", "10 Hours M,Hours W; Wage W", "11 Wage M,Hours W; Hours M,Wage W")  
plotdf3$names2 <- substr(plotdf3$names2,3,nchar(plotdf3$names2))
plotdf3$names2 <- gsub(";"," |", plotdf3$names2)

plotdf3 <- reshape2::melt(plotdf3, id.vars=names(plotdf3)[1:4],
                measure.vars=names(plotdf3)[5:ncol(plotdf3)],
                variable.name = "stat",
                value.name = "delta")

if(is.null(VCopObj_alt_order$VCdeco_se$aggregate)==FALSE){
  plotdf2 <- VCopObj_alt_order$VCdeco_se$aggregate$stats_deco_se
  plotdf2 <-  plotdf2[which(plotdf2$cat2=="direct"|plotdf2$cat1=="Observed"),]
  plotdf2 <- reshape2::melt(plotdf2, id.vars=names(plotdf2)[1:4],
                  measure.vars=names(plotdf2)[5:ncol(plotdf2)],
                  variable.name = "stat",
                  value.name = "se")
  plotdf3$se <- plotdf2$se[match(plotdf3$names1,plotdf2$names1)]
}else{
  plotdf3$se <- 0 
}

plotdf3$cat1 <- dplyr::recode_factor(plotdf3$cat1, 
                                     "Observed"="Observed",
                                     "Marginal"="Marginal",
                                     "Copula"="Dependence")

plotdf3$order1 <- rank(abs(plotdf3[,"delta"])) # 6:1  
plotdf3$names2 <- dplyr::recode(plotdf3$names2, 
                                "Total difference"="Observed difference",
                               " M Interactions > 0"="Marginal interactions",
                               " C Interactions > 0"="Dependence interactions")

## Decile Ratio
ggplot(plotdf3[which(plotdf3$stat=="P90/P10"),], aes(reorder(names2,order1),delta,fill=cat1)) + 
  geom_bar(stat="identity",  width=0.8) + #, position=position_dodge(0.9))+ 
  #geom_bar(stat="identity",  width=0.8, position=position_dodge2(0.9,reverse = TRUE)) +
  geom_hline(yintercept = 0, col="darkgrey") + 
  geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),color="black", width=0.7,
                #               position=position_dodge(0.9)) + 
                # position=position_dodge2(1.5, padding=0.5, reverse = TRUE))+ 
  )+
  coord_flip() +  labs(x="",y="Difference in log(D9/D1)", fill="") 

# Save
fn <- paste0(pfolder,"summary_deco_alt_order",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

### Compare decomposition with different reference groups
plotdf$reference <- "D-vine 2-1-3-4\n(baseline)"
plotdf3$reference <- "D-vine 1-2-3-4"

n1 <- setdiff(unique(plotdf$names2),unique(plotdf3$names2))
n2 <- setdiff(unique(plotdf3$names2),unique(plotdf$names2))
n1 <- n1[c(2,1,4,3)]
cbind(n1,n2)
n3 <- c(n2[1],n1[2],"Hours men,Hours women", "Wage men,Hours women")
cbind(n1,n2,n3)
cbind(n1,n3)
cbind(n2,n3)

plotdf1 <- plotdf
plotdf2 <- plotdf3

plotdf1$condcopula <- ifelse(plotdf1$names2 %in% n1[c(1,3,4)],"","Conditional copula")
plotdf2$condcopula <- ifelse(plotdf2$names2 %in% n2[2:4],"","Conditional copula")

plotdf1$names2 <- dplyr::recode(plotdf1$names2,
                                "Hours M,Wage W | Wage M"= "Hours men,Wage women", 
                                "Wage men,Wage women"="Wage men,Wage women" , 
                                 " Hours M,Hours W | Wage M,Wage W"="Hours men,Hours women",
                                " Wage M,Hours W | Wage W"="Wage men,Hours women" )
plotdf2$names2 <- dplyr::recode(plotdf2$names2,
                                "Hours men,Wage women"="Hours men,Wage women",
                                "Wage M,Wage W | Hours M"="Wage men,Wage women",  
                                " Hours M,Hours W | Wage W"="Hours men,Hours women",
                                " Wage M,Hours W | Hours M,Wage W"="Wage men,Hours women")

plotdf4 <- rbind(plotdf1,plotdf2)

## Decile Ratio
ggplot(plotdf4[which(plotdf4$stat=="P90/P10"&plotdf4$cat1!="Total"),], 
       aes(reorder(names2,order1),delta,
           fill=reference,
           alpha=condcopula)) + 
  #geom_bar(stat="identity",  width=0.8) + #, position=position_dodge(0.9))+ 
  geom_bar(stat="identity",  width=0.8, position=position_dodge2(0.9,reverse = TRUE)) +
  geom_hline(yintercept = 0, col="darkgrey") + 
  geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),color="black", width=0.7,
                #               position=position_dodge(0.9)
                position=position_dodge2(1.5, padding=0.5, reverse = TRUE)
  )+
  scale_alpha_discrete(range=c(0.5,1), guide = FALSE)+ 
  coord_flip() +  labs(x="",y="Difference in log(D9/D1)", fill="Order") 

# Save
fn <- paste0(pfolder,"summary_deco_compare_with_alt_order",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)


#######################################################################
## Plot summary bar plot for parametrically estimated copulas (log(D90/D10))
#######################################################################

plotdf3 <-  VCopObj_parest$VCdeco$aggregate$stats_deco
plotdf3 <-  plotdf3[which(plotdf3$cat2=="direct"|plotdf3$cat1=="Observed"),]

# Rename variables
plotdf3$names2 <- as.character(plotdf3$names2)
plotdf3$names2 <- gsub("WAGE_m","Wage men", plotdf3$names2)
plotdf3$names2 <- gsub("WAGE_f","Wage women", plotdf3$names2)
plotdf3$names2 <- gsub("ANNUALHOURS_f","Hours women", plotdf3$names2)
plotdf3$names2 <- gsub("ANNUALHOURS_m","Hours men", plotdf3$names2)
plotdf3$names2 <- substr(plotdf3$names2,3,nchar(plotdf3$names2))
plotdf3$names2 <- paste0(1:length(plotdf3$names2)," ",plotdf3$names2)
plotdf3$names2[c(9,10,11)]<- c("9 Hours M,Wage W; Wage M", "10 Wage M,Hours W; Wage W", "11 Hours M,Hours W; Wage M,Wage W")  
plotdf3$names2 <- substr(plotdf3$names2,3,nchar(plotdf3$names2))
plotdf3$names2 <- gsub(";"," |", plotdf3$names2)

plotdf3 <- melt(plotdf3, id.vars=names(plotdf3)[1:4],
                measure.vars=names(plotdf3)[5:ncol(plotdf3)],
                variable.name = "stat",
                value.name = "delta")

if(is.null(VCopObj_ref0$VCdeco_se$aggregate)==FALSE){
  plotdf2 <- VCopObj_ref0$VCdeco_se$aggregate$stats_deco_se
  plotdf2 <-  plotdf2[which(plotdf2$cat2=="direct"|plotdf2$cat1=="Observed"),]
  plotdf2 <- melt(plotdf2, id.vars=names(plotdf2)[1:4],
                  measure.vars=names(plotdf2)[5:ncol(plotdf2)],
                  variable.name = "stat",
                  value.name = "se")
  plotdf3$se <- plotdf2$se[match(plotdf3$names1,plotdf2$names1)]
}else{
  plotdf3$se <- 0 
}

plotdf3$cat1 <- dplyr::recode_factor(plotdf3$cat1, 
                                     "Observed"="Observed",
                                     "Marginal"="Marginal",
                                     "Copula"="Dependence")

plotdf3$order1 <- rank(abs(plotdf3[,"delta"])) # 6:1  
plotdf3$names2 <- dplyr::recode(plotdf3$names2, 
                                "Total difference"="Observed difference",
                               " M Interactions > 0"="Marginal interactions",
                               " C Interactions > 0"="Dependence interactions")

## Decile Ratio
ggplot(plotdf3[which(plotdf3$stat=="P90/P10"),], aes(reorder(names2,order1),delta,fill=cat1)) + 
  geom_bar(stat="identity",  width=0.8) + #, position=position_dodge(0.9))+ 
  #geom_bar(stat="identity",  width=0.8, position=position_dodge2(0.9,reverse = TRUE)) +
  geom_hline(yintercept = 0, col="darkgrey") + 
  geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),color="black", width=0.7,
                #               position=position_dodge(0.9)) + 
                # position=position_dodge2(1.5, padding=0.5, reverse = TRUE))+ 
  )+
  coord_flip() +  labs(x="",y="Difference in log(D9/D1)", fill="") 

# Save
fn <- paste0(pfolder,"summary_deco_parametric",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

### Compare decomposition with different estimation methods
plotdf$reference <- "non-parametric\n(baseline)"
plotdf3$reference <- "parametric"
plotdf4 <- rbind(plotdf,plotdf3)

## Decile Ratio
ggplot(plotdf4[which(plotdf4$stat=="P90/P10"&plotdf4$cat1!="Total"),], aes(reorder(names2,order1),delta,fill=reference)) + 
  #geom_bar(stat="identity",  width=0.8) + #, position=position_dodge(0.9))+ 
  geom_bar(stat="identity",  width=0.8, position=position_dodge2(0.9,reverse = TRUE)) +
  geom_hline(yintercept = 0, col="darkgrey") + 
  geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),color="black", width=0.7,
                #     position=position_dodge(0.9)
                position=position_dodge2(1.5, padding=0.2, reverse = TRUE)
  )+
  coord_flip() +  labs(x="",y="Difference in log(D9/D1)", fill="Estimation") 

# Save
fn <- paste0(pfolder,"summary_deco_compare_parametric_non_parametric",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)


#######################################################################
## Summary bar plot for decomposition conditional on education & potex
#######################################################################

plotdf3 <- VCopObj_cond$VCdeco$composition_effect$stats_deco
plotdf4 <- VCopObj_cond$VCdeco$wage_structure_effect$stats_deco
plotdf3$effect <- "Education"
plotdf4$effect <- "Other"

plotdf3 <- rbind(plotdf3,plotdf4)
plotdf3 <-  subset(plotdf3,cat2=="direct"|cat1=="Observed")

# Rename variables
plotdf3$names2 <- as.character(plotdf3$names2)
plotdf3$names2 <- gsub("WAGE_m","Wage Men", plotdf3$names2)
plotdf3$names2 <- gsub("WAGE_f","Wage Women", plotdf3$names2)
plotdf3$names2 <- gsub("ANNUALHOURS_f","Hours Women", plotdf3$names2)
plotdf3$names2 <- gsub("ANNUALHOURS_m","Hours Men", plotdf3$names2)
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

if(is.null(VCopObj_cond$VCdeco_se$aggregate)==FALSE){
  plotdf2a <- VCopObj_cond$VCdeco_se$composition_effect$stats_deco
  plotdf2b <- VCopObj_cond$VCdeco_se$wage_structure_effect$stats_deco
  plotdf2a$effect <- "Education"
  plotdf2b$effect <- "Other"
  plotdf2 <- rbind(plotdf2a,plotdf2b)
  plotdf2 <-  subset(plotdf2,cat2=="direct"|cat1=="Observed")
  plotdf2 <- plotdf2[, c(1:4,ncol(plotdf2),5:(ncol(plotdf2)-1))]
  plotdf2 <- melt(plotdf2, id.vars=names(plotdf2)[1:5],
                  measure.vars=names(plotdf2)[6:ncol(plotdf2)],
                  variable.name = "stat",
                  value.name = "se")
  plotdf3a <- subset(plotdf3,effect =="Education")
  plotdf2a <- subset(plotdf2,effect =="Education")
  plotdf3b <- subset(plotdf3,effect !="Education")
  plotdf2b <- subset(plotdf2,effect !="Education")
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

plotdf3 <- plotdf3 %>% dplyr::group_by(names1,cat1,cat2,stat) %>% 
            dplyr::mutate(total=sum(delta),
                          composition=delta[which(effect=="Education")])  %>% as.data.frame()
plotdf3$effect <- factor(plotdf3$effect, levels=c("Other","Education"))
plotdf3$effect <- dplyr::recode_factor(plotdf3$effect, 
                          "Education"="Edu. & Exper.\n(Composition)",
                          "Other"="Other\n(Structure)")
plotdf3$names2 <- dplyr::recode(plotdf3$names2, 
                                "Total difference"="Observed difference",
                                "M Interactions > 0"="Marginal interactions",
                                "C Interactions > 0"="Dependence interactions")

## Decile Ratio
plotdf4 <- plotdf3[which(plotdf3$stat=="P90/P10"),]
plotdf4$order1 <- rank(abs(plotdf4[,"total"])) # 6:1  
plotdf4$order1 <- rank(abs(plotdf4[,"composition"])) # 6:1 


ggplot(plotdf4, 
       aes(reorder(names2,order1),
           delta,
           fill=cat1, 
           alpha=effect)) + 
  geom_bar(stat="identity", 
           #position = "stack", width=0.8) + 
           width=0.95, position=position_dodge2(0.9,reverse = TRUE)) +
  geom_hline(yintercept = 0, col="darkgrey") + 
  geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),
                    color="black", width=0.75,
                              position=position_dodge2(1.5, padding=0.1, reverse = TRUE)
  )+
  coord_flip() +  
  labs(x="",y="Difference in log(D9/D1)", fill="Copula deco", alpha="Observables") +
  scale_alpha_manual(values=c(1,0.4)) 

# Save
fn <- paste0(pfolder,"summary_deco_education_decile_ratio",ptype)
ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)

# ## Gini coefficient
# plotdf4 <- plotdf3[which(plotdf3$stat=="Gini"),]
# plotdf4$order1 <- rank(abs(plotdf4[,"total"])) # 6:1  
# plotdf4$order1 <- rank(abs(plotdf4[,"composition"])) # 6:1 
# 
# ggplot(plotdf4, 
#        aes(reorder(names2,order1),
#            delta,
#            fill=cat1, 
#            alpha=effect)) + 
#   geom_bar(stat="identity", position = "stack", width=0.8) + 
#   geom_hline(yintercept = 0, col="darkgrey") + 
#   # geom_errorbar(aes(ymin=delta-1.96*se, ymax=delta+1.96*se),color="black", width=0.7,
#   #               #               position=position_dodge(0.9)) + 
#   #               # position=position_dodge2(1.5, padding=0.5, reverse = TRUE))+ 
#   # )+
#   coord_flip() +  
#   labs(x="",y="Difference in log(D9/D1)", fill="Copula deco", alpha="Observables") +
#   scale_alpha_manual(values=c(1,0.4)) 
# 
# # Save
# fn <- paste0(pfolder,"summary_deco_education_gini",ptype)
# ggsave(fn, units=punits, scale=pscale, width=pwidth, height=pheight)
