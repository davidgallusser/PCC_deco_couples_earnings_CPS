##############################################################
## 5 - CPS - compute descriptive statistics
##############################################################

##############################################################
# Estimate descriptive statistics for all waves

# Setwd
setwd(paste0(wd,datafiles_d))

# 
stats_all_years <- NULL
stats_pos_all_years <- NULL
corrs_all_years <- NULL
years <- c(1976:2020)


# Run loop  
tic("Run descriptive stats")
for(i in years){
print(i)
load(paste0("CPScouples",i,".rda"))
stats_res <- cbind(wtd.sum(cdata[,"EARNINGS"],weight=cdata[,"ASECWTH_m"], trim=0.005, log_trans=FALSE),
                   apply(cdata[,c("WAGE_m","HOURS_m")],2,
                         function(x) wtd.sum(x,weight=cdata[,"ASECWT_m"], trim=0.005, log_trans=FALSE)),
                   apply(cdata[,c("WAGE_f","HOURS_f")],2,
                         function(x) wtd.sum(x,weight=cdata[,"ASECWT_f"], trim=0.005, log_trans=FALSE))
                   )
stats_pos_res <- cbind(apply(cdata[which(cdata$LFP_m==1),c("WAGE_m","HOURS_m")],2,
                                  function(x) wtd.sum(x,weight=cdata[which(cdata$LFP_m==1),"ASECWT_m"], trim=0.005, log_trans=FALSE)),
                       apply(cdata[which(cdata$LFP_f==1),c("WAGE_f","HOURS_f")],2,
                                  function(x) wtd.sum(x,weight=cdata[which(cdata$LFP_f==1),"ASECWT_f"], trim=0.005, log_trans=FALSE)))
corrs_res <- wtd.cor(cdata[,c("WAGE_m","HOURS_m","WAGE_f","HOURS_f")],
                  weights=cdata[,"ASECWTH_m"],trim=c(0.005,0,0.005,0))

stats_res  <- as.data.frame(stats_res)
stats_res$year <- i
stats_res$stat <- rownames(stats_res)
stats_pos_res  <- as.data.frame(stats_pos_res)
stats_pos_res$year <- i
stats_pos_res$stat <- rownames(stats_pos_res)
corrs_res$year <- i

stats_all_years <- rbind(stats_all_years,stats_res)
stats_pos_all_years <- rbind(stats_pos_all_years,stats_pos_res)
corrs_all_years <- rbind(corrs_all_years,corrs_res)

}
toc("")

# Save data 
save(stats_all_years, corrs_all_years, stats_pos_all_years, file="couple_stats_CPS_1976-2019.rda")

##############################################################
# Compare two years

# # Setwd
# setwd(paste0(wd,datafiles_d))
# years <- c(1976,2016)
# cdataall <- NULL
# for(i in years){
#   load(paste0("CPScouples",i,".rda"))
#   cdata <- cdata[sample(1:nrow(cdata),2000),]
#   cdataall <- rbind(cdataall,cdata)
# }
# cdata <- cdataall
# rm("cdataall")
# 
# # Set fomula 
# f <- EARNINGS ~ WAGE_m:HOURS_m + WAGE_f:HOURS_f
# 
# # Estimate stats of marginals
# tic("marginals")
# marginals_stats <- VC_stats_marginals(f,cdata,weights=ASECWTH_m, group=YEAR_m,
#                                       kendall = FALSE, trimgini=0.005, trimcor=c(0.005,0,0.005,0),
#                                       bootstrap = FALSE, it=100, ncore=4)
# toc("marginals")
# 
# setwd(paste0(wd,datafiles_d))
# save(marginals_stats, years, file=paste0("marginal_stats_",years[1],"-",years[2],".rda"))
