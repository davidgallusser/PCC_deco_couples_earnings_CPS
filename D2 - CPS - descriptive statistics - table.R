##############################################################
## D2 - CPS - Descriptive statistics tables
##############################################################

##############################################################
# Load data
setwd(paste0(wd,datafiles_d))
load("couple_stats_CPS_1976-2019.rda")
setwd(wd)

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
                                        "Weekly hours" <- c("HOURS_m","HOURS_f"),    
                                        "Hourly wages" <- c("WAGE_m","WAGE_f"),
                                        "Earnings" <- c("EARNINGS"))                                           
stats_all_years2$year <- stats_all_years2$year-1 

##############################################################
# Create table
####              %obs.=0, P10, P50, P90, Gini
#### Earnings 
     # 1976
     # 2016
     # % change 
#### Wages Men
#### Wages Women
#### Hours Men
#### Hours Women

sely <- c(1976,2016)
selstat <- c("%obs.=0")
tab1 <- subset(stats_all_years2, 
               year %in% sely & stat %in% selstat & cat %in%  "Employed & non-active")

selstat <- c("%obs.=0","P10","P50","P90", "Gini")
tab2 <- subset(stats_all_years2, 
               year %in% sely & stat %in% selstat & var2 %in%  "Earnings")

selstat <- c("P10","P50","P90", "Gini")
tab3 <- subset(stats_all_years2, 
               year %in% sely & stat %in% selstat & cat %in%  "Employed")

tab <- rbind(tab1, tab2, tab3)
tab <- tab[,-c(3,5)]

# Replace Unemployed Rate in Wage Variable
tab[which(tab$gender=="Men"&tab$stat=="%obs.=0"&tab$var2=="Hourly wages"), "value"] <- tab[which(tab$gender=="Men"&tab$stat=="%obs.=0"&tab$var2=="Weekly hours"), "value"]
tab[which(tab$gender=="Women"&tab$stat=="%obs.=0"&tab$var2=="Hourly wages"), "value"] <- tab[which(tab$gender=="Women"&tab$stat=="%obs.=0"&tab$var2=="Weekly hours"), "value"]

# Adjust year
tab$year <- tab$year-1

# Add growth rate
add <- tab  %>% dplyr::group_by(stat,var2,gender) %>% dplyr::summarise(year="Growth 1975-2015",
                                                                       value=(value[2]-value[1])/value[1]*100) %>% as.data.frame()
tab <- rbind(tab,add[,names(tab)])
tab$var <- paste0(tab$gender," ", tab$var2)
tab <- tab[,-c(4,5)]
tab <- dcast(tab, var + year ~ stat, mean)

# Put into right format
sely <- sely -1 
tab[which(tab$year %in% sely == FALSE),3:ncol(tab)] <- apply(tab[which(tab$year %in% sely==FALSE),3:ncol(tab)],2,scales::comma,suffix="\\%",accuracy=1)
tab[which(tab$year %in% sely),c("%obs.=0","Gini")] <- apply(tab[which(tab$year %in% sely),c("%obs.=0","Gini")],2,function(x) scales::comma(as.numeric(x), accuracy = 0.001))
sel <- which(tab$year %in% sely & tab$var %in% c("Men Weekly hours","Women Weekly hours"))
tab[sel,5:7] <- apply(tab[sel,5:7],2,function(x) scales::comma(as.numeric(x), accuracy = 0.1))
sel <- which(tab$year %in% sely & tab$var %in% c("Men Hourly wages","Women Hourly wages"))
tab[sel,5:7] <- apply(tab[sel,5:7],2,function(x) scales::comma(as.numeric(x), accuracy = 0.01, prefix="\\$"))
sel <- which(tab$year %in% sely & tab$var %in% c("NA Earnings"))
tab[sel,5:7] <- apply(tab[sel,5:7],2,function(x) scales::comma(as.numeric(x), accuracy = 1, prefix="\\$", big.mark = ","))


  
tab <- tab[,c(1,2,3,5,6,7,4)]
tab$var <- dplyr::recode_factor(tab$var,
                         "NA Earnings"="Couples' annual earnings",     
                         "Men Hourly wages"="Male hourly wages",
                         "Women Hourly wages"="Female hourly wages",
                         "Men Weekly hours"="Male weekly hours", 
                         "Women Weekly hours"="Female weekly hours")
tab$`%obs.=0` <- tidyr::replace_na(tab$`%obs.=0`, 0)

tab$year <- paste0("\\quad ",tab$year)
tab$year <- factor(tab$year, levels=unique(tab$year)[c(1,2,3)])
tab <- tab %>% arrange(var,year)

uvars<- unique(tab$var)
tab <- tab[,-1]
# add <- tab[1:length(uvars),]
# add$year <- uvars
# add[,-1] <- NA

## 
# tab <- rbind(add,tab)
# tab <- tab[c(1,6:8,
#        2,9:11,
#        3,12:14,
#        4,15:17,
#        5,18:20),]

## Title
fn <- paste("summary_marginals.tex")
title <- "Descriptive statistics of couples' earnings, wages, and hours \\label{tab:desc}"
names(tab) <- c("","Share $X=0$","1st decile","Median","9th decile", "Gini")
names(tab) <- paste("\\multicolumn{1}{c}{",names(tab), "}")
#add_first_row <- "\\toprule\n & \\multicolumn{1}{c}{\\textbf{All}} &\\multicolumn{4}{c}{\\textbf{Employed}}  \\\\  \n "


add_first_row <- "\\toprule\n & \\multicolumn{1}{c}{\\textbf{All}}   &
                          \\multicolumn{4}{c}{\\textbf{Employed}}   \\\\ \\cmidrule(r){2-2}\\cmidrule{3-6} \n "
add_first_row <- paste0(add_first_row,paste0(names(tab), collapse="&"),"\\\\\n",collapse="")
add_rows <- paste0("\\multicolumn{6}{l}{\\textbf{", uvars,"}} \\\\ \n")
add_source <- "\\bottomrule \n 
               \\multicolumn{6}{l}{\\footnotesize{Source: March CPS 1976 and 2016, author's calculations.}}"

add_rows <-c(add_first_row,add_rows,add_source)


# Digits
tab.xtab <- xtable(tab, 
                   caption=title, 
                   type="latex",
                   booktabs = FALSE,
                   auto=TRUE)
digits(tab.xtab) <- matrix(c(rep(rep(c(3,3,0),5),3),
                             rep(c(rep(c(0,0,0),1),rep(c(2,2,0),2),rep(c(2,2,0),2)),3),
                             rep(c(3,3,0),5)), ncol=7, byrow=FALSE)
align(tab.xtab) <- c("X","X",rep("R{1.8cm}",5))


print(tab)

# Save table
setwd(pfolder)
mod_xtable(fn,tab.xtab,
           booktabs = TRUE, 
           caption.placement = "top",
           include.rownames = FALSE, include.colnames=FALSE,
           format.args=list(big.mark = ",",  decimal.mark = "."),
           tabular.environment = 'tabularx', width="0.95\\textwidth",
           sanitize.text.function = function(x){x},
           size = "small",
           #sanitize.colnames.function = identity, 
           #sanitize.rownames.function = identity,
           hline.after=c(-1),
           add.to.row = list(pos = list(-1,0,3,6,9,12,15),
                             command = add_rows)
           )
setwd(wd)


