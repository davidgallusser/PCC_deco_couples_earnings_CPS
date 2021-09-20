##############################################################
## A4 CPS wage imputation
##############################################################

setwd(paste0(wd,datafiles_d))

##############################################################
## Models

# Selection model
# NCHILD       Number of own children in household  
# NCHLT5       Number of own children in household under age 5 
# INCUNERN     Unearned income 
# INIDR        Income from interests, dividends, rent (to be constructed for years after 1976)

# INCWAGE      Labor earnings (of spouse)
# NCHILD518      (constructed) NCHLTleq5 = NCHILD -  NCHLT5

# Outcome model
# METRO        METRO Typ: 0 not identifiable, 1 not in metro, 2 central city, 3 outside central city 4 central city status missing
# STATEIP      State
# REGION       Region
# RACE         1=white, 2=black; there are other categories but they are not consistent over time
# HISPAN     
# EDUC

sf <- LFP_f ~ NCHLT5_f + NCHILD518_f + asinh(abs(INCIDR_hh)) + asinh(INCWAGE_m) + EDUC_f + AGE_f + I(AGE_f^2) + METRO_f  + REGION_f
wf <- log(WAGE_f) ~ EDUC_f*(EXP_f + I(EXP_f^2)) + METRO_f + REGION_f + RACE_f + HISPAN_f 

sm <- LFP_m ~ NCHLT5_m + NCHILD518_m + asinh(abs(INCIDR_hh)) + asinh(INCWAGE_f) + EDUC_m + AGE_m + I(AGE_m^2) + METRO_m  + REGION_m
wm <- log(WAGE_m) ~ EDUC_m*(EXP_m + I(EXP_m^2)) + METRO_m + REGION_m + RACE_m + HISPAN_m 

# selection <- sf
# outcome <- wf
# #data=cdata
# weights="ASECWT_f" 
# group="YEAR_f"
# method="ml"
# short_info=TRUE

# Hours imputation model
# https://www.econstor.eu/bitstream/10419/20270/1/dp1035.pdf
# sf <- LFP_f ~ NCHLT5_f + NCHILD518_f + asinh(abs(INCIDR_hh)) + asinh(INCWAGE_m) + EDUC_f + AGE_f + I(AGE_f^2) + METRO_f  + REGION_f
# hf <- HOURS_f ~ NCHLT5_f + NCHILD518_f + log(WAGE_f) + EDUC_m




##############################################################
## Impute wages: Loop

years <- 1976:2020

tic("Impute all")
for(i in years){
  print(i)
  
  #Load data
  load(paste0("CPScouples",i,".rda"))
  
  if(i!=1993){
    sm <- LFP_m ~ NCHLT5_m + NCHILD518_m + asinh(abs(INCIDR_hh)) + asinh(INCWAGE_f) + EDUC_m + AGE_m + I(AGE_m^2) + METRO_m  + REGION_m
  }else{
    sm <- LFP_m ~ NCHLT5_m + NCHILD518_m + asinh(abs(INCIDR_hh)) + asinh(INCWAGE_f) + EDUC2_m + AGE_m + I(AGE_m^2) + METRO_m  + REGION_m
  }
  
  # Impute wages
  cat("--------------------------------------------------\n")
  cat("Men\n")
  set.seed(123)
  cdata <- impute_wages(sm,
                        wm,
                        data=cdata, weights = "ASECWT_m", group="YEAR_m", method="ml",
                        short_info=TRUE)
  cat("\nWomen\n")
  set.seed(123)
  cdata <- impute_wages(sf,
                        wf,
                        data=cdata, weights = "ASECWT_f", group="YEAR_f", method="ml",
                        short_info=TRUE)
  
  # Create earnings variable
  cdata$EARNINGS <- (cdata$WAGE_m*cdata$ANNUALHOURS_m + cdata$WAGE_f*cdata$ANNUALHOURS_f)
  
  # Save data
  save(cdata, file=paste0("CPScouples",i,".rda"))
}
toc()

##############################################################
## Heckman model in main years 1975 and 2015

i <- 1976
load(paste0("CPScouples",i,".rda"))
sel_men_1976 <- selection(sm,wm,cdata,weights=cdata$ASECWT_m, method="ml")
sel_women_1976 <- selection(sf,wf,cdata,weights=cdata$ASECWT_f, method="ml")

i <- 2016
load(paste0("CPScouples",i,".rda"))
sel_men_2016 <- selection(sm,wm,cdata,weights=cdata$ASECWT_m, method="ml")
sel_women_2016 <- selection(sf,wf,cdata,weights=cdata$ASECWT_f, method="ml")

# Create Output Table
tab1 <- cbind(summary(sel_men_1976)$estimate[,c(1,2)],
              summary(sel_women_1976)$estimate[,c(1,2)])
tab2 <- cbind(summary(sel_men_2016)$estimate[,c(1,2)],
              summary(sel_women_2016)$estimate[,c(1,2)])
tab <-  cbind(tab1[match(rownames(tab2),rownames(tab1)),], tab2)
rownames(tab) <- gsub("_m","",rownames(tab2))
rownames(tab)[4] <- "asinh(abs(INCIDR hh))" 
rownames(tab)[5] <- "asinh(INCWAGE spouse)" 
rownames(tab)[14] <- "I(AGE2)" 
rownames(tab) <- gsub("EXP\\^2","EXP2",rownames(tab)) 
rownames(tab)[1:65] <- paste0("\\texttt{",rownames(tab)[1:65],"}")
rownames(tab)[66:67] <- c("$\\hat{\\sigma}$","$\\hat{\\rho}$")
tab1 <- format(round(tab[,c(1,3,5,7)],3),scientific=FALSE)
tab2 <- format(round(tab[,c(2,4,6,8)],3),scientific=FALSE)
tab <- cbind(tab1[,1],tab2[,1],tab1[,2],tab2[,2],tab1[,3],tab2[,3],tab1[,4],tab2[,4])
tab <- gsub(" ","",tab)
tab[,c(2,4,6,8)] <- paste0("(",tab[,c(2,4,6,8)],")")
tab_add <- cbind(c(#as.numeric(summary(sel_men_1976)$loglik),
                unlist(summary(sel_men_1976)$param[c("nObs","N0","N1")])),
                 rep(NA,3),
                 c(#as.numeric(summary(sel_women_1976)$loglik),
                  unlist(summary(sel_women_1976)$param[c("nObs","N0","N1")])),
                 rep(NA,3),
                 c(#as.numeric(summary(sel_men_2016)$loglik),
                   unlist(summary(sel_men_2016)$param[c("nObs","N0","N1")])),
                 rep(NA,3),
                 c(#as.numeric(summary(sel_women_2016)$loglik),
                   unlist(summary(sel_women_2016)$param[c("nObs","N0","N1")])),
                 rep(NA,3))

rownames(tab_add) <- c(#"Log-Likelihood",
                       "Observations",
                       "Censored",
                       "Observed")
tab_add <- round(tab_add,0)
tab <- rbind(tab,tab_add)

tab <- gsub("(NA)","",tab)
tab <- gsub("NA","",tab)
tab <- gsub("\\(\\)","",tab)
tab <- cbind(rownames(tab),tab)
colnames(tab) <- c("",rep(c("Men","","Women",""),2))

## Title
fn <- paste("Heckman_estimates_1975-2015-selection.tex")
title <- "Estimates of Heckman model: selection equation (s.e. in parentheses)\\label{tab:heckman1}"
add_first_row <- "\\toprule\n & \\multicolumn{4}{c}{\\textbf{1975}}   &
                          \\multicolumn{4}{c}{\\textbf{2015}}   \\\\ \\cmidrule(r){2-5}\\cmidrule{5-9} \n 
                              &  \\multicolumn{2}{c}{Men} & \\multicolumn{2}{c}{Women} &  \\multicolumn{2}{c}{Men} & \\multicolumn{2}{c}{Women}  \\\\\n" 

add_rows <- c("\\multicolumn{9}{l}{\\textbf{Selection equation}} \\\\ \n")
add_source <- "\\bottomrule \n 
               \\multicolumn{9}{l}{\\footnotesize{Source: March CPS 1976 and 2016, author's calculations.}}"

add_rows <-c(add_first_row,add_rows,add_source)

# Digits
tab.xtab <- xtable(tab[1:26,], 
                   caption=title, 
                   type="latex",
                   booktabs = FALSE,
                   auto=TRUE)
align(tab.xtab) <- c("X","X",rep(c("R{1.1cm}","R{1cm}"),4))

# Save table
setwd(pfolder)
mod_xtable(fn,tab.xtab,
           booktabs = TRUE, 
           caption.placement = "top",
           include.rownames = FALSE, include.colnames=FALSE,
           format.args=list(big.mark = ",",  decimal.mark = "."),
           tabular.environment = 'tabularx', width="1\\textwidth",
           sanitize.text.function = function(x){x},
           size = "footnotesize",
           #sanitize.colnames.function = identity, 
           sanitize.rownames.function = identity,
           hline.after=c(-1),
           add.to.row = list(pos = list(-1,0,26),
                             command = add_rows)
)
setwd(wd)


## Title
fn <- paste("Heckman_estimates_1975-2015-outcome.tex")
title <- "Estimates of Heckman model: outcome equation (s.e. in parentheses)\\label{tab:heckman2}"
add_first_row <- "\\toprule\n & \\multicolumn{4}{c}{\\textbf{1975}}   &
                          \\multicolumn{4}{c}{\\textbf{2015}}   \\\\ \\cmidrule(r){2-5}\\cmidrule{5-9} \n 
                              &  \\multicolumn{2}{c}{Men} & \\multicolumn{2}{c}{Women} &  \\multicolumn{2}{c}{Men} & \\multicolumn{2}{c}{Women}  \\\\\n" 

add_rows <- c("\\multicolumn{9}{l}{\\textbf{Outcome equation}} \\\\ \n")
add_source <- "\\bottomrule \n 
               \\multicolumn{9}{l}{\\footnotesize{Source: March CPS 1976 and 2016, author's calculations.}}"

add_rows <-c(add_first_row,add_rows,add_source)

# Digits
tab.xtab <- xtable(tab[27:70,], 
                   caption=title, 
                   type="latex",
                   booktabs = FALSE,
                   auto=TRUE)
align(tab.xtab) <- c("X","X",rep(c("R{1.1cm}","R{1cm}"),4))

# Save table
setwd(pfolder)
mod_xtable(fn,tab.xtab,
           booktabs = TRUE, 
           caption.placement = "top",
           include.rownames = FALSE, include.colnames=FALSE,
           format.args=list(big.mark = ",",  decimal.mark = "."),
           tabular.environment = 'tabularx', width="1\\textwidth",
           sanitize.text.function = function(x){x},
           size = "footnotesize",
           #sanitize.colnames.function = identity, 
           sanitize.rownames.function = identity,
           hline.after=c(-1,39,41),
           add.to.row = list(pos = list(-1,0,44),
                             command = add_rows)
)
setwd(wd)