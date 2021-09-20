##############################################################
## D7 - CPS - deco GoF statistics and plots 1976 - 2016
##############################################################

setwd(paste0(wd,datafiles_d))
years <- c(1976,1996)
load(paste0("copula_decomposition_main_",years[1],"-",years[2],".rda"))
VCopObj7696 <- VCopObj
years <- c(1976,2016)
load(paste0("copula_decomposition_main_",years[1],"-",years[2],".rda"))
load(paste0("copula_best_models_main_",years[1],"-",years[2],".rda"))
load(paste0("copula_decomposition_parametric_",years[1],"-",years[2],".rda"))
setwd(wd)

##############################################################
# Q-Q plots
logscale=FALSE
fit_qqplot <- VC_deco_qqplot(VCopObj,formula=NULL, type=1, logscale=logscale, group=c("1975","2015"), tau=seq(1,99,2)/100) + 
  labs(title=NULL,y="Simulated earnings in 2019 US-$", x="Observed earnings in 2019 US-$")

# Save
fn <- paste0(pfolder,"qqplot_fit",ptype)
ggsave(fn, plot=fit_qqplot, units=punits, scale=pscale, width=pwidth, height=pheight)
  
##############################################################
# Kolmogorov-Smirnov and Cramér-von-Mises-Statistics

cat("p-values: observed dist. = fitted dist.")
stats_1 <- VCopObj$VCdeco_se$aggregate$quants_sim_vs_obs_KS_CMS_test
stats_2 <- VCopObj7696$VCdeco_se$aggregate$quants_sim_vs_obs_KS_CMS_test

tab <- data.frame(test=c("Kolmogorov-Smirnov","Cram\\'er-von-Mises-Smirnov"),
           `1976`=stats_1[3:4],
           `1996`=stats_2[1:2],
           `2016`=stats_1[1:2])
rownames(tab) <- tab$test
names(tab) <- c("test",1976,1996,2016)
tab <- tab[,-1]
title <- "p-values of bootstrapped goodness of fit tests \\\\($H_0$: observed earnings dist.=simulated earnings dist.) \\label{tab:ks}"
fn <- "KS-CMS-tests.tex"
add_rows <-  "\\bottomrule \n 
              \\multicolumn{4}{l}{\\multirow{3}{*}{
                                 \\footnotesize{\\shortstack[l]{Test statistics are evaluated at every ventile and at the 99th percentile, respectively. We \\\\ bootstrap the distributions of the test statistics as implemented by \\cite{Chen.2017}. \\\\ The bootstrap is based on 200 replications.}}
                                                  }} \\\\ \n"

# Digits
tab.xtab <- xtable(tab, 
                   caption=title, 
                   type="latex",
                   booktabs = TRUE,
                   auto=TRUE)
digits(tab.xtab) <- matrix(rep(3,4*2), ncol=4, byrow=FALSE)
align(tab.xtab) <- c("X",rep("R{2cm}",3))

# Save table
setwd(pfolder)
mod_xtable(fn,tab.xtab,
           booktabs = TRUE, 
           caption.placement = "top",
           include.rownames = TRUE, include.colnames=TRUE,
           format.args=list(big.mark = ",",  decimal.mark = "."),
           tabular.environment = 'tabularx', width="0.8\\textwidth",
           sanitize.text.function = function(x){x},
           size = "small",
           #sanitize.colnames.function = identity, 
           #sanitize.rownames.function = identity,
           hline.after=c(-1,0),
           add.to.row = list(pos = list(2),
                             command = add_rows)
)
setwd(wd)

##############################################################
# AIC/BIC different models
best_models <- VC_best_models$vine_IC %>% mutate(rank0 = rank(AIC0),
                                                 rank1 = rank(AIC1)) 

best_models <- best_models[,c(1,6,2,3,7,4,5)]
tab <- best_models  %>% arrange(rank1)

rownames(tab) <- tab[,1]
rownames(tab)[8] <- paste0("\\textbf{",rownames(tab)[8],"}", collapse="")
tab <- tab[,-1]
cn <- c("",rep(c("AIC rank","AIC","BIC"),2))
title <- "Information criteria of best fitting pair-copula constructions \\label{tab:ic}"
fn <- "ic-table.tex"
add_first_row <- "\\toprule\n \\textbf{Vine structure} & \\multicolumn{3}{c}{\\textbf{1975}}   &
                          \\multicolumn{3}{c}{\\textbf{2015}}   \\\\ \\cmidrule(r){2-4}\\cmidrule{5-7} \n "
add_colnames <- paste0(add_first_row,paste0(cn, collapse="&"),"\\\\\n",collapse="")
add_caption <- "\\bottomrule \n 
                \\multicolumn{7}{l}{\\multirow{3}{*}{
                                 \\footnotesize{\\shortstack[l]{Vine structure defines PCC by edges in first tree. Hyphns indicate edges. We enumerate variables as \\\\ follows: 1=Wages men, 2=Hours men, 3=Wages women, and 4=Hours women. For instance, D-vine \\\\ 2-1-3-4 is the PCC  that includes pair-copulas between 2\\&1, 1\\&3,  and 3\\&4 in first tree.}}
                                                  }} \\\\ \n"
  
  
# Digits
tab.xtab <- xtable(tab, 
                   caption=title, 
                   type="latex",
                   booktabs = TRUE,
                   auto=TRUE)
digits(tab.xtab) <- matrix(rep(0,7*16), ncol=7, byrow=FALSE)
align(tab.xtab) <- c("X",rep("R{1.5cm}",6))


setwd(pfolder)
mod_xtable(fn,tab.xtab,
           booktabs = TRUE,
           tabular.environment = 'tabularx', width="0.95\\textwidth",
           caption.placement = "top",
           size="small",
           include.rownames = TRUE, include.colnames=FALSE,
           sanitize.text.function = function(x){x},
           format.args=list(big.mark = ",",  decimal.mark = "."),
           hline.after=c(-1),
           add.to.row = list(pos = list(-1,16),
                             command = c(add_colnames,add_caption))
            )
setwd(wd)


##############################################################
# Table with  paramertic copula estimation
res0 <- summary(VCopObj_parest$VCfit$VCfit0)
res1 <- summary(VCopObj_parest$VCfit$VCfit1)

res <- rbind(res0,res1) 
res$conditioned <- do.call("c", lapply(res$conditioned, function(x) paste0(x, collapse=",")))
res$conditioning <- do.call("c", lapply(res$conditioning, function(x) paste0(x, collapse=",")))
res$parameters <- do.call("c", lapply(res$parameters , function(x) paste0(round(x,2), collapse=", ")))
res$tau <- round(res$tau,2)
res$loglik <- round(res$loglik,0)
res <- res %>% as.data.frame

res$family <- recode(res$family,
       "bb1"="Clayton-Gumbel (BB1)",
       "t"="Student t",
       "joe"="Joe",
       "gumbel"="Gumbel",
      "clayton"="Clayton",
     "bb8"="Joe-Frank (BB8)")
res <- res[,-c(2,5,9)]
names(res) <- c("Tree","Var.", "Cond. Var.", "Model","Rotation", "Par. est.", "$\\tau$","loglik")

add_rows <- c("\\midrule \n 
               \\multicolumn{7}{l}{\\textbf{1975}}\\\\ \n",
              "\\multicolumn{7}{l}{\\textbf{2015}}\\\\ \n")
add_caption <- "\\bottomrule \n 
                \\multicolumn{7}{l}{\\multirow{1}{*}{
                                 \\footnotesize{\\shortstack[l]{Variables: 1=Wages men, 2=Hours men, 3=Wages women, 4=Hours women}}
                                                  }} \\\\ \n"
add_rows <- c(add_rows,add_caption)
title <- "Selected copula models and estimates of parametric estimation \\label{tab:parest}"
fn <- "parametric_estimates.tex"

# Digits
tab.xtab <- xtable(res, 
                   caption=title, 
                   type="latex",
                   booktabs = TRUE,
                   auto=TRUE)
#digits(tab.xtab) <- matrix(rep(0,7*16), ncol=7, byrow=FALSE)
align(tab.xtab) <- c(rep("L{0.7cm}",2),rep("L{0.7cm}",1),rep("L{1.7cm}",1),"X",rep("R{1.5cm}",1),rep("R{1.8cm}",1),rep("R{0.9cm}",2))


setwd(pfolder)
mod_xtable(fn,tab.xtab,
           booktabs = TRUE,
           tabular.environment = 'tabularx', width="0.95\\textwidth",
           caption.placement = "top",
           size="small",
           include.rownames = FALSE, include.colnames=TRUE,
           sanitize.text.function = function(x){x},
           format.args=list(big.mark = ",",  decimal.mark = "."),
           hline.after=c(-1),
           add.to.row = list(pos = list(0,6,12),
                             command = add_rows)
)
setwd(wd)
