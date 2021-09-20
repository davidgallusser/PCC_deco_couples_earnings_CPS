##############################################################
## D8 - CPS - deco summary table 1976-2016
##############################################################

years <- c(1976,2016)
setwd(paste0(wd,datafiles_d))
#load(paste0("copula_decomposition_",years[1],"-",years[2],".rda"))
load(paste0("copula_decomposition_main_",years[1],"-",years[2],".rda"))
load(paste0("copula_decomposition_ref0_",years[1],"-",years[2],".rda"))
load(paste0("copula_decomposition_alt_order_",years[1],"-",years[2],".rda"))
load(paste0("copula_decomposition_parametric_",years[1],"-",years[2],".rda"))
load(paste0("copula_decomposition_cond_on_edu_",years[1],"-",years[2],".rda"))
setwd(wd)


##############################################################
# Decomposition table
create_stats_table(VCopObj, type="aggregate",
                    stats="inequality",
                    caption=NULL,
                    varnames=c("Wages men","Hours men","Wages women","Hours women"))

# Save tables
caption <- "Decomposition of changes in couples' earnings inequality, 1975-2015\\label{tab:deco_main}"
col.width <- 1.1

tab_xtab <- create_stats_table(VCopObj, type="aggregate",
                   stats="inequality",
                   caption=caption,
                   varnames=c("Wages men","Hours men","Wages women","Hours women"),
                   tex=TRUE, 
                   col.width=col.width )

cn <- colnames(tab_xtab)[seq(1,9,2)]
add_row <- paste("\\toprule\n",paste0(paste0("&\\multicolumn{2}{c}{",cn,"}"),collapse=""),"\\\\\n", collapse="")
add_source <- "\\bottomrule \n 
               \\multicolumn{11}{l}{\\footnotesize{Decomposition based on pair-copula construction discussed in \\ref{subsec:mainresults}. Standard errors in parenthesis based on 200 bootstrap replications.}}"

rownames(tab_xtab)[c(14,15)] <- c("\\quad 5 Wages men,Hours women;Wages w.",          
                                  "\\quad 6 Hours m.,Hours w.;Wages m.,Wages w.")

fn <- paste0(pfolder,"deco_1976_2016_main_results.tex")
setwd(wd)
mod_xtable(fn,tab_xtab,
           booktabs = TRUE,
           floating = TRUE, floating.environment = "sidewaystable",
           tabular.environment = 'tabularx', width="0.9\\textwidth",
           caption.placement = "top",
           size="small",
           include.rownames = TRUE, include.colnames=FALSE,
           sanitize.colnames.function = function(x){x}, 
           sanitize.rownames.function = function(x){x},
           format.args=list(big.mark = ",",  decimal.mark = "."),
           hline.after=c(0),
           add.to.row = list(pos = list(0,18),
                             command = c(add_row,add_source)))
setwd(wd)

