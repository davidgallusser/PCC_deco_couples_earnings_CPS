##############################################################
##
## Accounting for Inequality Trends in 
## Couples' Earnings with Pair-Copula Constructions:
## R code for the empirical analysis of CPS ASEC data
##
## David Gallusser, david.gallusser@unibas.ch
## September 2021
##
##############################################################
##
## To replicate results, set 'wd' to the directory 
## where the replication files are saved. 
## 
## Warning: 
## The bootstrap is time-consuming (~ 21 hours for 
## every decomposition). VC_deco() allows to switch
## off bootstrapping.
##
##############################################################
### Set working directories

# Uni
wd <- "C:/Users/gallusse/switchdrive/Uni Basel/Dissprojekte/6 Couples earnings decomposition/R code and data/PCC_deco_couples_earnings_CPS/"
setwd(wd)
# Laptop
wd<- "C:/Users/David Gallusser/switchdrive/Uni Basel/Dissprojekte/6 Couples earnings decomposition/R code and data/PCC_deco_couples_earnings_CPS/"
setwd(wd)
wd <- getwd()

# Data directory: CPS data is stored here; results will be save here
datafiles_d <- "/../PCC_deco_couples_earnings_CPS_data/" 

# Directory for plots and tables: Output tables and plots will be save here
pfolder <- "../../../_Thesis/Manuscript/Chapter3/plots/" #"../../plots/PCC_deco_couples_earnings_CPS/"

##############################################################
### Install required packages, load relevant libraries

# Packages to be loaded
package_Depends <- c("rvinecopulib",    # Fast estimation of vine copulas
                     "sampleSelection", # Sample selection/Heckman models
                     "Hmisc",           # weighted statistics
                     "kdecopula",       # Copula kernel density estimates
                     "parallel",        # Parallel computing
                     "pbapply",         # Progress bar for apply function
                     "dplyr",           # data wrangling
                     "reshape2",        # long-wide data transformation
                     "xtable",          # Export tables
                     "tictoc",          # time measurement
                     "ggplot2",         #plotting
                     "ggthemes",        #plotting 
                     "grid",            #plotting
                     "gridExtra",       #plotting
                     "egg"              #plotting
                     )

# Packages to be installed
package_Imports <-   c("memisc",
                       "tidyr")
packages <- c(package_Imports,package_Depends)
  
# Check if installed. If not, install binary version.
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  packages_to_install <- setdiff(packages, rownames(installed.packages()))
  cat(paste0("Packages to install:\n",paste0(packages_to_install,collapse=", ")))
  install.packages(packages_to_install ,
                   type="binary")
}

# Load packages
for(i in package_Depends){
  require(i,character.only = TRUE)
}

###############################################################
### Set parameters for output plots and tables
source("functions/0 - set graphical parameters for output plots.R")

###############################################################
### Load functions
setwd(wd)

### VC_deco() function to perform pcc decomposition
# 1: wrapper function to perform pcc decompostion; calls rvinecopulib::vinecop() to fit pcc 
# 2: functions to simulate counterfactuals and compute decomposition terms
# 3: functions to bootstrap decomposition procedure
# 4: summary function for VCdeco results (plots, tables)
# 5: DFL deco function for decomposition that controls for observables
source("functions/1 - VCdeco - estimation using rvinecopulib.R")
source("functions/2 - VCdeco - decomposition using rvinecopulib.R")
source("functions/3 - VCdeco - bootstrap.R")
source("functions/4 - VCdeco - summary functions for VC_deco.R")
source("functions/5 - VCdeco - DFL-deco-2-2.R")

### Additional statistic functions
# 6: Functions to analyze marginal distribution and dependencies
# 7: User written stats function (e.g., normalized weighted rank, Gini coefficient)
source("functions/6 - VCdeco - analyze marginals dependencies.R")
source("functions/7 - VCdeco - misc stats functions.R")

### Data preparation functions 
# 8: Functions to impute unobserved wages with a Heckman sample selection model
# 9: Functions to prepare data set with couple observations
# 10: Functions to impute topcoded CPS incomes
source("functions/8 - VCdeco - wage imputation.R")
source("functions/9 - VCdeco - data prep functions.R")
source("functions/10 - VCdeco - topcoded income imputation functions.R")

###############################################################
### A Prepare data
# 1 Load data from zip file and create annual data sets
source("A1 - CPS - import data.R")

# 2 Impute top coded wages 
setwd(wd)
source("A2 - CPS - impute top coded wages.R")

# 3 Construct couple data set 
setwd(wd)
source("A3 - CPS - prepare data.R")

# 4 Impute non-observed wages
setwd(wd)
source("A4 - CPS - wage imputation.R")

###############################################################
### Descriptive statistics
# B1 Compute descriptive statistics for all waves
setwd(wd)
source("B1 - CPS - descriptive statistics.R")

# B2 Estimate single bivariate copula densities
setwd(wd)
source("B2 - CPS - single bivariate copula densities.R")

###############################################################
# C Pair-copula decomposition

# C1 Decomposition for 1976-2016 
# - Baseline specification
# - Controlling for observables
# - Alternative specification for robustness checks
setwd(wd)
source("C1 - CPS - decomposition 1976-2016.R")

# C2 Decomposition for 1976-1996
setwd(wd)
source("C2 - CPS - decomposition 1976-1996.R")

# C3 Decomposition for 1996-2016
setwd(wd)
source("C3 - CPS - decomposition 1996-2016.R")

# C4 Check GoF for different PCC for 1976-2016
setwd(wd)
source("C4 - CPS - decomposition 1976-2016 - GoF for different PCC.R")

###############################################################
# D Plots and tables

# D1 Plot descriptive statistics
# (Figures 1, C.1, C.2, and C.3)
setwd(wd)
source("D1 - CPS - descriptive statistics - plots.R")

# D2 Descriptive statistics table
# (Table C.1)
setwd(wd)
source("D2 - CPS - descriptive statistics - table.R")

# D3 Decomposition results: "Quantile treatment effect" plots, baseline decomposition 1976-2016
# (Figures 3 and E.1)
setwd(wd)
source("D3 - CPS - deco quantile plots 1976-2016 main results.R")

# D4 Decomposition results: D9/D1 summary plots 1976-2016
# (Figures 5, 6, 7, 8, and E.2)
setwd(wd)
source("D4 - CPS - deco summary plots 1976-2016.R")

# D5 Decomposition results: D9/D1 summary plots ubperiods 1976-1996-2016
# (Figure 4)
setwd(wd)
source("D5 - CPS - deco summary plots subperiods 1976-1996-2016.R")

# D6 vine graph and copula density plots 1976-2016
# (Figures 2, D.1, D.3, and D.4)
setwd(wd)
source("D6 - CPS - copula density plots 1976-2016.R")

# D7 GoF statistics and plots
# (Tables D.1, D.2, and D.3, Figure D.2)
setwd(wd)
source("D7 - CPS - deco GoF statistics and plots 1976-2016.R")

# D8 Decomposition results: Summary table 1976-2016
# (Table E.1)
source("D8 - CPS - deco summary table.R")

