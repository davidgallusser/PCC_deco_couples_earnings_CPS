##############################################################
## A1 Import data sets
##############################################################

### Load IPUMS-CPS data set 1962-2020
# Data citation: 
# Sarah Flood, Miriam King, Renae Rodgers, Steven Ruggles and J. Robert Warren.
# Integrated Public Use Microdata Series, Current Population Survey: Version 7.0 [dataset]. 
# Minneapolis, MN: IPUMS, 2020.
# https://doi.org/10.18128/D030.V7.0

# Set path to working directory
setwd(paste0(wd,datafiles_d))

#########################################################
### Load all years (takes some time!)
memory.limit(size=16000)
CPSall <- read.table(gzfile("cps_00004.csv.gz"), sep=",", header=TRUE)
save(CPSall, file="CPSall.rda")

# Save years into single .rda files 
for(i in  unique(CPSall$YEAR)){
  print(i)
  CPS <- CPSall[which(CPSall$YEAR==i),]
  save(CPS, file=paste0("CPS",i,".rda"))
}

### Important variables

## Weighting variables
# ASECWT      Person level CPS-ASEC weight
# ASECWTH     HH level CPS-ASEC weight
# ASECFWT     Family CPS-ASEC weight

## ID
# CPSID        HH ID
# CPSIDP       Personal ID
# SPLOC        Spouse/Partner (Person Number!)
# PERNUM       Together with YEAR , MONTH, and SERIAL, PERNUM uniquely identifies each person within IPUMS-CPS samples, though not across IPUMS-CPS samples.
# YEAR         
# MONTH
# SERIAL

# PID          #constructed person id: paste(YEAR,MONTH,SERIAL,PERNUM, sep="-") 
# SID          #constructed spouse id: paste(YEAR,MONTH,SERIAL,SPLOC, sep="-") 

## Earnings 
# CPI99        Consumer Price Index (base 1999)

# INCWAGE      wage income = OINCWAGE + if(SRCEARN==1, INCLONGJ, 0)
# INCLONGJ     income longest job (top coded; needs to be imputed)
# OINCWAGE     other wage income (top coded; needs to be imputed))
# SRCEARN      Source of income for the job the individual holds longest (1=wage & salary, 2=self-employment, 3=farm self-employment, 4=working without pay)
# TINCLONGJ    top code flag  (0=not top coded, 1=top coded)
# TOINCWAGE    top code flag


# FULLPART    Worked full or part time last year (1=FULL TIME, 2=PART TIME, 9=UNKOWN, 0=NIU)
# UHRSWORKLY  usual weekly working hours last year (since 1976)
# WKSWORK1    weeks workd last year (since 1976)  
# WKSWORK2    weeks worked last year (intervalled, since 1962)

# AHRSWORKT   Hours worked last week (since 1992, = 0 if ABSENT)
# ABSENT      Absent from work last week 1=No,2=yes, laid off, 3=yes, other reason (vacation, illness, labor dispute) 
# UHRSWORKT 	Hours usually worked per week at all jobs (since 1994)
# UHRSWORK1   Hours usually worked per week at main job (since 1982)

# Zusätzliche Variablen:
# WKSWORK2    weeks worked last year (intervalled, since 1962)
# INCIDR      income interests, dividendes, rents (1968-1975)
# INCDRT      income from dividends, rents, and trusts (1976-1987)
# INCINT      income from interest (since 1976)
# INCDIVID    income from dividends (since 1988)
# INCRENT     income from rents (since 1988)
# WKSUNEM1    weeks unemployed last year (since 1976)  
# WKSUNEM2    weeks unemployed last year (intervalled, since 1962)
# FULLPART    fulltime-partime 