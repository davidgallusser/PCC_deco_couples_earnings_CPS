##############################################################
## A3 Prepare CPS data
##############################################################

setwd(paste0(wd,datafiles_d))

##############################################################
### Run loop to construct couple data set for every year

obs_per_year <- NULL

# UHRSWORKLY and WKSWORK1 available since 1976
years <- 1976:2020
for(i in years){
  print(i)
  load(file=paste0("CPS",i,".rda"))

##############################################################
## Sample selection 

# Total observations
obs_tot <- nrow(CPS)
  
#Select only individuals between 21 and 55
CPS <- CPS[which(CPS$AGE>=21&CPS$AGE<=55),]  
obs_age <- nrow(CPS)

#Select only adult civilians (1); armed forces (2) are excluded
CPS <- CPS[which(CPS$POPSTAT==1),]  
obs_civil <- nrow(CPS)

#Include  "Works for wages or salary" (20-28)
#and individuals who did not work (00, "not in universe"))
#Exclude self-employed and unpaid family worker
CPS <- CPS[which(is.element(CPS$CLASSWLY,c(0,20:28))),] 
obs_dependent <- nrow(CPS)


##############################################################
## Create own ID Variable to match spouses
CPS$PID <- CPS$YEAR*100000*100*100+CPS$MONTH*100000*100+CPS$SERIAL*100+CPS$PERNUM
CPS$HID <- CPS$YEAR*100000*100*100+CPS$MONTH*100000*100+CPS$SERIAL*100
CPS$SPOUSEID <- ifelse(CPS$SPLOC==0,NA,CPS$YEAR*100000*100*100+CPS$MONTH*100000*100+CPS$SERIAL*100+CPS$SPLOC)

##############################################################
## Define hourly weekly wages

## Working time variables
# Actually, there is only one way to construct hours per year:
#1)  ANNUALHOURS = UHRSWORKLY*WKSWORK1, only possible since 1976
# The second avenue would be to construct hours using the actual working hours
# (AHRSWORKT) variable that exists longer
# However, reference period of AHRSWORKT (and absent) is not "last year"
# but the interview week in March
#2)  if(ABSENT==1){
#         ANNUALHOURS <- AHRSWORKT * MEAN(WKSWORK2 by WKSWORK2-interval)
#      }else{
#         ANNUALHOURS <- MEAN(AHRSOWRKT by FULLPART category) * MEAN(WKSWORK2 by WKSWORK2-interval) 
#      } 


# weekly hours define NA 
#CPS$UHRSWORKLY <- recode_na(CPS$UHRSWORKLY, n=999)
CPS[which(CPS$UHRSWORKLY==999),"UHRSWORKLY"] <- 0

# Average weekly hours
CPS$ANNUALHOURS <- CPS$UHRSWORKLY*CPS$WKSWORK1
CPS$HOURS <- CPS$ANNUALHOURS/52

# Deflated hourly income
CPI99_2019 <- 0.663
CPS$WAGE <- (CPS$INCWAGE*(CPS$CPI99/CPI99_2019))/CPS$ANNUALHOURS
CPS[which(is.infinite(CPS$WAGE)|is.na(CPS$WAGE)),"WAGE"] <- 0     

# 'Top code' extreme hourly wages above 2000$
CPS[which(CPS$WAGE>2000),"WAGE"] <- 2000

summary(CPS$HOURS)
summary(CPS$WAGE)

#Define labor force participation variable 
CPS$LFP <- as.numeric(CPS$WAGE>0)
summary(CPS$LFP)

##############################################################
## Variables for wage imputation model
CPS[,"INCUNERN"] <- recode_na(CPS[,"INCUNERN"],n=c(99998,99999))
CPS[,"INCTOT"] <- recode_na(CPS[,"INCTOT"],n=c(999999999,999999998))
CPS[,"INCBUS"] <- recode_na(CPS[,"INCBUS"],n=c(99999999,99999998))
CPS[,"INCFARM"] <- recode_na(CPS[,"INCFARM"],n=c(99999999,99999998 ))

## Create unearned income for years after 1967
sel <- which(CPS$YEAR>1967)
CPS[sel,"INCUNERN"] <- CPS[sel,"INCTOT"] - CPS[sel,"INCWAGE_NI"]  - CPS[sel,"INCBUS"]  - CPS[sel,"INCFARM"]

## Create income from dividends, interests, and rents 
# INCIDR      income interests, dividendes, rents (1968-1975)
# INCDRT      income from dividends, rents, and trusts (1976-1987)
# INCINT      income from interest (since 1976)
# INCDIVID    income from dividends (since 1988)
# INCRENT     income from rents (since 1988)

# Recode na
CPS[,"INCIDR"] <- recode_na(CPS[,"INCIDR"],n=99999)
CPS[,"INCDRT"] <- recode_na(CPS[,"INCDRT"],n=99999)
CPS[,"INCINT"] <- recode_na(CPS[,"INCINT"],n=9999999)
CPS[,"INCDIVID"] <- recode_na(CPS[,"INCDIVID"],n=9999999)
CPS[,"INCRENT"] <- recode_na(CPS[,"INCRENT"],n=9999999)

# Construct consistent capital income variable
sel <- which(is.element(CPS$YEAR,1976:1987))
if(length(sel)>0){
CPS[sel,"INCIDR"] <-  CPS[sel,"INCDRT"] + CPS[sel,"INCINT"]
}
sel <- which(CPS$YEAR>1987)
if(length(sel)>0){
CPS[sel,"INCIDR"] <-  CPS[sel,"INCRENT"] + CPS[sel,"INCDIVID"] + CPS[sel,"INCINT"]
}

# Recode race variable
CPS[which(CPS$RACE>200),"RACE"] <- 300 #Recode all individual not Black/White to be other
CPS$RACE <- as.factor(CPS$RACE)

# Create Hispanic dummy
CPS[,"HISPAN"] <- as.numeric(is.element(CPS$HISPAN,100:412)) 

# Potential Experience
CPS[which(CPS$EDUC==999),"EDUC"] <- NA
CPS$EXP <-  memisc::recode(CPS$EDUC, 
                    0<-0:40,    #No HS to first year HS
                    1<-50, 
                    2<-60,
                    3<-70:73,   #High school
                    4<-80:81,
                    5<-90:92,
                    7<-100:111, #Bachelor's degree 
                    9<-120:123, #Master's degree
                    10<-124,    #Prof School
                    11<-125)    #PhD
CPS$EXP <-  CPS$AGE-14-CPS$EXP
CPS$EXP <- ifelse(CPS$EXP<0,0,CPS$EXP)

# Potential Experience 2
breaks <- c(0,5,10,15,20,25,30,100)
CPS$EXP2 <- cut(CPS$EXP, breaks=breaks, right=FALSE)

## Define factors
CPS$EDUC <- memisc::recode(CPS$EDUC, 
                           1<-0:60,    #No HS
                           2<-70:73,   #High school
                           3<-80:81,   #Some college 1
                           4<-90:92,   #Some college 2
                           5<-100:111, #Bachelor's degree 
                           6<-120:123, #Master's degree
                           7<-124,     #Prof School
                           8<-125)     #PhD

# Create second education variable to avoid perfect seperation in selection models
CPS$EDUC2 <- CPS$EDUC
CPS[which(CPS$EDUC2==6|CPS$EDUC2==7|CPS$EDUC2==8),"EDUC2"] <- 5

# Create factors
CPS$EDUC <- relevel(as.factor(CPS$EDUC),ref="2") # EDUC: Set high school as reference category
CPS$EDUC2 <- relevel(as.factor(CPS$EDUC2),ref="2") # EDUC: Set high school as reference category



# Geographical variables
CPS[which(CPS$METRO==9&CPS$YEAR==1976),"METRO"] <- NA
CPS$METRO <- relevel(as.factor(CPS$METRO),ref="1")  #METRO type: 0 not identifiable, 1 not in metro, 2 central city, 3 outside central city 4 central city status missing
CPS$STATEFIP <- relevel(as.factor(CPS$STATE),ref="6") #6=California 
CPS$REGION <- relevel(as.factor(CPS$REGION),ref="12") #12=Middle Atlantic Division

# Number of children older than 5
CPS$NCHILD518 <- CPS$NCHILD - CPS$NCHLT5 

##############################################################
## Merge couples and retrieve different sex couples only
cdata <- lapply(split(CPS,CPS$YEAR), function(x) merge_couples(pdata=x,
                                                               pid="PID",
                                                               parid="SPOUSEID", 
                                                               gendervar="SEX", 
                                                               level_male=1)[[1]]
                )
cdata <- do.call("rbind",cdata)

# Add household level data
cdata$INCUNERN_hh <- cdata$INCUNERN_m + cdata$INCUNERN_f  
cdata$INCIDR_hh <-  cdata$INCIDR_m + cdata$INCIDR_f 

##############################################################
# N Couples
obs_couples <- nrow(cdata)*2
obs_per_year <- rbind(obs_per_year,data.frame(YEAR=i,
                              TOTAL=obs_tot,
                              AGE=obs_age,
                              CIVILIAN=obs_civil,
                              DEPENDENT=obs_dependent,
                              COUPLES=obs_couples))

##############################################################
## Save data set
save(cdata, file=paste0("CPScouples",i,".rda"))

}

save(obs_per_year, file="observations_per_year.rda")
obs_per_year



