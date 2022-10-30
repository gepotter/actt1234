# Gail Potter
# MAY 2021
# ACTT1/2/3 post hoc analysis- create analysis data sets
# R version 4.0.4
# inputs: adtte, adsl, adcm, adlb, advs from ACTT 1, 2, and 3
# outputs: analysis.data.Rda (only includes RDV arm), 
# allarms.Rda (same thing but all arms- use for checking against CSR) 

library(haven)
library(r2rtf)
library(dplyr)
library(sessioninfo)
library(reshape)
library(readxl)
library(pROC)


rm(list=ls())
indata1="N:/COVID-19/ACTT/ACTT-1/Data/"
indata2="N:/COVID-19/ACTT/ACTT-2/Data/"
indata3="N:/COVID-19/ACTT/ACTT-3/Data/2021_02_24_Lock/"
in.crp.updated ="N:/COVID-19/ACTT/ACTT-1/Data/24MAY2022 CRP"
mb_update = 'N:/COVID-19/ACTT/ACTT-1/Data/MB_update_31MAR2021'

outdir="N:/COVID-19/ACTT/Combined ACTT Analyses/ACTT123 comparison of RDV arm/output"



actt1=read_sas(paste(indata1,"adtte.sas7bdat",sep=""))
actt2=read_sas(paste(indata2,"adtte.sas7bdat",sep=""))
actt3=read_sas(paste(indata3,"adtte.sas7bdat",sep=""))
select=dplyr::select
ttedat = full_join(actt1, full_join(actt2, actt3)) %>%
  select(AVAL, CNSR, USUBJID, PARAMCD)

# subject level analysis data set 
rename=dplyr::rename

adsl1 = read_xpt("N:/COVID-19/ACTT/ACTT-1/Data/adsl.xpt")%>%
  mutate(cohort = recode(as.character(TRT01A),'Remdesivir' = "ACTT1"), phase=1)
adsl2 = read_sas("N:/COVID-19/ACTT/ACTT-2/Data/adsl.sas7bdat")%>% 
  mutate(cohort = recode(as.character(TRT01A),'Placebo + RDV' = "ACTT2"),phase=2) %>%
           rename(HYPFL=HYPTFL, CANCERFL=CANCFL, CORFL=CORQFL)
adsl3 = read_sas("N:/COVID-19/ACTT/ACTT-3/Data/2021_02_24_Lock/adsl.sas7bdat")%>% 
  mutate(cohort = recode(as.character(TRT01A),'Placebo + RDV' = "ACTT3"),phase=3) %>%
  rename(HYPFL=HYPTFL, CANCERFL=CANCFL, CORFL=CORQFL)

adsl3$CLDFL='N'  ## chronic liver disease was an exclusion criterion.
select=dplyr::select
adsl=full_join(adsl1, full_join (adsl2,adsl3))%>%
  mutate(BMIFL=BMI>=35&!is.na(BMI), BMI30FL = BMI>=30 &!is.na(BMI), 
         DIABFL=as.numeric((DIAB1FL=='Y' & !is.na(DIAB1FL))| 
         (DIAB2FL=='Y')&!is.na(DIAB2FL))) %>%
select(AGE, SEX, RACE, ETHNIC, BMIFL, USUBJID, ITTFL,  TRTFL, DIABFL,phase,
       DIAB1FL, DIAB2FL,BMI30FL, CKDFL, CLDFL, CRDFL, ASTHMAFL,
       BDURSYMP, CORFL,  CADFL,HYPFL,  CHFFL,  CANCERFL, IMMDFL,BCRP, DSCHDT,
       cohort, BCSOSN, SITEID, SITENAME, TRTFL, BMI, RANDDT) %>%
  mutate(CKDFL=!is.na(CKDFL)&as.numeric(CKDFL=='Y'), CLDFL=!is.na(CLDFL)&as.numeric(CLDFL=='Y'), 
         CRDFL=!is.na(CRDFL)&as.numeric(CRDFL=='Y'), ASTHMAFL=!is.na(ASTHMAFL)&as.numeric(ASTHMAFL=='Y'), 
         CORFL=!is.na(CORFL)&as.numeric(CORFL=='Y'),  
         CADFL=!is.na(CADFL)&as.numeric(CADFL=='Y'),HYPFL=!is.na(HYPFL)&as.numeric(HYPFL=='Y'),  
         CHFFL=!is.na(CHFFL)&as.numeric(CHFFL=='Y'),  CANCERFL=!is.na(CANCERFL)&as.numeric(CANCERFL=='Y'), 
         IMMDFL=!is.na(IMMDFL)&as.numeric(IMMDFL=='Y')) %>%
mutate(num.comorbidities = 
   as.numeric(BMI30FL) + (DIABFL) + CKDFL + CLDFL + CRDFL + ASTHMAFL +CORFL +
   CADFL + HYPFL + CHFFL+ CANCERFL+ IMMDFL )

# ADEF data set- for NEWS
adef1=read_xpt(paste(indata1, 'adef.xpt', sep='')) %>% 
  filter(PARAMCD=="NEWSTOT" & AVISITN==0)
adef2=read_sas(paste(indata2, 'adef.sas7bdat', sep=''))%>% 
  filter(PARAMCD=="NEWSTOT" & AVISITN==0)
adef3=read_sas(paste(indata3, 'adef.sas7bdat', sep=''))%>% 
  filter(PARAMCD=="NEWSTOT" & AVISITN==0)
adef=full_join(full_join(adef1,adef2),adef3) %>% mutate(news = AVAL) %>%
  select(USUBJID, news)

# ADEF data set- for NEWS
c1=read_xpt(paste(indata1, 'adef.xpt', sep='')) %>% 
  filter(PARAMCD=="NEWSCONS" & AVISITN==0)
c2=read_sas(paste(indata2, 'adef.sas7bdat', sep=''))%>% 
  filter(PARAMCD=="NEWSCONS" & AVISITN==0)
c3=read_sas(paste(indata3, 'adef.sas7bdat', sep=''))%>% 
  filter(PARAMCD=="NEWSCONS" & AVISITN==0)
consc=full_join(full_join(c1,c2),c3) %>% mutate(consciousness = AVAL) %>%
  select(USUBJID, consciousness)


# ADEF data set- for hospital stay duration
h1=read_xpt(paste(indata1, 'adef.xpt', sep='')) %>% 
  filter(PARAMCD=="DAYSHOSP")
h2=read_sas(paste(indata2, 'adef.sas7bdat', sep=''))%>% 
  filter(PARAMCD=="DAYSHOSP")
h3=read_sas(paste(indata3, 'adef.sas7bdat', sep=''))%>% 
  filter(PARAMCD=="DAYSHOSP") 
durations=full_join(full_join(h1,h2),h3) %>% mutate(dayshosp = AVAL) %>%
  select(USUBJID, dayshosp)

# HO data set- for time from symptom onset to hospitalization
ho1=read_xpt(paste(indata1, 'HO.xpt', sep='')) %>% 
  filter(HOCAT=="INITIAL HOSPITALIZATION") 
ho2=read_xpt(paste(indata2, 'HO.xpt', sep=''))%>% 
  filter(HOCAT=="INITIAL HOSPITALIZATION")
ho3=read_sas(paste(indata3, 'HO.sas7bdat', sep=''))%>% 
  filter(HOCAT=="INITIAL HOSPITALIZATION")
ho=full_join(full_join(ho1,ho2),ho3) %>% 
  mutate(hospdate = as.Date(HOSTDTC)) %>%
  select(USUBJID, hospdate) 




## grab baseline vital signs:

advs1=read_xpt(paste(indata1, "VS.xpt", sep='')) 
advs2=read_xpt(paste(indata2, "VS.xpt", sep='')) 
advs3=read_xpt(paste(indata3, "VS.xpt", sep='')) 


vs1 = full_join(full_join(advs1, advs2), advs3) %>%
  filter(VISITDY==1)


## One subject had duplicate results; we analyzed the second:
table(table(vs1$USUBJID, vs1$VSTEST))
vs1$idtest=paste(vs1$USUBJID, vs1$VSTEST)
# View(vs1[which(duplicated(vs1$idtest)),])
# View(vs1[which(vs1$USUBJID=="COV.03114"),])
vs = vs1 %>% arrange(idtest, desc(VSDTC))
vs=vs[!duplicated(vs$idtest) & !is.na(vs$VSSTRESN),] %>% 
  mutate(test=paste(VSTEST, " (", VSSTRESU,")",sep=''))

subj.vitals = cast(vs, USUBJID~VSTEST, value='VSSTRESN')


## grab baseline lab values:
lb1=read_xpt(paste(indata1,"adlb.xpt",sep=''))%>% filter(AVISIT=="Baseline")
lb2=read_sas(paste(indata2,"adlb.sas7bdat",sep=''))%>% filter(AVISIT=="Baseline")
lb3=read_sas(paste(indata3,"adlb.sas7bdat",sep=''))%>% filter(AVISIT=="Baseline")

labs0 = full_join(full_join(lb1, lb2), lb3) %>%
  filter( !PARAM %in% c('ABO Blood Group' , 
                        'Lab Data','Coagulation Data',
                        'Prothrombin Intl. Normalized Ratio',
                        'Prothrombin Time (sec)',
                        'Prothrombin Intl. Normalized Ratio (RATIO)',
                        'Choriogonadotropin Beta',
                        'Creatinine Clearance, Cockcroft-Gault Method, (mL/min)',
                        'Glucose (mg/dL)',
                        "C Reactive Protein (mg/L)",  ## baseline CRP will be obtained from adsl
                        "D-Dimer (mg/L FEU)"))


## Grab CRP values for ACTT1
select=dplyr::select

## ADDED 24 MAY 2022
crp = read_sas(paste(in.crp.updated, '/actt1_crp_retest.sas7bdat',sep='')) %>%
  rename(USUBJID=usubjid, BCRP=CRP) %>% select(USUBJID, BCRP)

crp.actt1=left_join(crp, adsl%>% select(USUBJID,TRTFL,  cohort, phase))
crp.all = full_join(crp.actt1, adsl%>%filter(phase %in% 2:3) %>%
          select(BCRP, USUBJID, cohort, TRTFL)) %>% filter(TRTFL=='Y')


original.labs = full_join(labs0, adsl %>% select(USUBJID,TRTFL,  cohort)) %>% 
  filter(!is.na(PARAM) & !is.na(AVAL))  

## Recode outliers
labs1=original.labs %>% 
  mutate(AVAL=ifelse(PARAMCD=='CREAT'  & AVAL>5, NA, AVAL),
         AVAL=ifelse(PARAMCD=='HGB' & AVAL>30, NA, AVAL),
         AVAL=ifelse(PARAMCD %in% c('NEUT','WBC', 'LYM', 'MONO') & AVAL>100, AVAL/1000, AVAL),
         AVAL=ifelse(PARAMCD %in% c('NEUT', 'WBC', 'LYM', 'MONO') & AVAL<.01, AVAL*1000,AVAL),
         AVAL=ifelse(PARAMCD =='PLAT' & AVAL<1, AVAL*1000,AVAL),
         AVAL=ifelse(PARAMCD =='EOS' & AVAL>1, AVAL/1000,AVAL),
         AVAL=ifelse(PARAMCD =='BASO' & AVAL>6, AVAL/1000,AVAL),
         PARAM=ifelse(PARAM=='Glomerular Filtration Rate, Estimated (mL/min/1.73m2)',
                      'eGFR (mL/min/1.73m2)', PARAM))

subj.labs = cast(labs1, USUBJID~PARAMCD, value='AVAL') %>%
  mutate(BASOLE = BASO/WBC,
         EOSLE = EOS/WBC,
         LYMLE = LYM/WBC,
         MONOLE = MONO/WBC,
         NEUTLE = NEUT/WBC)

## Recode ratios > 1 to missing as well as their components
neutle.ind = which(subj.labs$NEUTLE>1)
subj.labs$NEUTLE[neutle.ind]=NA
subj.labs$NEUT[neutle.ind]=NA
subj.labs$WBC[neutle.ind]=NA

which(subj.labs$EOSLE>1) ## NONE

lym.ind = which(subj.labs$LYMLE>1)
subj.labs$LYMLE[lym.ind]=NA
subj.labs$LYM[lym.ind]=NA
subj.labs$WBC[lym.ind]=NA

mon.ind = which(subj.labs$MONOLE>1)
subj.labs$MONOLE[mon.ind]=NA
subj.labs$MONO[mon.ind]=NA
subj.labs$WBC[mon.ind]=NA

which(subj.labs$BASOLE>1)  ## NONE

newratios.long = melt(subj.labs, id.vars='USUBJID')
names(newratios.long)[2]='AVAL'


labs2 = full_join(labs1 %>% select(!'AVAL'),
                  newratios.long %>% filter(!is.na(AVAL) )) 

## need to grab parameter names
param.mapping = labs2 %>% select(PARAM, PARAMCD) %>%
  filter(!duplicated(PARAM) & !is.na(PARAM))

labs3=left_join(labs2 %>% select(,!PARAM),param.mapping, by='PARAMCD') 
dim(labs3)

## Need to grab cohort
labs=left_join(labs3 %>% select(!cohort), adsl %>% select(USUBJID, cohort))



## For all parameters with values 10^9/L, change units to 10^3/mcL (same value)
flag = substr(labs$PARAM, nchar(labs$PARAM)-7, nchar(labs$PARAM))=='(10^9/L)'
param.name = substr(labs$PARAM, 1,nchar(labs$PARAM)-9)
labs$PARAM[which(flag)] = paste(param.name[which(flag)], '(10^3/mcL)')

tapply(labs$AVAL, labs$PARAM, summary)

## log scale
labs$AVAL_nozero = labs$AVAL
labs$AVAL_nozero[labs$AVAL==0]=NA

mins=labs %>% group_by(PARAM) %>% 
  summarise(min=min(AVAL_nozero, na.rm=TRUE))
labs.transformed=left_join(labs,mins)
labs.transformed$log.AVAL_nozero=log(labs.transformed$AVAL_nozero)


labs.transformed$log.AVAL_nozero[which(labs.transformed$AVAL==0)]=
  log(labs.transformed$min[which(labs.transformed$AVAL==0)]/2)


subj.labs2 = cast(labs.transformed, USUBJID~PARAMCD, value='log.AVAL_nozero') 
names(subj.labs2)[2:length(names(subj.labs2))]=
  paste("LOG.", names(subj.labs2)[2:length(names(subj.labs2))],sep='')


cov1 = full_join(adef,adsl %>% select(-'BCRP'))
cov2 = full_join(cov1, ho)
cov3 = full_join(cov2, subj.vitals)
cov4 = full_join(cov3, subj.labs)
cov5 = full_join(cov4, subj.labs2)
cov6 = full_join(cov5, consc)
cov7 = full_join(cov6, durations)
cov8 = full_join(cov7, crp.all %>% select(USUBJID, BCRP))


cov8$'Respiratory Rate'[which(cov8$BCSOSN==7)]=30
cov8$'Oxygen Saturation'[which(cov8$BCSOSN==7)]=90
cov8$ETHNIC[cov8$ETHNIC =='NOT REPORTED'] = 'NOT REPORTED/UNKNOWN'
cov8$ETHNIC[cov8$ETHNIC =='UNKNOWN'] = 'NOT REPORTED/UNKNOWN'

covariates = cov8 %>% select(-c(Weight, Height, BMI)) %>%
  mutate(symp.onset.to.hospitalization = BDURSYMP + hospdate-RANDDT) %>%
  mutate(bcsos=paste(BCSOSN, ": ", 
                     recode(BCSOSN,
                            '4'="hospitalized, not requiring supplemental oxygen but requiring ongoing medical care (related to Covid-19 or to other medical conditions)",
                            '5'="hospitalized, requiring any supplemental oxygen",
                            '6'="hospitalized, requiring noninvasive ventilation or use of high-flow oxygen devices",
                            '7'="hospitalized, receiving invasive mechanical ventilation or extracorporeal membrane oxygenation (ECMO)"),sep=''),
         hisp=factor(ETHNIC=='HISPANIC OR LATINO'), HYPFL=factor(HYPFL), IMMDFL=factor(IMMDFL),
         RACE=factor(RACE),
         logcrp=log(BCRP),
         Black=as.numeric(RACE=='BLACK OR AFRICAN AMERICAN'), 
         Asian=as.numeric(RACE=='ASIAN'), symp.onset.to.hospitalization=as.numeric(symp.onset.to.hospitalization),
         male=as.numeric(SEX=='M')) %>%
  mutate(mort4c = 2*(AGE>=50)+2*(AGE>=60)+2*(AGE>=70)+1*(AGE>=80) + male + (num.comorbidities==1) +
           (num.comorbidities>1) + (`Respiratory Rate` >=20 | BCSOSN==7) + (`Respiratory Rate` >=30| BCSOSN==7) +
           2*(`Oxygen Saturation`<92 | BCSOSN>4) + 
           (consciousness==3) + (BCRP>=50) + (BCRP>=100),
         hosp.duration = dayshosp + as.numeric(covariates$RANDDT - covariates$hospdate))
    ## hosp.duration = as.numeric(DSCHDT-hospdate))  ## updated 20 OCT 2021
covariates$symp.onset.to.hospitalization[covariates$symp.onset.to.hospitalization<0 & 
                                           !is.na(covariates$symp.onset.to.hospitalization)]=0

all=left_join(ttedat, covariates) 

# Censor people at Day 28:  ## UPDATED from 29 to 28 on 19JUL2021
all$CNSR[all$AVAL>28]=1
all$AVAL[all$AVAL>28]=28


analysis.data = all %>%
filter(cohort %in% c("ACTT1",'ACTT2','ACTT3') & TRTFL=='Y') 
setwd(outdir)
save(analysis.data, file="analysis.data.Rda")
save(all, file="allarms.Rda")
save(mins, file='actt123_mins')
save.image("actt123data.rdata")
