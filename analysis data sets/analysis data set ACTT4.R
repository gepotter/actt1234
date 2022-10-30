# Gail Potter
# APR 2022
# Create analysis data set for ACTT 3 vs 4 comparison
# R version 4.0.4
# inputs: adtte, adsl, adcm, adlb, advs from ACTT 3 and 4
# outputs: 

library(arsenal)
library(haven)
library(r2rtf)
library(dplyr)
library(sessioninfo)
library(reshape)
library(readxl)
library(pROC)


rm(list=ls())
indata2="N:/COVID-19/ACTT/ACTT-2/Data/"
indata4.sdtm="N:/COVID-19/ACTT/ACTT-4/Data/ACTT-4 SDTM Data Transfer Package/"
indata4.adam="N:/COVID-19/ACTT/ACTT-4/Data/ACTT-4 ADaM Data Transfer Package/"

outdir="N:/COVID-19/ACTT/Combined ACTT Analyses/ACTT123 comparison of RDV arm/output"

actt2=read_sas(paste(indata2,"adtte.sas7bdat",sep=""))
actt4=read_xpt(paste(indata4.adam,"adtte.xpt",sep=""))
select=dplyr::select
ttedat = full_join(actt2, actt4) %>%select(AVAL, CNSR, USUBJID, PARAMCD)

setwd(outdir)
save.image('data24.rdata')

# subject level analysis data set 
rename=dplyr::rename

adsl2.0 = read_sas(paste(indata2,"adsl.sas7bdat",sep=''))
adsl2 = adsl2.0  %>% mutate(phase=2) %>%
  rename(HYPFL=HYPTFL, CANCERFL=CANCFL, CORFL=CORQFL)


select=dplyr::select

save.image('data24.rdata')

adsl4.0 = read_xpt(paste(indata4.adam,"adsl.xpt",sep=''))

adsl4= adsl4.0 %>% mutate(phase=4) %>%
  rename(HYPFL=HYPTFL, CANCERFL=CANCFL) %>% rename(CORFL=PCOVSOFL)
  

adsl=full_join (adsl2,adsl4)%>%
  mutate(BMIFL=BMI>=35&!is.na(BMI), BMI30FL = BMI>=30 &!is.na(BMI), 
               DIABFL=as.numeric((DIAB1FL=='Y' & !is.na(DIAB1FL))| 
                                   (DIAB2FL=='Y')&!is.na(DIAB2FL))) %>%
  select(AGE, SEX, RACE, ETHNIC, BMIFL, USUBJID, ITTFL,  TRTFL, DIABFL,phase,
         DIAB1FL, DIAB2FL,BMI30FL, CKDFL, CRDFL, ASTHMAFL,CORFL,
         DEXMTHFL, ## added 05 MAY 2022
         BDURSYMP, CADFL,HYPFL,  CHFFL,  CANCERFL, IMMDFL,BCRP, DSCHDT,
         phase, BCSOSN, SITEID, SITENAME, TRTFL, BMI, RANDDT) %>%
  mutate(CKDFL=!is.na(CKDFL)&as.numeric(CKDFL=='Y'), 
         CRDFL=!is.na(CRDFL)&as.numeric(CRDFL=='Y'), 
         ASTHMAFL=!is.na(ASTHMAFL)&as.numeric(ASTHMAFL=='Y'), 
         CADFL=!is.na(CADFL)&as.numeric(CADFL=='Y'),
         HYPFL=!is.na(HYPFL)&as.numeric(HYPFL=='Y'),  
         CHFFL=!is.na(CHFFL)&as.numeric(CHFFL=='Y'),  
         CANCERFL=!is.na(CANCERFL)&as.numeric(CANCERFL=='Y'), 
         IMMDFL=!is.na(IMMDFL)&as.numeric(IMMDFL=='Y'),
         CORFL = !is.na(CORFL) & as.numeric(CORFL=='Y')) %>%
mutate(num.comorbidities = 
   as.numeric(BMI30FL) + (DIABFL) + CKDFL + CRDFL + ASTHMAFL +
   CORFL + CADFL + HYPFL + CHFFL+ CANCERFL+ IMMDFL )

save.image('data24.rdata')


# ADEF data set- for NEWS
adef2.0=read_sas(paste(indata2, 'adef.sas7bdat', sep=''))
adef4.0=read_xpt(paste(indata4.adam, 'adef.xpt', sep=''))



adef2 = adef2.0 %>% filter(PARAMCD=="NEWSTOT" & AVISITN==0)
adef4= adef4.0 %>% filter(PARAMCD=="NEWSTOT" & AVISITN==0)
save.image('data24.rdata')


adef=full_join(adef2,adef4) %>% mutate(news = AVAL) %>%
  select(USUBJID, news)

# RS data set- for consciousness

rs2.0=read_xpt(paste(indata2, 'RS.xpt', sep=''))
rs2 = rs2.0 %>% 
  filter(RSTESTCD == 'NEWS0107' & VISITNUM==101) %>% 
  mutate(consciousness = ifelse(RSORRES=='Alert', 0, 3)) %>% 
  select(USUBJID, consciousness) %>%
  filter(!duplicated(USUBJID)) ## one person had duplicate, identical baseline measurements

rs4.0=read_xpt(paste(indata4.sdtm, 'RS.xpt', sep=''))

rs4 = rs4.0 %>% 
  filter(RSTESTCD == 'NEWS0107' & VISITNUM==101) %>% 
  mutate(consciousness = ifelse(RSORRES=='Alert', 0, 3)) %>% 
  select(USUBJID, consciousness)

consc=full_join(rs2,rs4) 


# ADEF data set- for hospital stay duration
h2=adef2.0%>%  filter(PARAMCD=="DAYSHOSP") 
h4= adef4.0 %>%   filter(PARAMCD=="DAYSHOSP") 
durations=full_join(h2,h4) %>% mutate(dayshosp = AVAL) %>%
  select(USUBJID, dayshosp)


save.image('data24.rdata')



# HO data set- for time from symptom onset to hospitalization
ho2=read_sas(paste(indata2, 'HO.sas7bdat', sep=''))%>% 
  filter(HOCAT=="INITIAL HOSPITALIZATION")
ho4=read_xpt(paste(indata4.sdtm, 'HO.xpt', sep=''))%>% 
  filter(HOCAT=="INITIAL HOSPITALIZATION")
ho=full_join(ho2, ho4) %>% 
  mutate(hospdate = as.Date(HOSTDTC)) %>%
  select(USUBJID, hospdate) 

save.image('data24.rdata')


## grab baseline vital signs:
advs2=read_xpt(paste(indata2, "VS.xpt", sep='')) 
advs4=read_xpt(paste(indata4.sdtm, "VS.xpt", sep='')) 


vs1 = full_join(advs2, advs4) %>% filter(VISITDY==1)
save.image('data24.rdata')

## Check for duplicate lab and vital signs
table(table(vs1$USUBJID, vs1$VSTEST))
vs1$idtest=paste(vs1$USUBJID, vs1$VSTEST)

vs = vs1 %>% arrange(idtest, desc(VSDTC))
vs=vs[!duplicated(vs$idtest) & !is.na(vs$VSSTRESN),] %>% 
  mutate(test=paste(VSTEST, " (", VSSTRESU,")",sep=''))

subj.vitals = cast(vs, USUBJID~VSTEST, value='VSSTRESN')


## grab baseline lab values:
lb2=read_sas(paste(indata2,"adlb.sas7bdat",sep=''))%>% filter(AVISIT=="Baseline")
lb4=read_xpt(paste(indata4.adam,"adlb.xpt",sep=''))%>% filter(AVISIT=="Baseline")
save.image('data24.rdata')


labs0 = full_join(lb2, lb4) %>%
  filter(PARAMCD %in% c('ALT','AST','BILI','GFRE','EOS','EOSLE',
                        'HGB','WBC','LYM','LYMLE','NEUT','NEUTLE','PLAT')) 


original.labs = full_join(labs0, adsl %>% select(USUBJID,TRTFL,  cohort)) %>% 
#  filter(cohort %in% c('ACTT1', 'ACTT2', 'ACTT3') & TRTFL=='Y') %>%
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

## check
ggplot(data=original.labs %>% filter(cohort %in% c('ACTT3','ACTT4')))+
  geom_boxplot(aes(x=cohort, y=log(AVAL+.01),fill=cohort))+
  facet_wrap(~PARAMCD)

subj.labs = cast(labs1, USUBJID~PARAMCD, value='AVAL') %>%
  mutate(EOSLE = EOS/WBC,
         LYMLE = LYM/WBC,
         NEUTLE = NEUT/WBC)

## Recode ratios > 1 to missing as well as their components
neutle.ind = which(subj.labs$NEUTLE>1)
subj.labs$NEUTLE[neutle.ind]=NA
subj.labs$NEUT[neutle.ind]=NA
subj.labs$WBC[neutle.ind]=NA

which(subj.labs$EOSLE>1) ## NONE
which(subj.labs$LYMLE>1) ## NONE

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

load('actt123_mins')
mins123 = mins[-c(3,4,6,14,15),]
mins=labs %>% group_by(PARAM) %>% 
  summarise(min=min(AVAL_nozero, na.rm=TRUE))
table(mins$PARAM==mins123$PARAM)

## If min was lower in ACTT 1 or 2, grab that one: 
mins$min[mins123$min<mins$min]=mins123$min[mins123$min<mins$min]

labs.transformed=left_join(labs,mins)
labs.transformed$log.AVAL_nozero=log(labs.transformed$AVAL_nozero)

labs.transformed$log.AVAL_nozero[which(labs.transformed$AVAL==0)]=
  log(labs.transformed$min[which(labs.transformed$AVAL==0)]/2)


subj.labs2 = cast(labs.transformed, USUBJID~PARAMCD, value='log.AVAL_nozero') 
names(subj.labs2)[2:length(names(subj.labs2))]=
  paste("LOG.", names(subj.labs2)[2:length(names(subj.labs2))],sep='')

save.image('data24.rdata')

cov1 = full_join(adef,adsl)
cov2 = full_join(cov1, ho)
cov3 = full_join(cov2, subj.vitals)
cov4 = full_join(cov3, subj.labs)
cov5 = full_join(cov4, subj.labs2)
cov6 = full_join(cov5, consc)
cov7 = full_join(cov6, durations)

cov7$'Respiratory Rate'[which(cov7$BCSOSN==7)]=30
cov7$'Oxygen Saturation'[which(cov7$BCSOSN==7)]=90
cov7$ETHNIC[cov7$ETHNIC =='NOT REPORTED'] = 'NOT REPORTED/UNKNOWN'
cov7$ETHNIC[cov7$ETHNIC =='UNKNOWN'] = 'NOT REPORTED/UNKNOWN'

covariates = cov7 %>% select(-c(Weight, Height, BMI)) %>%
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
           (consciousness==3) + (BCRP>=50) + (BCRP>=100)) %>%
  mutate(news2 = (`Respiratory Rate`<=8)*3 +
           (`Respiratory Rate`>=9  & `Respiratory Rate`<=11)*1 +  
           (`Respiratory Rate`>=21  & `Respiratory Rate`<=24)*2 +  
           (`Respiratory Rate`>=25)*3 +  
           (`Oxygen Saturation`<=91)*3+ 
           (`Oxygen Saturation`>=92  & `Oxygen Saturation`<=93)*2 + 
           (`Oxygen Saturation`>=94  & `Oxygen Saturation`<=95)*1 + 
           (BCSOSN>4)*2+
           (`Temperature`<=35)*3+ 
           (`Temperature`>=35.1  & `Temperature`<=36)*1 + 
           (`Temperature`>=38.1  & `Temperature`<=39)*1 + 
           (`Temperature`>=39.1)*2 + 
           (`Systolic Blood Pressure`<=90)*3+
           (`Systolic Blood Pressure`>=91  & `Systolic Blood Pressure`<=100)*2 + 
           (`Systolic Blood Pressure`>=101  & `Systolic Blood Pressure`<=110)*1 + 
           (`Systolic Blood Pressure`>=220)*3 + 
           (`Heart Rate`<=40)*3+
           (`Heart Rate`>=41  & `Heart Rate`<=50)*1 + 
           (`Heart Rate`>=91  & `Heart Rate`<=110)*1 + 
           (`Heart Rate`>=111  & `Heart Rate`<=130)*2 + 
           (`Heart Rate`>=131)*3 + 
           consciousness)

covariates$news[covariates$phase=='4'] = covariates$news2[covariates$phase=='4'] 
covariates$symp.onset.to.hospitalization[covariates$symp.onset.to.hospitalization<0 & 
                                           !is.na(covariates$symp.onset.to.hospitalization)]=0

covariates$DEXMTHFL[covariates$cohort=='ACTT4']='Y'
covariates2=covariates %>% mutate(dex=!is.na(DEXMTHFL)&DEXMTHFL=='Y' )

# Censor people at Day 28
ttedat$CNSR[ttedat$AVAL>28]=1
ttedat$AVAL[ttedat$AVAL>28]=28

all=left_join(ttedat, covariates2) %>% 
  filter(PARAMCD=='TTRECOV') %>% 
  mutate(White=RACE=='WHITE', 
         severity67 = BCSOSN>=6, 
         site=factor(SITEID),
         recovered=1-CNSR) %>% 
        rename(tte_recov=AVAL) %>% select(-news2)

mort = ttedat %>% filter(PARAMCD=='D29DTHEF') %>%
  mutate(died=1-CNSR) %>%  rename(tte_death=AVAL) %>%
  select(USUBJID, died, tte_death)

analysis.data24 = left_join(all, mort) %>% 
 filter(cohort %in% c("ACTT2",'ACTT4') & TRTFL=='Y' & BCSOSN>4 & dex==1) 

table(analysis.data24$cohort)

setwd(outdir)
save(analysis.data24, file="analysis.data24.Rda")
save.image("actt24data.rdata")

