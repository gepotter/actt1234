# Programmer: Gail Potter
# Date: NOV 2021
# Protocol: ACTT 1/2/3 combined analysis
# R version 4.0.4
# input: analysis.data.Rda, imputation_model2, imputation_model3


#library(twang)
library(dplyr)
library(stringr)
library(haven)
library(CBPS)
library(Hmisc)
library(cobalt)
library(tableone)
library(RISCA)

rm(list=ls())
indata1="N:/COVID-19/ACTT/ACTT-1/Data/"
indata2="N:/COVID-19/ACTT/ACTT-2/Data/"
indata3="N:/COVID-19/ACTT/ACTT-3/Data/2021_02_24_Lock/"

outdir="N:/COVID-19/ACTT/Combined ACTT Analyses/ACTT123 comparison of RDV arm/output"
setwd(outdir)
load("analysis.data.Rda")
load("imputation_model2")
load("imputation_model3")

select=dplyr::select


dat=analysis.data %>% filter(PARAMCD=='TTRECOV') %>%mutate(White=RACE=='WHITE')
dat$severity67 = dat$BCSOSN>=6
dat$site=factor(dat$SITEID)

names(dat) =  str_remove(names(dat), " ")
names(dat) =  str_remove(names(dat), ", ")
names(dat) =  str_remove(names(dat), "/")
names(dat) =  str_remove(names(dat), " ")
names(dat) =  str_remove(names(dat), " ")
names(dat) =  str_remove(names(dat), " ")
names(dat) = unlist(lapply(strsplit(names(dat), split='(', fixed=TRUE), '[[', 1))


select=dplyr::select

actt23.pscore.vars = dat %>% 
  filter(cohort %in% c("ACTT3","ACTT2") & BCSOSN!=7 & !CLDFL) %>% 
  select(male, Black, Asian, hisp, AGE, BCSOSN, news, symp.onset.to.hospitalization, BDURSYMP,
         num.comorbidities, mort4c, HeartRate, OxygenSaturation, RespiratoryRate,	
         SystolicBloodPressure, Temperature, logcrp, LOG.ALT , LOG.AST , LOG.BILI , LOG.EOS, LOG.EOSLE,
         LOG.GFRE , LOG.HGB , LOG.LYM , LOG.LYMLE , LOG.NEUT , LOG.NEUTLE ,
         LOG.PLAT, cohort)

sum(complete.cases(actt23.pscore.vars))
mean(complete.cases(actt23.pscore.vars))


## function to augment a single data set

augment = function(dat, mod) {
  miss.ind=which(is.na(dat$logcrp)) 
  dat$logcrp[miss.ind] = rnorm(n=length(miss.ind),
                               mean=predict(mod, dat[miss.ind,]), sd=sd(mod$resid) )
  dat$BCRP[miss.ind]=exp(dat$logcrp[miss.ind])
  return(dat %>% mutate(  mort4c = 2*(AGE>=50)+2*(AGE>=60)+2*(AGE>=70)+1*(AGE>=80) + male + 
                            (num.comorbidities==1) +
                            (num.comorbidities>1) + (RespiratoryRate >=20 | BCSOSN==7) + 
                            (RespiratoryRate >=30| BCSOSN==7) +
                            2*(OxygenSaturation<92 | BCSOSN>4) + 
                            (consciousness==3) + (BCRP>=50) + (BCRP>=100)) )
}


augment.meanvalue = function(dat, mod) {
  miss.ind=which(is.na(dat$logcrp)) 
  dat$logcrp[miss.ind] = predict(mod, dat[miss.ind,])
  dat$BCRP[miss.ind]=exp(dat$logcrp[miss.ind])
  return(dat %>% mutate(  mort4c = 2*(AGE>=50)+2*(AGE>=60)+2*(AGE>=70)+1*(AGE>=80) + male + 
                            (num.comorbidities==1) +
                            (num.comorbidities>1) + (RespiratoryRate >=20 | BCSOSN==7) + 
                            (RespiratoryRate >=30| BCSOSN==7) +
                            2*(OxygenSaturation<92 | BCSOSN>4) + 
                            (consciousness==3) + (BCRP>=50) + (BCRP>=100)) )
}

dat23.all = dat %>% filter(cohort %in% c("ACTT3","ACTT2") & BCSOSN!=7 & !CLDFL) %>%
  select(-c("DSCHDT"  ,"CREAT"     ,"LOG.CREAT"))
    ## These variables were not analyzed and they have some missing data

## grab analyzable cases only- those with complete obs for everything but CRP
dat23 = dat23.all[complete.cases(dat23.all%>% 
                                   select(male, Black, Asian, White, hisp, AGE, BCSOSN, news, symp.onset.to.hospitalization, BDURSYMP,
                                          num.comorbidities, HeartRate, OxygenSaturation, RespiratoryRate,	
                                          SystolicBloodPressure, Temperature, LOG.ALT , LOG.AST , LOG.BILI , LOG.EOS, LOG.EOSLE,
                                          LOG.GFRE , LOG.HGB , LOG.LYM , LOG.LYMLE , LOG.NEUT , LOG.NEUTLE ,
                                          LOG.PLAT, cohort)),] ## MORT4C AND CRP EXCLUDED SINCE WE'LL IMPUTE THESE
dim(dat23)
dim(dat23.all)

dat23$comparison.group.flag = dat23$cohort=='ACTT3'

write.csv(dat23 %>% select(AVAL, CNSR, cohort), file=paste(outdir, '/dat23_recovery.csv', sep=''))


dat23_meanvalue_imputed = 
  full_join(augment.meanvalue(dat23 %>% filter(phase=='2'), mod=step2),
            augment.meanvalue(dat23 %>% filter(phase=='3'), mod=step3))%>%
  arrange(USUBJID) %>%
  mutate(flag80=AGE>=80,
         flag70=AGE>=70 & AGE<80,
         flag60=AGE>=60 & AGE<70,
         flag50=AGE>=50 & AGE<60)

mod.observed <- CBPS(comparison.group.flag ~ male + Black	+ Asian	+	White + hisp + AGE	+ 
                       severity67 + news + mort4c +BDURSYMP + symp.onset.to.hospitalization + num.comorbidities+
                       HeartRate	+OxygenSaturation	+ RespiratoryRate	+	SystolicBloodPressure	+Temperature	+
                       LOG.ALT + LOG.AST + LOG.BILI + logcrp+  LOG.GFRE +LOG.EOS+ LOG.EOSLE+
                       LOG.HGB + LOG.WBC + LOG.LYM + LOG.LYMLE + LOG.NEUT + LOG.NEUTLE +
                       LOG.PLAT +(AGE>=50&AGE<60)+(AGE>=60&AGE<70)+(AGE>=70&AGE<80)+(AGE>=80), 
                     data = dat23_meanvalue_imputed, ATT=0, method='exact')

mod.interaction <- CBPS(comparison.group.flag ~ male + Black	+ Asian	+	White + hisp + AGE	+
                          BDURSYMP + symp.onset.to.hospitalization + num.comorbidities+
                          factor(severity67) *(news + mort4c +HeartRate	+OxygenSaturation	+ RespiratoryRate	+	SystolicBloodPressure	+Temperature	+
                          LOG.ALT + LOG.AST + LOG.BILI + logcrp+  LOG.GFRE +LOG.EOS+ LOG.EOSLE+
                          LOG.HGB + LOG.WBC + LOG.LYM + LOG.LYMLE + LOG.NEUT + LOG.NEUTLE +
                          LOG.PLAT +(AGE>=50&AGE<60)+(AGE>=60&AGE<70)+(AGE>=70&AGE<80)+(AGE>=80)), 
                        data = dat23_meanvalue_imputed, ATT=0, method='exact')

test.stat = mod.observed$deviance-mod.interaction$deviance
df=length(mod.interaction$coef)-length(mod.observed$coef)
pchisq(test.stat, df=df, lower.tail=FALSE)

Pr.comparison.group = mod.observed$fitted.values
wt= as.numeric(dat23_meanvalue_imputed$comparison.group.flag==TRUE)*(1/Pr.comparison.group) + 
  as.numeric(dat23_meanvalue_imputed$comparison.group.flag==FALSE)*(1/(1-Pr.comparison.group))
comparison.group.sum = sum(wt[dat23_meanvalue_imputed$comparison.group.flag])
reference.group.sum  = sum(wt[!dat23_meanvalue_imputed$comparison.group.flag])


pr.treat=mean(dat23_meanvalue_imputed$comparison.group.flag==1)
pr.control=1-pr.treat
pr.treat.vector = pr.treat*dat23_meanvalue_imputed$comparison.group.flag + 
  pr.control*(1-dat23_meanvalue_imputed$comparison.group.flag)
stabilized.weights = wt*pr.treat.vector

dat23_meanvalue_imputed$wt = stabilized.weights

## Add mortality data
mort = analysis.data %>% filter(PARAMCD=='D29DTHEF') %>%
  mutate(died=1-CNSR) %>%  rename(tte_death=AVAL) %>%
  filter(cohort %in% c('ACTT2','ACTT3'))%>%
  select(USUBJID, died, tte_death, logcrp, mort4c, BCRP)



################  Identify sites included in both ACTT 2 and ACTT3
# for sensitivity analysis
site2=as.character(dat$site[dat$cohort=='ACTT2'])
site3=as.character(dat$site[dat$cohort=='ACTT3'])
sites23 = unique(site2[which(site2%in% site3)])


dat23_with_weights = left_join(dat23_meanvalue_imputed %>% 
  rename(imputed.logcrp = logcrp, imputed.mort4c=mort4c,imputed.BCRP=BCRP) %>%
  mutate(recovered=1-CNSR) %>% 
  rename(tte_recov=AVAL), mort)%>%
  mutate(egfr.excl = GFRE < 30,
         neut.excl = abs(NEUT) < 1,
         lymph.excl = abs(LYM) < .20,
         wbc.excl = WBC < 1.5,
         plat.excl = PLAT < 50,
         lab.exclusions = egfr.excl | neut.excl | lymph.excl | wbc.excl | plat.excl,
         labflag  = !lab.exclusions,
         siteflag = site %in% sites23)

newnames=c(severity67="Ordinal score 6 or 7",
Temperature='Temperature',
hisp='Hispanic', 
LOG.LYMLE='log(lymphocyte/leukocyte ratio)', 
logcrp='log(CRP)', 
LOG.LYM='log(lymphocytes)',
LOG.EOSLE='log(eosinophil/leukocyte ratio)', 
LOG.EOS='log(eosinophils)', 
num.comorbidities='Number of comorbidities',
RespiratoryRate='Respiratory Rate',
HeartRate='Heart Rate', 
AGE='Age', 
LOG.NEUTLE='log(Neutrophil/leukocyte ratio', 
White='White', 
LOG.GFRE='log(eGFR)', 
male='Male', 
news='National Early Warning Score',
OxygenSaturation='Oxygen Saturation',
LOG.AST='log(AST)',
LOG.NEUT='log(neutrophils)',
'AGE >= 80'='Age 80 or higher',
'AGE >= 60 & AGE < 70'='Age 60-69',
LOG.BILI='log(bilirubin)',
SystolicBloodPressure='SystolicBloodPressure',
'AGE >= 70 & AGE < 80'='Age 70-79',
LOG.ALT='log(ALT)',
Black='Black',
mort4c='Modified 4c mortality scale',
Asian='Asian',
symp.onset.to.hospitalization='Symptom onset to hospitalization',
LOG.WBC='log(White blood cells)',
LOG.PLAT='log(platelets)',
'AGE >= 50 & AGE < 60' ='Age 50-59',
LOG.HGB='log(hemoglobin)',
BDURSYMP='Symptom onset to enrollment')


love.plot(mod.observed,drop.distance = TRUE, 
          var.order = "unadjusted", binary='std', 
          var.name=newnames,
          abs = TRUE, line = TRUE,thresholds = c(m = .1))

