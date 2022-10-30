# Programmer Name: Gail Potter
# JUL 2021
# R version:  4.0.4
# Purpose: imputation model for missing CRP
# Input data sets: analysis.data.Rda
# Output:  
# Relies on: analysis data sets.R (don't need to run this, just need to load analysis.data.Rda) 
# History:   Updated 08APR2022 to remove individual comorbidity flags

library(dplyr)
library(stringr)
library(MASS)
library(ggplot2)
library(scales)
library(gridExtra)

rm(list=ls())
outdir="N:/COVID-19/ACTT/Combined ACTT Analyses/ACTT123 comparison of RDV arm/output"
setwd(outdir)
load("analysis.data34.Rda")

analysis.data= analysis.data34
analysis.data$severity67 = analysis.data$BCSOSN>=6
analysis.data$site=factor(analysis.data$SITEID)

names(analysis.data) =  str_remove(names(analysis.data), " ")
names(analysis.data) =  str_remove(names(analysis.data), ", ")
names(analysis.data) =  str_remove(names(analysis.data), "/")
names(analysis.data) =  str_remove(names(analysis.data), " ")
names(analysis.data) =  str_remove(names(analysis.data), " ")
names(analysis.data) =  str_remove(names(analysis.data), " ")
names(analysis.data) = unlist(lapply(strsplit(names(analysis.data), split='(', fixed=TRUE), '[[', 1))

select=dplyr::select

subj.data = analysis.data%>% filter(!duplicated(USUBJID))%>%  select(logcrp,site, AGE,RACE,male,ETHNIC,BMIFL,DIABFL,phase,
                                        CKDFL,CRDFL,ASTHMAFL,BDURSYMP,CORFL,HYPFL,CADFL,
                                        CHFFL,CANCERFL,IMMDFL,num.comorbidities,news,HeartRate ,OxygenSaturation,RespiratoryRate,
                                        SystolicBloodPressure,Temperature,consciousness,BCSOSN, 
                                        LOG.ALT,LOG.AST,LOG.BILI,LOG.EOS,LOG.EOSLE,LOG.GFRE,LOG.HGB,
                                        LOG.LYM,LOG.LYMLE,LOG.NEUT,LOG.NEUTLE,LOG.PLAT,LOG.WBC,
                                        ALT,AST,BILI,EOS,EOSLE,GFRE,HGB,LYM,LYMLE,NEUT,NEUTLE,PLAT,WBC,
                                        symp.onset.to.hospitalization,severity67,cohort) 

dat4 = subj.data %>% filter(cohort=="ACTT4") %>% select(-cohort)

nmiss=NULL
for (c in 1:ncol(subj.data))nmiss[c]  = sum(is.na(subj.data[,c]))


## Fit separate model for each cohort
dat4.cc = dat4[complete.cases(dat4),]
mean(complete.cases(dat4))

### REMOVED CLDFL FOR ACTT 4
formula = 'logcrp~male+RACE+ETHNIC+AGE+factor(BCSOSN)+
  news+BDURSYMP+symp.onset.to.hospitalization+ 
  BMIFL+DIABFL+CKDFL+CRDFL+ASTHMAFL+CORFL+CADFL+HYPFL+
  CHFFL+CANCERFL+IMMDFL+
  HeartRate +OxygenSaturation+RespiratoryRate+
  SystolicBloodPressure+Temperature+
  LOG.ALT+LOG.AST+LOG.BILI+LOG.GFRE+sqrt(EOS)+sqrt(EOSLE)+
  HGB+ LOG.WBC+ LOG.LYM+
  LOG.LYMLE+LOG.NEUT+NEUTLE+LOG.PLAT'

mod4=lm(formula=formula,  data=dat4.cc)
mod4
# Based on inspection of scatter plots, a square root transformation was 
# used for eosinophils and eosinophil/leukocyte ratio, no transformation 
# for hemoglobin or neutrophils, and a natural log transformation for all other labs.

step4=step(mod4)
AIC(mod4)
AIC(step4)

plot(step4$coef)
plot(step4$coef[-1])

save(step4, file=paste(outdir,'/imputation_model4',sep=''))

