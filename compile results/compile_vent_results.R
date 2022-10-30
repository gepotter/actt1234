# Programmer Name: Gail Potter
# Date (month/year) program was first created:  NOV 2021
# R version:  4.0.4
# Compile bootstrap results from cluster, create output tables

###### Libraries.   
library(r2rtf)
library(dplyr)
library(survival)
library(haven)
library(CBPS)
library(cobalt)
library(scales)

rm(list=ls())

###### Directories 
outdir= 'N:/COVID-19/ACTT/Combined ACTT Analyses/ACTT123 comparison of RDV arm/output'
bootdir = 'N:/COVID-19/ACTT/Combined ACTT Analyses/ACTT123 comparison of RDV arm/output/cluster_results_ventilation'

or=NULL
seeds=5:9
for (s in 1:length(seeds)){
  load(paste(bootdir,'/boothrs_ventilation_50imputations_20_resamples_', seeds[s],sep=''))
  or=c(or,unlist(boothrs1))
}
length(or)
quantile(or, c(0.025, 0.975))
exp(quantile(or, c(0.025, 0.975)))

load(paste(outdir,'/dat12_with_weights', sep=''))
dat12.67 = dat12_with_weights %>% filter(BCSOSN %in% 6:7) %>% mutate(ventilated = BCSOSN==7)
mod.observed = CBPS(comparison.group.flag ~ male + Black	+ Asian	+	White + hisp + AGE	+ 
                      severity67 + news + imputed.mort4c +BDURSYMP + symp.onset.to.hospitalization + num.comorbidities+
                      HeartRate	+OxygenSaturation	+ RespiratoryRate	+	SystolicBloodPressure	+Temperature	+
                      LOG.ALT + LOG.AST + LOG.BILI + imputed.logcrp+  LOG.GFRE +LOG.EOS+ LOG.EOSLE+
                      LOG.HGB + LOG.WBC + LOG.LYM + LOG.LYMLE + LOG.NEUT + LOG.NEUTLE +
                      LOG.PLAT +(AGE>=50&AGE<60)+(AGE>=60&AGE<70)+(AGE>=70&AGE<80)+(AGE>=80), 
                    data = dat12.67, ATT=0, method='exact')

newnames=c(imputed.logcrp ='log(CRP)',
           imputed.mort4c='Modified 4c mortality scale',
           Temperature='Temperature',
           hisp='Hispanic', 
           LOG.LYMLE='log(lymphocyte/leukocyte ratio)', 
           imputed.logcrp='log(CRP)', 
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
           Asian='Asian',
           symp.onset.to.hospitalization='Symptom onset to hospitalization',
           LOG.WBC='log(White blood cells)',
           LOG.PLAT='log(platelets)',
           'AGE >= 50 & AGE < 60' ='Age 50-59',
           LOG.HGB='log(hemoglobin)',
           BDURSYMP='Symptom onset to enrollment')
setwd(outdir)
res=300
png("loveplot_ventilation.png", width=res*480/72, height=.9*res*480/72, res=300)
love.plot(mod.observed,drop.distance = TRUE, 
          var.order = "unadjusted", binary='std', 
          var.name=newnames,
          abs = TRUE, line = TRUE,thresholds = c(m = .1))
dev.off()

Pr.comparison.group = mod.observed$fitted.values
wt= as.numeric(dat12.67$comparison.group.flag==TRUE)*(1/Pr.comparison.group) + 
  as.numeric(dat12.67$comparison.group.flag==FALSE)*(1/(1-Pr.comparison.group))

pr.treat=mean(dat12.67$comparison.group.flag==1)
pr.control=1-pr.treat
pr.treat.vector = pr.treat*dat12.67$comparison.group.flag + 
  pr.control*(1-dat12.67$comparison.group.flag)
dat12.67$stabilized.weights = wt*pr.treat.vector
dat12.67$ventilated = dat12.67$BCSOSN==7

mod = glm(ventilated ~ cohort, weights=stabilized.weights, 
          data = dat12.67,family = binomial)
est = mod$coef[2]

unwtd.mod = glm(ventilated ~ cohort,  data = dat12.67,family = binomial)
summary(unwtd.mod)$coef

est.unwtd = summary(unwtd.mod)$coef[2,1]
pval.unwtd=summary(unwtd.mod)$coef[2,4]
ci.unwtd = confint(unwtd.mod )[2,]

unwtd.ventilation = c(exp(est.unwtd), exp(ci.unwtd), pval.unwtd)
round(unwtd.ventilation,4)

ci.lim=NULL
alphas=seq(.95,.9999,.0001)
for (a in 1:length(alphas) ) ci.lim[a] = quantile(or, alphas[a])
ci.lim
ci.lim.alpha  = alphas[which(ci.lim==max(ci.lim[ci.lim<=0]))]
pval=2*(1-ci.lim.alpha)
pval

wtd.ventilation = c(exp(est), exp(quantile(or, c(0.025, 0.975))), pval)
wtd.ventilation

tab=rbind(unwtd.ventilation, wtd.ventilation)
tab[, 1:3]=round(tab[, 1:3],2)
tab[, 4]=pvalue(as.numeric(tab[,4]))
tab
# tab[,4]=signif(tab[,4],3)
# tab[3:4,4]=as.character(round(tab[3:4,4],3))
colnames(tab)=c("Odds Ratio", "CI low", "CI up", "P-value")
rownames(tab)=c("Unweighted ventilation", "Weighted ventilation")

save(tab, file=paste(outdir,'/venttab',sep=''))


dat = data.frame(Method=rownames(tab), hr = tab[,1], ci=paste(tab[,2], ',', tab[,3]), pval=tab[,4])

dat %>% rtf_colheader(
    #colheader = " | ACTT-1 | ACTT-2 | ACTT-2 | ACTT-3"
    colheader = "Outcome | Odds Ratio | 95% Confidence Interval | P-value",
    text_font_size=11, col_rel_width = c(10,6,6,6)
  )  %>% rtf_body( text_font_size=11, col_rel_width = c(10,6,6,6))%>%
    rtf_page(orientation='portrait', margin = rep(0.5,6),
           border_first = "single", border_last = "single",nrow=60)%>%
  rtf_encode() %>%
  write_rtf(paste(outdir,"/vent12.rtf",sep=''))




