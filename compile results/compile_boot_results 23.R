# Programmer Name: Gail Potter
# Date (month/year) program was first created:  NOV 2021
# R version:  4.0.4
# Compile bootstrap results from cluster, create output tables

###### Libraries.   
library(r2rtf)
library(dplyr)
library(survival)
library(haven)
library(scales)

###### Input directories 
rm(list=ls())
indir1='N:/COVID-19/ACTT/Combined ACTT Analyses/ACTT123 comparison of RDV arm/output/cluster_results_recovery/ACTT 2 vs 3'
indir2='N:/COVID-19/ACTT/Combined ACTT Analyses/ACTT123 comparison of RDV arm/output/cluster_results_mortality/ACTT 2 vs 3'

###### Output directories 
outdir= 'N:/COVID-19/ACTT/Combined ACTT Analyses/ACTT123 comparison of RDV arm/output'


setwd(indir1)
loghr=NULL
seeds=5:9
for (s in 1:length(seeds)){
  load(paste(indir1,'/boothrs_recovery23_50imputations_20_resamples_', seeds[s],sep=''))
  loghr=c(loghr,unlist(boothrs1))
}
length(loghr)
quantile(loghr[1:5000], c(0.025, 0.975))
exp(quantile(loghr[1:5000], c(0.025, 0.975)))


## Get parameter estimate from mean-imputed data
load( paste(outdir,'/dat23_with_weights', sep=''))
fit2=coxph(Surv(tte_recov, recovered) ~ comparison.group.flag+ cluster(USUBJID), 
           data = dat23_with_weights, 
           weights=dat23_with_weights$wt)
est=fit2$coef

## bootstrap CI values
ci.low = quantile(loghr[1:5000], 0.025)
ci.up = quantile(loghr[1:5000], 0.975)
z=est/sd(loghr)

## bootstrap p-value
pval= 2*pnorm(abs(z), lower.tail=FALSE)



wtd.recovery = c(exp(est), exp(ci.low), exp(ci.up), pval)
wtd.recovery

fit.unweighted =coxph(Surv(tte_recov, recovered) ~ comparison.group.flag+ cluster(USUBJID), 
           data = dat23_with_weights)
est=fit.unweighted$coef
se=abs(est)/sqrt(fit.unweighted$score)


## score CI 
score.results = read_sas(paste(outdir, 'scoreci_recovery.sas7bdat', sep='/'))
ci.low = as.numeric(score.results[2,4])
ci.up = as.numeric(score.results[2,5])

# ci.low = est - 1.96*se
# ci.up = est + 1.96*se 
pval =  summary(fit.unweighted)$sctest[3]
unwtd.recovery = c(exp(est), exp(ci.low), exp(ci.up), pval)



###########    MORTALITY 
seeds= 5:9
loghr=NULL
for (s in 1:length(seeds)){
  load(paste(indir2, '/boothrs_mortality_stabilized_catage_50imputations_20_resamples_', seeds[s],sep=''))
  loghr=c(loghr,unlist(boothrs1))
}
length(loghr)
exp(quantile(loghr, c(0.025, 0.975)))
quantile(exp(loghr), c(0.025, 0.975))


## Get estimate from mean-imputed data
fit2=coxph(Surv(tte_death, died) ~ comparison.group.flag+ cluster(USUBJID), 
           data = dat23_with_weights, 
           weights=dat23_with_weights$wt)
est=fit2$coef
ci.low = quantile(loghr, 0.025)
ci.up = quantile(loghr, 0.975)
z=est/sd(loghr)
# pval= 2*pnorm(abs(z), lower.tail=FALSE)

ci.lim=NULL
alphas=seq(.94,.99,.00001)
for (a in 1: length(alphas) ) ci.lim[a] = quantile(loghr, alphas[a])
min(ci.lim[ci.lim>=0])
ci.lim.alpha  = alphas[which(ci.lim==min(ci.lim[ci.lim>=0]))]
pval=2*(1-ci.lim.alpha)
pval

wtd.mortality = c(exp(est), exp(ci.low), exp(ci.up), pval)

fit.unweighted =coxph(Surv(tte_death, died) ~ comparison.group.flag,
                      data = dat23_with_weights)
est=fit.unweighted$coef
# se=abs(est)/sqrt(fit.unweighted$score)
# ci.low = est - 1.96*se
# ci.up = est + 1.96*se 

score.results = read_sas(paste(outdir, 'scoreci_mortality.sas7bdat', sep='/'))
ci.low = as.numeric(score.results[2,4])
ci.up = as.numeric(score.results[2,5])


pval =  summary(fit.unweighted)$sctest[3]
unwtd.mortality = as.numeric(c(exp(est), exp(ci.low), exp(ci.up), pval))

tab=rbind(unwtd.recovery, wtd.recovery, unwtd.mortality,wtd.mortality)
tab[c(1,2,4), 1:3]=round(tab[c(1,2,4), 1:3],2)
tab[3, 1:2]=round(tab[3, 1:2],2)
tab[3, 3]=round(tab[3, 3],3)
tab[,4]=pvalue(as.numeric(tab[,4]),accuracy=0.001)
tab[3:4,4]=as.character(round(tab[3:4,4],3))
colnames(tab)=c("Hazard Ratio", "CI low", "CI up", "P-value")
rownames(tab)=c("Unweighted recovery", "Weighted recovery", "Unweighted mortality",
                "Weighted mortality")
tab
save(tab, file=paste(outdir,'/actt23hrtab',sep=''))


dat = data.frame(Method=rownames(tab), hr = tab[,1], ci=paste(tab[,2], ',', tab[,3]), pval=tab[,4])

dat%>%
  rtf_colheader(
    #colheader = " | ACTT-1 | ACTT-2 | ACTT-2 | ACTT-3"
    colheader = "Outcome | Hazard Ratio | 95% Confidence Interval | P-value",
    text_font_size=11, col_rel_width = c(10,6,6,6)
  )  %>% rtf_body( text_font_size=11, col_rel_width = c(10,6,6,6))%>%
    rtf_page(orientation='portrait', margin = rep(0.5,6),
           border_first = "single", border_last = "single",nrow=60)%>%
  rtf_encode() %>%
  write_rtf(paste(outdir,"/results23.rtf",sep=''))




