# Programmer: Gail Potter
# Date: NOV 2021
# Protocol: ACTT 1/2/3 combined analysis
# Is the odds of baseline intubation different in ACTT 2 vs ACTT 1?
# R version 4.0.4
# input: dat12_with_weights, imputation_model1, imputation_model2

library(dplyr)
library(haven)
library(CBPS)
library(parallel)
library(survey)
library(getopt)

cmd_options= commandArgs(TRUE);
seed=getopt(matrix(c('seed','s', 1, "integer"), byrow=TRUE, ncol=4))$seed


load('dat12_with_weights')
load("imputation_model1")
load("imputation_model2")
select=dplyr::select
dat12=dat12_with_weights %>% 
  filter(BCSOSN %in% 6:7) %>% mutate(ventilated = BCSOSN==7) %>% select(-wt)

####  Calculate uncertainty
nboot=20 ## number of bootstraps per task.   

n.imputations=50

detectBatchCPUs <- function() { 
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
  if (is.na(ncores)) { 
    ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
  } 
  if (is.na(ncores)) { 
    return(2)
  } 
  return(ncores) 
}

mc=detectBatchCPUs()
cl=makeCluster(mc)
clusterSetRNGStream(cl, seed)
clusterExport(cl, varlist=list('dat12', 'step1', 'step2', 
                               'nboot', 'n.imputations'))
  boothrs1 = clusterEvalQ(cl, {
    library(dplyr)
    library(CBPS)
    library(parallel)
    library(survey)
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
      
    bootfun = function(dat, nboot, n.imputations) {
      or.bootstrap=NULL
      
      for (b in 1:nboot){
        boot.indices = sample(1:nrow(dat), replace=TRUE)
        bootdat = dat[boot.indices,]
        
        or.imputed=NULL
        for (i in 1:n.imputations) {
          augmented.dat = rbind(augment(bootdat %>% filter(phase==1), mod=step1),
                                augment(bootdat %>% filter(phase==2), mod=step2)) %>% arrange(USUBJID)

          mod <- CBPS(comparison.group.flag ~ male + Black	+ Asian	+	White + hisp + AGE	+ 
                        severity67 + news + mort4c +BDURSYMP + symp.onset.to.hospitalization + num.comorbidities+
                        HeartRate	+OxygenSaturation	+ RespiratoryRate	+	SystolicBloodPressure	+Temperature	+
                        LOG.ALT + LOG.AST + LOG.BILI + logcrp+  LOG.GFRE +LOG.EOS+ LOG.EOSLE+
                        LOG.HGB + LOG.WBC + LOG.LYM + LOG.LYMLE + LOG.NEUT + LOG.NEUTLE +
                        LOG.PLAT +(AGE>=50&AGE<60)+(AGE>=60&AGE<70)+(AGE>=70&AGE<80)+(AGE>=80) , 
                      data = augmented.dat, ATT=0, method='exact') 
          Pr.comparison.group = mod$fitted.values
          wt= as.numeric(augmented.dat$comparison.group.flag==TRUE)*(1/Pr.comparison.group) + 
            as.numeric(augmented.dat$comparison.group.flag==FALSE)*(1/(1-Pr.comparison.group))
          
          pr.treat=mean(augmented.dat$comparison.group.flag==1)
          pr.control=1-pr.treat
          pr.treat.vector = pr.treat*augmented.dat$comparison.group.flag + 
            pr.control*(1-augmented.dat$comparison.group.flag)
          augmented.dat$stabilized.weights=wt*pr.treat.vector
          
          mod = glm(ventilated ~ cohort, weights=stabilized.weights, 
                    family = binomial, data=augmented.dat)
          or.imputed[i] = mod$coef[2]
        }
        or.bootstrap[b] = mean(or.imputed)
      }
      return(or.bootstrap)
    }
    
    bootfun(dat12, nboot=nboot, n.imputations=n.imputations)
  })
stopCluster(cl)

save(boothrs1, file=
  paste('boothrs_ventilation', '_', n.imputations,'imputations_',nboot, '_resamples','_',seed, sep=''))
