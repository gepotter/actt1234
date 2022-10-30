# Programmer: Gail Potter
# Date: NOV 2021
# Protocol: ACTT 1/2/3 combined analysis
# Survival curves for ACTT 1 vs. 2
# R version 4.0.4
# Input: dat23_with_weights, actt23hrtab

library(dplyr)
library(ggplotify)
library(ggplot2)
library(ggbeeswarm)
library(CBPS)
library(Hmisc)
library(cobalt)
library(RISCA)
library(ggpubr)
library(patchwork)

rm(list=ls())
outdir="N:/COVID-19/ACTT/Combined ACTT Analyses/ACTT123 comparison of RDV arm/output"
setwd(outdir)

load(paste(outdir,'/dat23_with_weights', sep=''))
load(paste(outdir,'/actt23hrtab',sep=''))

dat23 = dat23_with_weights 

get.curves = function(tte, cnsr, phase, wt) {
  unweighted.curves=ipw.survival(times=tte, 
                                 failures=(1-cnsr), 
                                 variable=phase, 
                                 weights=rep(1, length(wt)))
  
  weighted.curves=ipw.survival(times=tte, 
                               failures=(1-cnsr), 
                               variable=phase, 
                               weights=wt)
  
  
  curvdat = data.frame(days=c(weighted.curves$table.surv$times, unweighted.curves$table.surv$times),
                       survival=c(weighted.curves$table.surv$survival, unweighted.curves$table.surv$survival),
                       Method=factor(c(rep('Weighted',nrow(weighted.curves$table.surv)),
                                       rep('Unweighted', nrow(unweighted.curves$table.surv)))),
                       Stage=factor(paste("ACTT-",c(weighted.curves$table.surv$variable, 
                                                    unweighted.curves$table.surv$variable), sep='')))
  return(curvdat)
}

recov=get.curves(tte=dat23$tte_recov, cnsr=(1-dat23$recovered), phase=dat23$phase, 
                wt=dat23$wt)
mort=get.curves(tte=dat23$tte_death, cnsr=(1-dat23$died), 
                phase=dat23$phase, wt=dat23$wt)

leg=as.ggplot(get_legend(ggplot(data=mort, aes(x=days, y=1-survival, 
                                    group=interaction(Method, Stage),
                                    colour=Stage, 
                                    linetype=Method))+
  geom_line(size=1)+ylim(0,.1)+
    theme_classic()+
  scale_colour_manual(values=c("navyblue", "red"))+xlab("Study Day")+
  scale_linetype_manual(values=c("dotted", "solid"))+
theme(legend.position='bottom')))
  
# res=300
# setwd(outdir)
# png("ipw_curves23_both.png", width=1.2*res*480/72, height=.6*res*480/72, res=300)
# (p.recov.both | p.mort.both) / leg + plot_layout(heights=c(10,2))
# dev.off()


hr.recov = paste('Hazard ratio ', tab[2,1], " [", tab[2,2], ', ', tab[2,3], ']',sep='')
pval.recov = paste('p=',tab[2,4],sep='')
hr.mort = paste('Hazard ratio ',tab[4,1], " [", 
                format(as.numeric(tab[4,2]), nsmall=2), ', ', tab[4,3], ']',sep='')
pval.mort  = paste('p=',tab[4,4],sep='')


wtd.recov = ggplot(data=recov %>% filter(Method=='Weighted'), 
       aes(x=days, y=1-survival, 
           group=Stage,colour=Stage))+
  geom_line(size=1)+ylim(0,1)+theme_classic()+
  theme(text=element_text(family="serif"))+
  geom_text(x=18, y=.2, label=hr.recov, col='black', family='serif')+
  geom_text(x=18, y=.125, label=pval.recov, col='black', family='serif')+
  geom_text(x=18, y=.6, label='ACTT-1', col='navyblue', family='serif')+
  geom_text(x=11, y=.78, label='ACTT-2', col='red', family='serif')+
  theme(legend.position='none')+
  scale_colour_manual(values=c("navyblue", "red"))+xlab("Study Day")+
  scale_x_continuous(breaks=c(0,7,14,21,28), limits = c(0, 28))+
  ylab('Probability of recovery')+
  ggtitle('Time to recovery')


wtd.mort = ggplot(data=mort %>% filter(Method=='Weighted'), 
                     aes(x=days, y=1-survival, 
                     group=Stage,colour=Stage))+
  geom_line(size=1)+ylim(0,.105)+
  theme_classic()+
  theme(text=element_text(family="serif"))+
  theme(legend.position='none')+
  scale_colour_manual(values=c("navyblue", "red"))+
  xlab("Study Day")+
  geom_text(x=12, y=.09, label=hr.mort, col='black', family='serif')+
  geom_text(x=12, y=.0825, label=pval.mort, col='black', family='serif')+
  geom_text(x=8, y=.05, label='ACTT-1', col='navyblue', family='serif')+
  geom_text(x=14, y=.022, label='ACTT-2', col='red', family='serif')+
  scale_x_continuous(breaks=c(0,7,14,21,28), limits = c(0, 28))+
  ylab('Probability of death')+
  ggtitle('Time to death')

# res=300
# setwd(outdir)
# png("ipw_curves12.png", width=1.2*res*480/72, height=.6*res*480/72, res=300)
# wtd.recov | wtd.mort
# dev.off()

p.recov.both = ggplot(data=recov, aes(x=days, y=1-survival, 
                                      group=interaction(Method, Stage),
                                      colour=Stage, 
                                      linetype=Method))+
  geom_line(size=1)+ylim(0,1)+theme_classic()+
  theme(legend.position='none',
        text=element_text(family='serif'))+
  geom_text(x=14, y=.62, label='ACTT-2', col='navyblue', family='serif')+
  geom_text(x=7, y=.85, label='ACTT-3', col='red', family='serif')+
  scale_colour_manual(values=c("navyblue", "red"))+xlab("Study Day")+
  scale_linetype_manual(values=c("dotted", "solid"))+
  scale_x_continuous(breaks=c(0,7,14,21,28), limits = c(0, 28))+
  ylab('Probability of recovery')+
  geom_segment(x=16, y=.15, xend=18, yend=.15, linetype='dotted', colour='black',size=1)+
  geom_text(x=19, y=.15, label='Unweighted', hjust=0, col='black', family='serif')+
  geom_segment(x=16, y=.05, xend=18, yend=.05,  colour='black',size=1)+
  geom_text(x=19, y=.05, label='Weighted', hjust=0, col='black', family='serif')+
  ggtitle('Time to recovery')


estci= format(as.numeric(tab[1,1:3], nsmall=2))
pval=tab[1,4]#format(as.numeric(tab[1,4], nsmall=3))
hr.recov.unadj = c('Unweighted',estci[1], paste("[", estci[2], ', ', estci[3], ']', sep=''),pval)

estci= format(as.numeric(tab[2,1:3], nsmall=2))
pval=tab[2,4]#format(as.numeric(tab[2,4], nsmall=3))
hr.recov.adj = c('Weighted',estci[1], paste("[", estci[2], ', ', estci[3], ']', sep=''),pval)

hr.recov.unadj
hr.recov.adj


estci= c(tab[3,1:2], format(as.numeric(tab[3,3], nsmall=2)))
pval=format(as.numeric(tab[3,4], nsmall=3))
hr.mort.unadj = c('Unweighted',estci[1], paste("[", estci[2], ', ', estci[3], ']', sep=''),pval)

estci= format(as.numeric(tab[4,1:3], nsmall=2))
pval=format(as.numeric(tab[4,4], nsmall=3))
hr.mort.adj = c('Weighted',estci[1], paste("[", estci[2], ', ', estci[3], ']', sep=''),pval)

hr.mort.unadj
hr.mort.adj

p.mort.both = ggplot(data=mort, aes(x=days, y=1-survival, 
                                    group=interaction(Method, Stage),
                                    colour=Stage, 
                                    linetype=Method))+
  # geom_hline(yintercept=.1, col='gray')+
  # geom_hline(yintercept=.05, col='gray')+
  geom_line(size=1)+
  ylim(0,.3)+
  theme_classic()+
  theme(legend.position='none',
        text=element_text(family='serif'))+
  scale_colour_manual(values=c("navyblue", "red"))+xlab("Study Day")+
  scale_linetype_manual(values=c("dotted", "solid"))+
  scale_x_continuous(breaks=c(0,7,14,21,28), limits = c(0, 28))+
  ylab('Probability of death')+
  geom_text(x=8, y=.035, label='ACTT-2', col='navyblue', family='serif')+
  geom_text(x=14, y=.0, label='ACTT-3', col='red', family='serif')+
  geom_segment(x=14, y=.3, xend=16, yend=.3, linetype='dotted', colour='black',size=1)+
  geom_text(x=17, y=.3, label='Unweighted', hjust=0, col='black', family='serif')+
  geom_segment(x=14, y=.27, xend=16, yend=.27,  colour='black',size=1)+
  geom_text(x=17, y=.27, label='Weighted', hjust=0, col='black', family='serif')+
  ggtitle('Time to death')



tab.recov=rbind(c("Model", "Hazard ratio", "95% CI", "P-value"),
                hr.recov.unadj,hr.recov.adj)
tab.recov=rbind(c("Model", "Hazard ratio", "95% CI", "P-value"),
                hr.recov.unadj,hr.recov.adj)[,1:3]
tab.recov

nrow=nrow(tab.recov)
ncol=ncol(tab.recov)
row.locations = c(1, 2.25, 3.5)
tabdat <- data.frame(V0 = rep(nrow:1,ncol), 
                      V05 = rep(row.locations,each=nrow),
                      V1 = as.vector(tab.recov))

tablegrob = ggplot(tabdat,  aes(x = V05, y = V0, label = V1))+
  geom_text(size = 4.5, hjust=0.5, vjust=0.5,family='serif') + 
  theme_classic() +
  geom_hline(aes(yintercept=3.5)) + 
  geom_hline(aes(yintercept=0.5)) + 
  theme(panel.grid.major = element_blank(), 
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(colour="white"),
        axis.text.y = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_line(colour="white")) + 
  labs(x="",y="") +
  coord_cartesian(xlim=c(0.5,4), ylim=c(-.5,7))+
  annotate("text", x=1.5,y=4.5,label='Recovery Estimates', family='serif',size=5)


tablegrob

tab.mort=rbind(c("Model", "Hazard ratio", "95% CI", "P-value"),
               hr.mort.unadj,hr.mort.adj)[,1:3]
rownames(tab.mort)=NULL
tab.mort

nrow=nrow(tab.mort)
ncol=ncol(tab.mort)
row.locations = c(1, 2.25, 3.5)
tabdat2 <- data.frame(V0 = rep(nrow:1,ncol), 
                      V05 = rep(row.locations,each=nrow),
                      V1 = as.vector(tab.mort))
tabdat2

tablegrob2 = ggplot(tabdat2, aes(x = V05, y = V0, label = V1))+
  geom_text(size = 4.5, hjust=.5, vjust=0.5,family='serif') + 
  theme_classic() +
  geom_hline(aes(yintercept=3.5)) + 
  geom_hline(aes(yintercept=0.5)) + 
  theme(panel.grid.major = element_blank(), 
        legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x = element_text(colour="white"),
        axis.text.y = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks = element_line(colour="white")) + 
  labs(x="",y="") +
  coord_cartesian(xlim=c(0.5,4), ylim=c(-.5,7))+
  annotate("text", x=1.5,y=4.5,label='Mortality Estimates', family='serif',size=5)

tablegrob2


res=300
setwd(outdir)
png("ipw_curves23_both_slides.png", width=1.3*res*480/72, height=.8*res*480/72, res=300)
(p.recov.both + p.mort.both) / (tablegrob + tablegrob2) + plot_layout(heights=c(10, 7.5))
dev.off()


