# In this script we will re-analyse the data behind the
# recent GARD study with a certain question I have been
# curious about for a long time - dose RSI capture enough of
# the variance in outcome such that a dose-response becomes
# visible empirically.

require(survival)
require(prodlim)
require(ggplot2)

load("data_df.rda")
summary(rsi.all)
# Notice one patient, see later, has a negative survival time
# We shall discuss data quality and other statsitical
# analysis issues in another analysis
summary(factor(rsi.all$Site))

sum(rsi.all$Time_OS<0,na.rm=T)

# We only want all positive survival times and only those that had RT
# for our analysis
rsi.all<-rsi.all[is.na(rsi.all$Time_OS)==F &
                   rsi.all$Time_OS>0 & rsi.all$Received_RT==1,]

# Include GARD 
alpha=log(rsi.all$RSI)/(-2) - 0.05*2
rsi.all$gard=rsi.all$n*rsi.all$d*(alpha + 0.05*rsi.all$d)

# Explore each tumour site
summary(factor(rsi.all$Site))
# What is the event rate like in each study
summary(factor(rsi.all$Event_OS[rsi.all$Site=="breast_TN"]))
# too few events - think about this how comfortable do you 
# feel about 9 events! 
km0<-prodlim(Hist(Time_OS,Event_OS)~1,data=rsi.all[rsi.all$Site=="breast_TN",])
plot(km0,confint=T,percent=F,atrisk.at=seq(0,12,2),background.horizontal=NULL,
     axis1.at=seq(0,12,2),xlab="Time (Years)",legend=F,logrank=T,
     atrisk.title="No. at Risk",ylab="Survival Fraction",xlim=c(0,12),col=c(1),
     marktime=T)
# I think when you ahve such low number of events its good to see the uncertainty
# in the survival probability - you will easily con yourself otherwise

summary(factor(rsi.all$Event_OS[rsi.all$Site=="endometrial"]))
summary(factor(rsi.all$Event_OS[rsi.all$Site=="glioma"]))
summary(factor(rsi.all$Event_OS[rsi.all$Site=="lung"]))
summary(factor(rsi.all$Event_OS[rsi.all$Site=="pancreas"]))

# remove the TNBC cohort and explore variance in total dose
rsi.all<-rsi.all[rsi.all$Site!="breast_TN",]
ggplot()+geom_jitter(data=rsi.all,aes(Site,TD,col=Site),height=0,width=0.2,alpha=0.3)+
  theme_bw(base_size=16)+xlab("Tumour Site")+ylab("Total RT Dose (Gy)")
# we can see a range of total doses, with the highest dose overall
max(rsi.all$TD)# 70Gy
# If we assume 2Gy per fraction and dose given 5 days a week, we would need
# 14 weeks of observation to have then collected enough data
# Are there any events between 0 and 14 weeks or patients who
# become right-censored
sum(rsi.all$Time_OS<=98/365)
# no events so need to worry about time-dependent bias. If we had seen 
# any we should have probably considered a landamrk analysis 

# now we perform the assessments of interest
m0<-(coxph(Surv(Time_OS,Event_OS)~RSI,data=rsi.all[rsi.all$Site=="endometrial",]))
summary(m0)#0.55 (0.05)
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD,data=rsi.all[rsi.all$Site=="endometrial",]))
summary(m1)# 
anova(m0,m1)#0.671

# glioma where they also adjusted for MGMT expression - we shall first look
# without it - its important to note that for some reason this was the only
# tumour type where any adjustment was done
m0<-(coxph(Surv(Time_OS,Event_OS)~RSI,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m0)#0.51 (0.03)
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m1)# 
anova(m0,m1)#0.664

m2<-(coxph(Surv(Time_OS,Event_OS)~MGMT_Expression,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m2)#0.59 (0.03)
m3<-(coxph(Surv(Time_OS,Event_OS)~MGMT_Expression+RSI,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m3)#0.60 (0.03)
anova(m3,m2)#0.050
m4<-(coxph(Surv(Time_OS,Event_OS)~MGMT_Expression+RSI+TD,data=rsi.all[rsi.all$Site=="glioma",]))
summary(m4)# 
anova(m3,m4)#0.0.814
# After adjusting for MGMT we do see RSI explain a little of the variance but
# still not elucidating a dose-response - you could argue there isn't enough
# dose variation to see one but other results and prior knowledge tells us
# thats unlikely to be the case. This in itself tells a story there are far more
# important factors than small variations in dose that drive the response variable

m0<-(coxph(Surv(Time_OS,Event_OS)~RSI,data=rsi.all[rsi.all$Site=="lung",]))
summary(m0)#0.55 (0.04)
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD,data=rsi.all[rsi.all$Site=="lung",]))
summary(m1)# 
anova(m0,m1)#0.538

m0<-(coxph(Surv(Time_OS,Event_OS)~RSI,data=rsi.all[rsi.all$Site=="pancreas",]))
summary(m0)#0.55 (0.06)
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD,data=rsi.all[rsi.all$Site=="pancreas",]))
summary(m1)# 
anova(m0,m1)#0.585


m0<-(coxph(Surv(Time_OS,Event_OS)~RSI+strata(Site),data=rsi.all))
summary(m0)#0.52 (0.02)
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI*Site+strata(Site),data=rsi.all))
anova(m0,m1)#0.455
# We can relax with the interaction assumption
m1<-(coxph(Surv(Time_OS,Event_OS)~RSI+TD+strata(Site),data=rsi.all))
summary(m1)# 
anova(m0,m1)#0.584

# interestingly...
m0<-(coxph(Surv(Time_OS,Event_OS)~gard+strata(Site),data=rsi.all))
summary(m0)#0.52 (0.02)
# GARD is useless too.  Its likely that the breast data-set may have been
# driving everything, melanoma may have too, but thats not been provided.
m1<-(coxph(Surv(Time_OS,Event_OS)~gard*Site+strata(Site),data=rsi.all))
summary(m1)# 
anova(m0,m1)#0.823



