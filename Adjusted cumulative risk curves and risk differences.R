# The objective of this R script is to present two cumulative risk plots in a 
# regression discontinuity design adjusting for a calendar time trend.
# We then estimate the risk difference at 10 years of follow up.

# Authors: Andreas Rieckmann, Jessica Bengtsson, Emilie Courtin, Vahe Nafilyan
# Date: May 20, 2022

# # To assess consistency across multiple runs, start (comments, #, can be removed)
results = NA
for (i in 1:200) {
  print(i)
# # To assess consistency across multiple runs, end (comments, #, can be removed)  

#### Load packages ####
library(survival)
library(RColorBrewer)

#### Data simulation ####
# We simulate a data set with individuals who are followed for 10 years
n = 100000
id = 1:n
t = rnorm(n,10, 4)
# They are sampled from the years 2000 to 2010
inclusion_time = runif(n, min = 2000, max = 2010)
inclusion_year <- floor(inclusion_time)
nat_exp_exposed = ifelse (inclusion_time > 2005, 1 , 0) # One is exposed to the natural experiment if the calendar year is more than 2005
time_difference_to_natural_experiment = inclusion_time - 2005 # Re-coding the inclusion_time to the years since the natural experiment 
# The outcome, d, depends on a time effect and the natural experiment
p_all = 0.05 # all have an absolute baseline risk
p_time = 0.005 # the relative risk increased with 0.01 per calendar year. 
p_nat_exp = 0 # the effect the natural experiment is none. Note that an increased risk of the natural experiment will not be directly reflected in the final risk difference because of censoring at 10 years. 
d = rbinom(n, size = 1, prob = p_all  + p_time * time_difference_to_natural_experiment +  p_nat_exp * nat_exp_exposed)
# Censoring at 12 years of follow up
t = ifelse (t < 0 , 0 , t) # To ensure nobody leaves the study before baseline (t=0)
d = ifelse (t >= 10 , 0 , d) # To censor individuals at t = 10, such that events occurring after t = 10 are not counted.
t = ifelse (t >= 10 , 10 , t)
# The below plot shows the cumulative risk according to the inclusion_year. Blue/green are early years, red/black are late years
par(mfrow=c(1,1))
plot(survfit(Surv(t,d) ~ as.factor(inclusion_year) + cluster(id)), fun = "event",xlim = c(0,11),col=rev(brewer.pal(10,"Spectral")))

#### Plotting crude and adjusted cumulative risk ####
par(mfrow=c(1,2))

# Crude cumulative risk
fit_crude <- survfit(coxph(Surv(t,d) ~ strata(nat_exp_exposed) + cluster(id)))
plot(fit_crude, ylim=c(0,0.06), xlim = c(0,15), main = "Unadjusted", fun = "event", yaxs = 'i', xaxs = 'i', col = 1:2, lwd = 2, conf.int = T, xlab = "Years of follow up", ylab = "Risk")

# Adjusted for the calendar time effect cumulative risk
time_difference_to_natural_experiment_2 <- time_difference_to_natural_experiment * time_difference_to_natural_experiment
fit <- coxph(Surv(t,d) ~ as.numeric(time_difference_to_natural_experiment) + as.numeric(time_difference_to_natural_experiment_2) + strata(nat_exp_exposed) + cluster(id))
fit_adj <- survfit(fit,newdata = data.frame(time_difference_to_natural_experiment=0, time_difference_to_natural_experiment_2 = 0))
plot(fit_adj, ylim=c(0,0.06), xlim = c(0,15), main = "Adjusted for a calendar time effect", fun = "event", yaxs = 'i', xaxs = 'i', col = 1:2, lwd = 2, conf.int = T, xlab = "Years of follow up", ylab = "Risk")
abline(h=seq(0,0.08,0.001),col=adjustcolor("grey",0.2))

#### Risk difference after 10 years of follow up adjusting for a linear effect of calendar time #####
# We extract the survival functions for each intervention group
fit_adj_0 <- data.frame(time = fit_adj$time[1:fit_adj$strata[1]],
                     surv = fit_adj$surv[1:fit_adj$strata[1]],
                     lower = fit_adj$lower[1:fit_adj$strata[1]],
                     upper = fit_adj$upper[1:fit_adj$strata[1]])
fit_adj_1 <- data.frame(time = fit_adj$time[c(fit_adj$strata[1]+1):c(sum(fit_adj$strata))],
                     surv = fit_adj$surv[c(fit_adj$strata[1]+1):c(sum(fit_adj$strata))],
                     lower = fit_adj$lower[c(fit_adj$strata[1]+1):c(sum(fit_adj$strata))],
                     upper = fit_adj$upper[c(fit_adj$strata[1]+1):c(sum(fit_adj$strata))])

# Risk difference at 10 years of follow up (we find the se on a log scale and simulate the difference)
p_0<-fit_adj_0$surv[which.max(fit_adj_0$time[fit_adj_0$time<10])]
l_0<-fit_adj_0$lower[which.max(fit_adj_0$time[fit_adj_0$time<10])]
u_0<-fit_adj_0$upper[which.max(fit_adj_0$time[fit_adj_0$time<10])]
se_0 <- mean(c(log(p_0)-log(l_0),log(u_0)-log(p_0))) / 1.96
risk_0 <- quantile(1-exp(rnorm(10000,log(p_0),se_0)),c(0.025,0.5,0.975))
risk_0 <- format(round(risk_0,3),nsmall=3)
text(10,1-p_0,paste0(risk_0[2]," (",risk_0[1],";",risk_0[3],")"),pos=4, cex = 0.5)

p_1<-fit_adj_1$surv[which.max(fit_adj_1$time[fit_adj_1$time<10])]
l_1<-fit_adj_1$lower[which.max(fit_adj_1$time[fit_adj_1$time<10])]
u_1<-fit_adj_1$upper[which.max(fit_adj_1$time[fit_adj_1$time<10])]
se_1 <- mean(c(log(p_1)-log(l_1),log(u_1)-log(p_1))) / 1.96
risk_1 <- quantile(1-exp(rnorm(10000,log(p_1),se_1)),c(0.025,0.5,0.975))
risk_1 <- format(round(risk_1,3),nsmall=3)
text(10,1-p_1,paste0(risk_1[2]," (",risk_1[1],";",risk_1[3],")"),pos=4, cex = 0.5, col = "red")

(ests<-quantile((1-exp(rnorm(100000,log(p_1),se_1)))-(1-exp(rnorm(100000,log(p_0),se_0))),c(0.025,0.5,0.975))) # RR inc 95%CI
rd <- ests[2] # For assessing consistency
ests <- format(round(ests,3),nsmall=3)
text(1,0.05,paste0("Adjusted risk difference at 10 y\n",ests[2],"\n(",ests[1],";",ests[3],")"), pos = 4)
# This is the risk difference after 10 years of follow up among those exposed to the natural experiment compared to those who were not exposed to the natural experiment
# Risk ratios can also be estimated from the risks in each group.

# # To assess consistency across multiple runs, start (comments, #, can be removed)
results[i] <- rd
print(quantile(results,c(0.025,0.5,0.975))) # For assessing consistency
}
par(mfrow=c(1,1))
hist(results) # The distribution of the results from repeating the simulation multiple times.
# # To assess consistency across multiple runs, end (comments, #, can be removed)
