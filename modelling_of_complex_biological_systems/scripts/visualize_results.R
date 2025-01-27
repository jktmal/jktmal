library(tidyverse)
library(rstan)
library(lubridate)
library(stats4)
library(DescTools)
#install.packages("DescTools")

plot(1:3, 1:3)

#load("/home/jasiek/Dropbox/mcbs/project/results/mu_ifr.RData")

# Indiana
pop.IN <- 6785644
IN_fit <- readRDS("/home/jasiek/Dropbox/mcbs/project/results/IN_fit.rds")  
IN_sims <- rstan::extract(IN_fit)
load("/home/jasiek/Dropbox/mcbs/project/data/IN_ctp.RData")

median_nu = apply(IN_sims$nu, 2, median)
lwr.nu = apply(IN_sims$nu, 2, sort)[j,]
upr.nu = apply(IN_sims$nu, 2, sort)[k,]
median_ifr = median(IN_sims$ifr)
ifr = MeanCI(IN_sims$ifr)

a = IN_sims$nu[,375]
j=0.025*length(a)
k=0.975*length(a)

cc[,375]/(ctp_data$cum_cc[375]/pop.IN)
cc[,375]

#plot 1
plot(seq(1:375), median_nu, 'l', ylim = c(0, 20000), xaxt = "n", xlab = "", ylab="New infections", col='black') 
polygon(c(rev(seq(1:375)), seq(1:375)), c(rev(upr.nu), lwr.nu), col = "grey80", border = NA, density=200) + lines(seq(1:375), median_nu) + lines(seq(1:375), ctp_data$d/median_ifr, col = 'red')
axis(1, at=c(yday("2020-03-01")-57, yday("2020-06-01")-57, yday("2020-09-01")-57, yday("2020-12-01")-57, yday("2021-03-01")+364-57), labels = c("3/20", "6/20", "9/20", "12/20", "3/21"))


#yday("2020-02-27")

#plot 5
cc = rbind(cumsum(median_nu), cumsum(lwr.nu), cumsum(upr.nu)) / pop.IN

plot(seq(1:375), cc[1,]*100, 'l', ylim = c(0, 25), xlab = "", xaxt = "n", ylab = "Cumulative incidence (% of pop.)")
polygon(c(rev(seq(1:375)), seq(1:375)), c(rev(cc[3,]*100), cc[2,]*100), col = 'grey80', border = NA, density=200) + lines(seq(1:375), cc[1,]*100)
axis(1, at=c(yday("2020-03-01")-57, yday("2020-06-01")-57, yday("2020-09-01")-57, yday("2020-12-01")-57, yday("2021-03-01")+364-57), labels = c("3/20", "6/20", "9/20", "12/20", "3/21"))


#plot 4
#median_rt = apply(IN_sims$beta, 2, median)*median(IN_sims$gamma_inv)
ci_rt = apply(IN_sims$beta, 2, median)*median(IN_sims$gamma_inv)
lwr.rt = apply(IN_sims$beta, 2, sort)[j,]*median(IN_sims$gamma_inv)
upr.rt = apply(IN_sims$beta, 2, sort)[k,]*median(IN_sims$gamma_inv)

plot(seq(1:375), ci_rt, 'l', xlab="", xaxt="n", ylab="Reproductive number r(t)")
polygon(c(rev(seq(1:375)), seq(1:375)), c(rev(upr.rt), lwr.rt), col = 'grey80', border = NA, density=200) + lines(seq(1:375), ci_rt)
axis(1, at=c(yday("2020-03-01")-57, yday("2020-06-01")-57, yday("2020-09-01")-57, yday("2020-12-01")-57, yday("2021-03-01")+364-57), labels = c("3/20", "6/20", "9/20", "12/20", "3/21"))

# yday("2020-06-01") - 96

#plot 3 - undercount factor
ci_ct = ctp_data$cum_cc
ci_R <- apply(IN_sims$y[, , 2] + IN_sims$y[, , 3], 2, median)[2:376]
lwr.RI <- apply(IN_sims$y[, , 2] + IN_sims$y[, , 3], 2, sort)[j, 2:376] / ci_ct
upr.RI <- apply(IN_sims$y[, , 2] + IN_sims$y[, , 3], 2, sort)[k, 2:376] / ci_ct

uc <- (ci_R / ci_ct)[96:375]
plot(seq(96:375), uc, 'l', xlab="", xaxt="n", ylab="Undercount factor") 
polygon(c(rev(seq(96:375)), seq(96:375)), c(rev(upr.RI[96:375]), lwr.RI[96:375]), col = 'grey80', border = NA, density=200) + lines(seq(96:375), uc)
axis(1, at=c(yday("2020-03-01")-57, yday("2020-06-01")-57, yday("2020-09-01")-57, yday("2020-12-01")-57), labels = c("6/20", "9/20", "12/20", "3/21"))

#plot extra
plot(seq(1:375), median_nu, 'l', ylim = c(0, 20000), xaxt = "n", xlab = "", ylab="New infections", col='blue') + lines(seq(1:375), ctp_data$cc, col='black') 
axis(1, at=c(yday("2020-03-01")-57, yday("2020-06-01")-57, yday("2020-09-01")-57, yday("2020-12-01")-57, yday("2021-03-01")+364-57), labels = c("3/20", "6/20", "9/20", "12/20", "3/21"))


#plot 5 - unused
median.ct = median(IN_sims$phi) * ci_R * sqrt(ctp_data$cum_tests/pop.IN)
plot(cc[1,], median.ct/pop.IN, 'l')#, xlim=c(0.0, 0.03))#, ylim=c(0,0.08))

#plot 6 -unused
median.phi_t = median(IN_sims$phi)*sqrt(ctp_data$cum_tests/pop.IN)
median.phi_t1 = c( median.phi_t[1], median(IN_sims$phi)*sqrt(ctp_data$cum_tests[1:374]/pop.IN))
median.RI_t1 = c(ci_R[2], ci_R[2:375])

mean.ct = median_nu*median.phi_t + (median.phi_t - median.phi_t1)*median.RI_t1
plot(mean.ct, ctp_data$cc, 'p', xlim=c(0, 7000))#, xlim=c(0.0, 0.03))#, ylim=c(0,0.08))
scatter.smooth(mean.ct, ctp_data$cc, xlim=c(0,7000))

# Ohio

pop.OH <- 11790587
OH_fit <- readRDS("/home/jasiek/Dropbox/mcbs/project/results/OH_fit.rds")  
OH_sims <- rstan::extract(OH_fit)
load("/home/jasiek/Dropbox/mcbs/project/data/OH_ctp.RData")


median_nu = apply(OH_sims$nu, 2, median)
lwr.nu = apply(OH_sims$nu, 2, sort)[j,]
upr.nu = apply(OH_sims$nu, 2, sort)[k,]
median_ifr = median(OH_sims$ifr)
ifr = MeanCI(OH_sims$ifr)

a = OH_sims$nu[,368]
j=0.025*length(a)
k=0.975*length(a)

#plot 1
plot(seq(1:368), median_nu, 'l', ylim = c(0, 20000), xaxt = "n", xlab = "", ylab="New infections", col='blue') + lines(seq(1:368), ctp_data$cc, col='black') 
axis(1, at=c(yday("2020-03-01")-57, yday("2020-06-01")-57, yday("2020-09-01")-57, yday("2020-12-01")-57, yday("2021-03-01")+364-57), labels = c("3/20", "6/20", "9/20", "12/20", "3/21"))

#plot 2
plot(seq(1:368), median_nu, 'l', ylim = c(0, 20000), xaxt = "n", xlab = "", ylab="New infections", col='black') 
polygon(c(rev(seq(1:368)), seq(1:368)), c(rev(upr.nu), lwr.nu), col = "grey80", border = NA, density=200) + lines(seq(1:368), median_nu) + lines(seq(1:368), ctp_data$d/median_ifr, col = 'red')
axis(1, at=c(yday("2020-03-01")-57, yday("2020-06-01")-57, yday("2020-09-01")-57, yday("2020-12-01")-57, yday("2021-03-01")+364-57), labels = c("3/20", "6/20", "9/20", "12/20", "3/21"))

#plot 3 - undercount factor
ci_ct = ctp_data$cum_cc
ci_R <- apply(OH_sims$y[, , 2] + OH_sims$y[, , 3], 2, median)[2:369]
lwr.RI <- apply(OH_sims$y[, , 2] + OH_sims$y[, , 3], 2, sort)[j, 2:369] / ci_ct
upr.RI <- apply(OH_sims$y[, , 2] + OH_sims$y[, , 3], 2, sort)[k, 2:369] / ci_ct

uc <- (ci_R / ci_ct)[96:369]
plot(seq(96:369), uc, 'l', xlab="", xaxt="n", ylab="Undercount factor") 
polygon(c(rev(seq(96:369)), seq(96:369)), c(rev(upr.RI[96:369]), lwr.RI[96:369]), col = 'grey80', border = NA, density=200) + lines(seq(96:369), uc)
axis(1, at=c(yday("2020-03-01")-57, yday("2020-06-01")-57, yday("2020-09-01")-57, yday("2020-12-01")-57), labels = c("6/20", "9/20", "12/20", "3/21"))

#plot 4
ci_rt = apply(OH_sims$beta, 2, median)*median(OH_sims$gamma_inv)
lwr.rt = apply(OH_sims$beta, 2, sort)[j,]*median(OH_sims$gamma_inv)
upr.rt = apply(OH_sims$beta, 2, sort)[k,]*median(OH_sims$gamma_inv)

plot(seq(1:368), ci_rt, 'l', xlab="", xaxt="n", ylab="Reproductive number r(t)")
polygon(c(rev(seq(1:368)), seq(1:368)), c(rev(upr.rt), lwr.rt), col = 'grey80', border = NA, density=200) + lines(seq(1:368), ci_rt)
axis(1, at=c(yday("2020-03-01")-57, yday("2020-06-01")-57, yday("2020-09-01")-57, yday("2020-12-01")-57, yday("2021-03-01")+364-57), labels = c("3/20", "6/20", "9/20", "12/20", "3/21"))

#plot 5
cc = rbind(cumsum(median_nu), cumsum(lwr.nu), cumsum(upr.nu)) / pop.OH

plot(seq(1:368), cc[1,]*100, 'l', ylim = c(0, 25), xlab = "", xaxt = "n", ylab = "Cumulative incidence (% of pop.)")
polygon(c(rev(seq(1:368)), seq(1:368)), c(rev(cc[3,]*100), cc[2,]*100), col = 'grey80', border = NA, density=200) + lines(seq(1:368), cc[1,]*100)
axis(1, at=c(yday("2020-03-01")-57, yday("2020-06-01")-57, yday("2020-09-01")-57, yday("2020-12-01")-57, yday("2021-03-01")+364-57), labels = c("3/20", "6/20", "9/20", "12/20", "3/21"))

cc[,368]/(ctp_data$cum_cc[368]/pop.OH)