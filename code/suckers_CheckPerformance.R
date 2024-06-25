###########################################################################
### Checking/plotting behavior of models
### Christina Hernandez
### Nov 2023
###########################################################################

## Notes from Chrissy to Carl: You might want to think about how to structure
#this code smartly. For example, we might want to name things as ws_Pmat and
#ss_Pmat for the various projection kernels. Actually, it might be useful to
#name everything in the white sucker script with "ws_" at the beginning
#(functions, variables, etc.) so that we can have all of the white sucker
#objects and the summer sucker objects active at the same time, to enable
#comparisons.

# source the functions:
source('code/suckers_functions.R')

# plot the egg production (does not differ between phenotypes):
par(mfrow = c(1,1))
par(mar = c(5,5,2,5))
plot(femaleLH$Len, femaleLH$fecundity, xlim=c(0,600),
     xlab='Female length (mm)', ylab='Annual egg production')
exes<- 1:600
# Note that if you haven't run the suckers_modeling script first, then this line won't run
ss_whys<- exp(ss_m_par$egg_logslope*log(exes) + ss_m_par$egg_logintercept)
ws_whys<- exp(ws_m_par$egg_logslope*log(exes) + ws_m_par$egg_logintercept)
# hard-coded version (note that if anything changes about the data, this line would be incorrect)
#whys<- exp(2.776441*log(exes) - 7.988805)
lines(exes, ws_whys, lwd = 3)


# eigenvalues and eigenvectors:
ws_eigz<- eigen(ws_Kmat)
ss_eigz<- eigen(ss_Kmat)

# population growth rate
paste("White Sucker lambda = ", max(Re(ws_eigz$values)))
paste("Summer Sucker lambda = ", max(Re(ss_eigz$values)))


# calculate average lifespan:
lifespan(ws_Pmat)
lifespan(ss_Pmat) # note that this is dependent on starting size.
# Small individuals (including YOY) live only one year on average.
# Larger individuals have an expectation of more years of life.

# plot the stable size distribution
ws_popvec<- ws_eigz$vectors[,1]
ws_popvec<- Re(ws_popvec)/sum(Re(ws_popvec))
ss_popvec<- ss_eigz$vectors[,1]
ss_popvec<- Re(ss_popvec)/sum(Re(ss_popvec))
plot(ws_meshpts, ws_popvec, type='l', xlab='Length (mm)', ylab='Probability density', ylim = c(0,.03), main = "Length distributions")
lines(ss_meshpts, ss_popvec, lty = 2)
legend(350,0.015, c("White Sucker", "Summer Sucker"), lty = c(1,2))

# plot the stable size distribution without any YOY:
ws_dist_noYOY<- stable_size_dist_noYOY(ws_Pmat, ws_Fmat)
ss_dist_noYOY<- stable_size_dist_noYOY(ss_Pmat, ss_Fmat)
plot(ws_meshpts, ws_dist_noYOY, type='l', xlab='Length (mm)', ylab='Probability density', main = "length distributions (noYOY)")
lines(ss_meshpts, ss_dist_noYOY, lty = 2)
legend(350,0.01, c("White Sucker", "Summer Sucker"), lty = c(1,2))

# plot the stable size distribution, assuming no individuals under 80 mm caught:
ws_dist_noCatch80mm<- stable_size_dist_sizeThreshold(ws_Pmat, ws_Fmat, ws_meshpts, threshold=80)
ss_dist_noCatch80mm<- stable_size_dist_sizeThreshold(ss_Pmat, ss_Fmat, ss_meshpts, threshold=80)
plot(ws_meshpts, ws_dist_noCatch80mm, ylim = c(0,0.001), type='l', xlab='Length (mm)', ylab='Probability density', main = "length distributions over 80mm",xlim = c(300,550))
lines(ss_meshpts, ss_dist_noCatch80mm, lty = 2)
legend(220,0.02, c("White Sucker", "Summer Sucker"), lty = c(1,2), cex = 2)

# calculate the age at which 99% of individuals are dead:
paste("Age at 99% White Suckers dead = ",age_most_individuals_dead(ws_Pmat, ws_Fmat, proportion=0.99))
paste("Age at 99% Summer Suckers dead = ",age_most_individuals_dead(ss_Pmat, ss_Fmat, proportion=0.99))
paste("Age at 99.9% White Suckers dead = ",age_most_individuals_dead(ws_Pmat, ws_Fmat, proportion=0.999))
paste("Age at 99.9% Summer Suckers dead = ",age_most_individuals_dead(ss_Pmat, ss_Fmat, proportion=0.999))

# Plot growthrate

plot(femaleLH$age, femaleLH$Len, xlim=c(1,20), ylim=c(100, 600), xaxt="n", xlab = "Age", ylab = "Length (mm)")
axis(1, at = c(1:20))
exes<- 1:20
whys<- ss_vbStarts$Linf*(1-exp(-coef(ss_fitTypical)[1]*(exes-coef(ss_fitTypical)[2])))
ws_whys <- ws_vbStarts$Linf*(1-exp(-coef(ws_fitTypical)[1]*(exes-coef(ws_fitTypical)[2])))
lines(exes, whys, lty = 2)
lines(exes, ws_whys, lty = 1)
legend(1,550, c("WS ecotype", "SS ecotype"), lty = c(1,2), cex = 1.4)

# Plot the survival-at-size curve from the kernel:
ws_surv_at_size<- colSums(ws_Pmat)
ss_surv_at_size<- colSums(ss_Pmat)
plot(ws_meshpts, ws_surv_at_size, type='l', lwd = 3, xlab='Length (mm)', ylab='Probability of survival to t+1', main = "Survival-at-size curves")
lines(ss_meshpts, ss_surv_at_size, lwd = 3, lty = 2)
legend(250,0.4, c("WS ecotype", "SS ecotype"), lty = c(1,2), cex = 1.4)

# Plot the maturity ogive
exes<- 1:600
whys<- ss_fitted_matur(exes)
ws_whys <- ws_fitted_matur(exes)
plot(exes, whys, type = "l", lty = 2, lwd = 3, xlab = "Length (mm)", ylab = "Proportion spawning per year")
lines(exes, ws_whys, lwd = 3)
legend(300,0.4, c("WS ecotype", "SS ecotype"), lty = c(1,2), cex = 1.7)

# Plot the mature length distribution from the kernel:
ws_prob_matur<- ws_pb_z(ws_meshpts, ws_m_par)
ws_length_dist_matur<- ws_popvec*ws_prob_matur
ss_prob_matur<- ss_pb_z(ss_meshpts, ss_m_par)
ss_length_dist_matur<- ss_popvec*ss_prob_matur
plot(ss_meshpts, ss_length_dist_matur, lty = 2, type = 'l', xlab='Length (mm)', ylab='Probability density', ylim = c(0,0.007), xlim = c(0,600), main = "Mature length distribution")
lines(ws_meshpts, ws_length_dist_matur, lty = 1)
legend(270,0.005, c("White Sucker", "Summer Sucker"), lty = c(1,2), cex = 1.4)

# Mature length dist zoom
ws_prob_matur<- ws_pb_z(ws_meshpts, ws_m_par)
ws_length_dist_matur<- ws_popvec*ws_prob_matur
ss_prob_matur<- ss_pb_z(ss_meshpts, ss_m_par)
ss_length_dist_matur<- ss_popvec*ss_prob_matur
plot(ss_meshpts, ss_length_dist_matur, lty = 2, type = 'l', xlab='Length (mm)', ylab='Probability density', ylim = c(0,0.0055), xlim = c(100,500), main = "Mature length distribution")
lines(ws_meshpts, ws_length_dist_matur, lty = 1)
legend(270,0.005, c("White Sucker", "Summer Sucker"), lty = c(1,2), cex = 1.4)
plot(ws_meshpts, ws_length_dist_matur, lty = 2, type = 'l', xlab='Length (mm)', ylab='Probability density', ylim = c(0,0.0018), xlim = c(220,500), main = "Mature length distribution")
