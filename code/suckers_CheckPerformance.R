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
plot(femaleLH$Len, femaleLH$fecundity, xlim=c(0,600),
     xlab='Female length (mm)', ylab='Annual egg production')
exes<- 1:600
# Note that if you haven't run the suckers_modeling script first, then this line won't run
whys<- exp(m_par$egg_logslope*log(exes) + m_par$egg_logintercept)
# hard-coded version (note that if anything changes about the data, this line would be incorrect)
#whys<- exp(2.776441*log(exes) - 7.988805)
lines(exes, whys)

# plot the stable size distribution without any YOY:
ws_dist_noYOY<- stable_size_dist_noYOY(Pmat, Fmat)
plot(meshpts, ws_dist_noYOY, type='l', xlab='Length (mm)', ylab='Probability density')

# plot the stable size distribution, assuming no individuals under 80 mm caught:
ws_dist_noCatch80mm<- stable_size_dist_sizeThreshold(Pmat, Fmat, meshpts, threshold=80)
plot(meshpts, ws_dist_noCatch100mm, type='l', xlab='Length(mm)', ylab='Probability density')

# calculate the age at which 95% of individuals are dead:
age_most_individuals_dead(Pmat, Fmat, proportion=0.95)
# calculate the age at which 99% of individuals are dead:
age_most_individuals_dead(Pmat, Fmat, proportion=0.99)


