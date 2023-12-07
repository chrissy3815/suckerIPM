###########################################################################
### Modeling life history of white suckers and summer suckers: summer sucker version
### Christina Hernandez
### July 2023
###########################################################################
## Change log
# survival model changed to 4 param with shallow slope
# egg viability changed from 0.1 to 0.0035 -> this led to bimodal size distribution

# Load necessary packages:
library(readxl)
library(here)
library(exactLTRE)
source(here("code", "MatrixImage.R"))

# Read in data:
observations<- read_xlsx(here("data", "sucker_model_data.xlsx"))
lifehistory<- read_xlsx(here("data", "life_history_model.xlsx"))
# get rid of blank rows in observations table:
I<- which(is.na(observations$year))
observations<- observations[-I,]
# remove outlier fecundity for EKL_107
lifehistory$fecundity[3] <- NA
# convert fecundity to numeric:
lifehistory$fecundity<- as.numeric(lifehistory$fecundity)
# fix the mislabeled male who produced eggs:
I<- which(lifehistory$sex=="M" & lifehistory$fecundity>0)
lifehistory$sex[I]<- "F"
# pull out females only:
females<- observations[observations$sex=="F",]
femaleLH<- lifehistory[lifehistory$sex=="F",]

###########################################################################
### Plot annual size distributions
###########################################################################
EKL_2020<- females[females$year==2020,]
EKL_2021<- females[females$year==2021,]

par(mfrow=c(2,1))
hist(EKL_2020$length, breaks = 20)
hist(EKL_2021$length, breaks = 30)

###########################################################################
### Growth model
###########################################################################

# linear model:
lineargrowth<- lm(Len~age, femaleLH)
summary(lineargrowth)

# scatter plot:
par(mfrow=c(1,1))
plot(femaleLH$age, femaleLH$Len)

## Fit a von Bertalanffy model using non-linear least squares and a fixed Linf

# Take Linf from the literature: 500 mm, and this is reasonable for the dataset because the largest individual was 465
ss_vbStarts<- list(Linf=500, K=0.1, t0=-3)
ss_vbTypical<-Len~Linf*(1-exp(-K*(age-t0)))
ss_data_forfitting<- femaleLH[,c("age", "Len")]
ss_data_forfitting$Linf<- ss_vbStarts$Linf
#data_forfitting$K<- vbStarts$K
ss_fitTypical<-nls(ss_vbTypical,data=ss_data_forfitting, start=list(K=ss_vbStarts$K, t0=-3))

plot(femaleLH$age, femaleLH$Len, xlim=c(1,20), ylim=c(100, 600), xaxt="n")
axis(1, at = c(1:20))
exes<- 1:20
whys<- ss_vbStarts$Linf*(1-exp(-coef(ss_fitTypical)[1]*(exes-coef(ss_fitTypical)[2])))
lines(exes, whys)

# sd about mean: Pierce et al. say that they use max(L_obs)-Linf, but I can't make it make sense.
ss_grow_sd<- 25 #abs(max(femaleLH$Len)-500)
ss_growth_params<- list(Linf=500, K=coef(ss_fitTypical)[1], t0=coef(ss_fitTypical)[2],
                     grow_sd = ss_grow_sd)

###########################################################################
### Survival model
###########################################################################
# point estimates from literature
ss_surv_points<- data.frame(len=c(10, 20, 30, 80, 220, 330),
                         surv = 1-c(0.997, 0.997, 0.997, 0.97, 0.4, 0.4))

# fit a logistic curve:
# Survival model a (3-parameter)
ss_surv_model<- nls(surv~Smax/(1+exp(-k*(len-x0))), data=ss_surv_points,
                 algorithm = "port",
                 start = list(Smax = 0.6, k = 0.1, x0 = 130),
                 lower = rep(0,3), upper = c(1, Inf, 600))
ss_fitted_surv<- function(x){coef(ss_surv_model)[1] / (1+exp(-coef(ss_surv_model)[2]*(x-coef(ss_surv_model)[3])))}
ss_varied_surv<- function(x){0.75 / (1+exp(-.03*(x-150)))}

# plot:
plot(ss_surv_points$len, ss_surv_points$surv)
exes<- 1:600
whys<- ss_fitted_surv(exes)
ss_varied_whys <- ss_varied_surv(exes)
lines(exes, whys)



###########################################################################
### EGG PRODUCTION model
###########################################################################
# linear model in log-log space (Carl's results)
ss_egg_model<- lm(formula = log(fecundity) ~ log(Len), data=femaleLH)

ss_egg_logslope = ss_egg_model$coefficients[2] # 3.1082
ss_test_logslope = 3.4
ss_egg_logintercept =  ss_egg_model$coefficients[1] # -9.7183
ss_test_logintercept = -11.35

# plot it:
plot(femaleLH$Len, femaleLH$fecundity, xlim=c(0,600),
     xlab='Female length (mm)', ylab='Annual egg production')
exes<- 1:600
ss_test_whys<- exp(ss_test_logslope*log(exes) + ss_test_logintercept)
whys<- exp(ss_egg_logslope*log(exes) + ss_egg_logintercept)
#lines(exes, whys)
lines(exes, ss_test_whys)



###########################################################################
### Maturity ogive
###########################################################################
# expectation: ~5% of age 2 spawn, over 50% at age 3, 90% from age 4 onwards
# point estimates from Elk Lake data
ss_matur_points<- data.frame(len=c(90, 132, 140, 160, 190, 420),
                          p_spawn = c(0, 0.01, 0.05, 0.5, 0.9, 0.9))

# fit a logistic curve:
ss_matur_model<- nls(p_spawn~Pmax/(1+exp(-k*(len-x0))), data=ss_matur_points,
                  algorithm = "port",
                  start = list(Pmax = 0.9, k = 0.1, x0 = 160),
                  lower = rep(0,3), upper = c(1, Inf, 600))
ss_fitted_matur<- function(x){coef(ss_matur_model)[1] / (1+exp(-coef(ss_matur_model)[2]*(x-coef(ss_matur_model)[3])))}

# plot:
plot(ss_matur_points$len, ss_matur_points$p_spawn)
exes<- 1:600
whys<- ss_fitted_matur(exes)
lines(exes, whys, lty = 2)
lines(exes, ws_whys)

legend(250,0.4, c("White Sucker", "Summer Sucker"), lty = c(1,2), cex = 1.7)


###########################################################################
### Mostly following Pierce et al. 2023, model building:
###########################################################################

# assign parameters:
ss_m_par <- list(
  ## Growth parameters
  grow_rate = ss_growth_params$K, # growth rate
  Linf  = ss_growth_params$Linf, # maximum length in mm
  grow_sd   = ss_growth_params$grow_sd,  # growth sd
  ## Survival parameters a
  surv_max = coef(ss_surv_model)[1], # maximum survival value
  surv_k = coef(ss_surv_model)[2], # rate of increase of survival
  surv_midsize = coef(ss_surv_model)[3], # size at which survival is halfway between upper and lower limit
  ## Size of age-1 individuals:
  recruit_mean = 112, # mean size of age-1 individuals
  recruit_sd = ss_growth_params$grow_sd, # same as grow_sd
  ## PLACEHOLDER:
  egg_viable = 0.03,
  ## Estimated from fecundity data
  egg_logslope = ss_test_logslope, #egg_model$coefficients[2], # 2.776441
  egg_logintercept = ss_test_logintercept, #egg_model$coefficients[1], # -7.988805
  ## Spawning Probability
  pb_max = coef(ss_matur_model)[1], # maximum probability of spawning
  pb_k = coef(ss_matur_model)[2], # rate of increase of spawning probability with size
  pb_midsize = coef(ss_matur_model)[3], # size at which 50% of individuals spawn
  ## YOY survival probability:
  s0= 0.1 # PLACEHOLDER
)

##########################
## Section 2: Model Set-up
##########################

## Growth function
# given you are size z now returns the pdf of size z1 next time
# computed from von Bertanaffy equation z(t) = L_inf(1-e^K(t-t0))
# to find z(t+1) = L_inf*(1-e^(-K)) + e^(-K)*z(t)

ss_g_z1z <- function(z1, z, ss_m_par) {
  mu <- ss_m_par$Linf * (1 - exp(- ss_m_par$grow_rate)) +
    exp(-ss_m_par$grow_rate) * z           # mean size next year
  sig <- ss_m_par$grow_sd                       # sd about mean
  p_den_grow <- dnorm(z1, mean = mu, sd = sig)    # pdf that you are size z1
  # given you were size z
  return(p_den_grow)
}

# Adult Survival function a, 3-parameter logistic
ss_s_z <- function(z, ss_m_par) {
  ss_m_par$surv_max / (1 + exp( -ss_m_par$surv_k * (z - ss_m_par$surv_midsize)))
}

## Adult Survival function b, 4-parameter logistic
# s_z <- function(z, ss_m_par) {
#   ss_m_par$surv_min + (ss_m_par$surv_max - ss_m_par$surv_min) /
#     (1 + exp(ss_m_par$surv_beta * (log(z) - log(ss_m_par$surv_alpha))))
# }

## Reproduction, log-linear
ss_eggs_z <- function(z, ss_m_par) { # Eggs produced (note: data are in thousands)
  eggz<- exp(ss_m_par$egg_logslope*log(z) + ss_m_par$egg_logintercept)
  return(eggz)
}

## Probability of spawning, 3-parameter logistic
ss_pb_z<- function(z, ss_m_par){
  ss_m_par$pb_max / (1 + exp(-ss_m_par$pb_k * (z - ss_m_par$pb_midsize)))
}

## Recruit size pdf
ss_c_1z1 <- function(z1, ss_m_par) {
  mu <- ss_m_par$recruit_mean
  sig <- ss_m_par$recruit_sd
  p_den_recruit <- dnorm(z1, mean = mu, sd = sig)
  return(p_den_recruit)
}

#####################################################
## Section 3 - Build IPM kernels F and P
#####################################################

## Fecundity Kernel
ss_f_z1z <- function(z1, z, ss_m_par) {
  age1_dist <- ss_pb_z(z, ss_m_par) * ss_eggs_z(z, ss_m_par) *
    ss_m_par$egg_viable * ss_m_par$s0
  #returns fecundity kernel (as a matrix). Recruits= F.dot(n*delta_z)
  return(outer(ss_c_1z1(z1, ss_m_par), age1_dist))
}

## Growth and Survival Kernel
ss_p_z1z <- function(z1, z, ss_m_par) {
  N<- length(z1)
  delta_z<- z1[2]-z1[1]
  g_matrix <- matrix(0, N, N)
  for (x in 1:N) {
    g_matrix[, x] <- ss_g_z1z(z, rep(z[x], times = N), ss_m_par)
    g_matrix[, x] <- g_matrix[, x] / (sum(g_matrix[, x]) * delta_z)
  }
  return(g_matrix %*% diag(ss_s_z(z, ss_m_par)))
}

## Build the deterministic kernels ##
m <- 300 # number of meshpoints: bins for the integration
L <- 0.00   # lower size limit in mm
U <- 600.00    # upper size limit in mm - must be larger than Linf
h <- (U-L)/m # integration bin width
ss_meshpts <-  L + (1:m)*h - h/2

ss_Pmat<- h*ss_p_z1z(ss_meshpts, ss_meshpts, ss_m_par)
ss_Fmat<- h*(ss_f_z1z(ss_meshpts, ss_meshpts, ss_m_par))
# P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
# F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
ss_Kmat<- ss_Pmat+ss_Fmat
m_ss_Kmat<- ss_Kmat^0.3

## Plot the kernels to check it looks okay
matrix.image(ss_Pmat, x=ss_meshpts, y=ss_meshpts, main='Growth+Survival')
matrix.image(ss_Fmat, x=ss_meshpts, y=ss_meshpts, main='Reproduction')
matrix.image(ss_Kmat, x=ss_meshpts, y=ss_meshpts, main='Projection Kernel')
matrix.image(m_ss_Kmat, x=ss_meshpts, y=ss_meshpts, main='Projection Kernel^0.3')

## Calculate a few metrics to see how the model is behaving:
ss_eigz<- eigen(ss_Kmat)
# population growth rate
ss_lambda<- max(Re(ss_eigz$values))

# calculate average lifespan:
lifespan(ss_Pmat) # note that this is dependent on starting size.
# Small individuals (including YOY) live only one year on average.
# Larger individuals have an expectation of more years of life.

# plot the population size distribution at stable growth
ss_popvec<- ss_eigz$vectors[,1]
ss_popvec<- Re(ss_popvec)/sum(Re(ss_popvec))
par(mfrow=c(1,1))
barplot(ss_popvec, names.arg = ss_meshpts)
ss_lambda
