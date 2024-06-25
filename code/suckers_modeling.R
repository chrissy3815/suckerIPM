###########################################################################
### Modeling life history of white suckers and summer suckers
### Christina Hernandez
### July 2023
###########################################################################
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
hist(EKL_2020$length)
hist(EKL_2021$length)

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
ws_vbStarts<- list(Linf=460, K=0.1, t0=-3)
ws_vbTypical<-Len~Linf*(1-exp(-K*(age-t0)))
ws_data_forfitting<- femaleLH[,c("age", "Len")]
ws_data_forfitting$Linf<- ws_vbStarts$Linf
#data_forfitting$K<- vbStarts$K
ws_fitTypical<-nls(ws_vbTypical,data=ws_data_forfitting, start=list(K=ws_vbStarts$K, t0=-3))
#ws_fitTypical_ssdata<-nls(ws_vbTypical,data=ss_data_forfitting, start=list(K=ws_vbStarts$K, t0=-3))

plot(femaleLH$age, femaleLH$Len, xlim=c(1,20), ylim=c(100, 600))
axis(1, at = c(1:20))
exes<- 1:20
whys<- ws_vbStarts$Linf*(1-exp(-coef(ws_fitTypical)[1]*(exes-coef(ws_fitTypical)[2])))
#ws_ssdata_whys <- ws_vbStarts$Linf*(1-exp(-coef(ws_fitTypical_ssdata)[1]*(exes-coef(ws_fitTypical_ssdata)[2])))
lines(exes, whys, lwd = 3)
#lines(exes, ws_ssdata_whys, lwd = 3, lty = 2)

# sd about mean: Pierce et al. say that they use max(L_obs)-Linf, but I can't make it make sense.
ws_grow_sd<- 25 #abs(max(femaleLH$Len)-500)
ws_growth_params<- list(Linf=ws_vbStarts$Linf, K=coef(ws_fitTypical)[1], t0=coef(ws_fitTypical)[2],
                     grow_sd = ws_grow_sd)

###########################################################################
### Survival model
###########################################################################
# point estimates from literature
ws_surv_points<- data.frame(len=c(130, 140, 150, 180, 250, 300),
                         surv = 1-c(0.99, 0.99, 0.99, 0.92, 0.25, 0.25))

# fit a logistic curve:
ws_surv_model<- nls(surv~Smax/(1+exp(-k*(len-x0))), data=ws_surv_points,
                 algorithm = "port",
                 start = list(Smax = 0.80, k = 0.05, x0 = 150),
                 lower = rep(0,3), upper = c(1, Inf, 600))
ws_fitted_surv<- function(x){coef(ws_surv_model)[1] / (1+exp(-coef(ws_surv_model)[2]*(x-coef(ws_surv_model)[3])))}

# Survival model b (4-parameter)
# surv_min <-  0.003
# surv_max <- 0.75
# surv_alpha <- 112
# surv_beta <- -17
# four_fitted_surv<- function(z){
#   surv_min + (surv_max - surv_min) /
#     (1 + exp(surv_beta * (log(z) - log(surv_alpha))))
# }
ws_surv_min <-  0.003
ws_surv_mid <- 0.62
ws_surv_max <- 0.75 # controlled by sliding table
ws_surv_alpha <- 112
ws_surv_alpha2 <- 220
ws_surv_beta <- -25
ws_surv_beta2 <- -25
ws_seven_fitted_surv<- function(z) {
  ws_surv_min + ((ws_surv_mid - ws_surv_min) /
                (1 + exp(ws_surv_beta * (log(z) - log(ws_surv_alpha)))))  +
    ((ws_surv_max - ws_surv_mid) /
       (1 + exp(ws_surv_beta2 * (log(z) - log(ws_surv_alpha2)))))
}
# plot:
plot(ws_surv_points$len, ws_surv_points$surv, xlim = c(0,500), ylim = c(0,1), main = "3-param vs 4-param survival curve")
exes<- 1:600
ws_whys<- ws_fitted_surv(exes)
four_whys<- four_fitted_surv(exes)
#lines(exes, ws_whys)
lines(exes, four_whys)
lines(exes, ws_whys, lty = 2)

###########################################################################
### EGG PRODUCTION model
###########################################################################
# linear model in log-log space (Carl's results)
ws_egg_model<- lm(formula = log(fecundity) ~ log(Len), data=femaleLH)

ws_egg_logslope = ws_egg_model$coefficients[2] # 3.1082
ws_test_logslope = 3.4
ws_egg_logintercept =  ws_egg_model$coefficients[1] # -9.7183
ws_test_logintercept = -11.35

# plot it:
plot(femaleLH$Len, femaleLH$fecundity, xlim=c(0,600),
     xlab='Female length (mm)', ylab='Annual egg production')
exes<- 1:600
ws_whys<- exp(ws_egg_logslope*log(exes) + ws_egg_logintercept)
ws_test_whys<- exp(ws_test_logslope*log(exes) + ws_test_logintercept)
lines(exes, ws_whys)
lines(exes, ws_test_whys, lty = 2)

###########################################################################
### Maturity ogive
###########################################################################
# expectation: ~5% of age 2 spawn, over 50% at age 3, 65% from age 4 onwards
# point estimates from literature
ws_matur_points<- data.frame(len=c(200, 225, 250, 280, 320, 490),
                          p_spawn = c(0, 0.01, 0.05, 0.5, 0.65, 0.65))

# fit a logistic curve:
ws_matur_model<- nls(p_spawn~Pmax/(1+exp(-k*(len-x0))), data=ws_matur_points,
                  algorithm = "port",
                  start = list(Pmax = 0.65, k = 0.1, x0 = 230),
                  lower = rep(0,3), upper = c(1, Inf, 600))
ws_fitted_matur<- function(x){coef(ws_matur_model)[1] / (1+exp(-coef(ws_matur_model)[2]*(x-coef(ws_matur_model)[3])))}
# plot:
plot(ws_matur_points$len, ws_matur_points$p_spawn)
exes<- 1:600
ws_whys<- ws_fitted_matur(exes)
lines(exes, ws_whys)

###########################################################################
### Mostly following Pierce et al. 2023, model building:
###########################################################################

# assign parameters:
ws_m_par <- list(
  ## Growth parameters
  grow_rate = ws_growth_params$K, # growth rate
  Linf  = ws_growth_params$Linf, # maximum length in mm
  grow_sd   = ws_growth_params$grow_sd,  # growth sd
  ## Survival parameters a
  #surv_max = coef(ws_surv_model)[1], # maximum survival value
  #surv_k = coef(ws_surv_model)[2], # rate of increase of survival
  #surv_midsize = coef(ws_surv_model)[3], # size at which survival is halfway between upper and lower limit
  ## Survival parameters b
  surv_min =  ws_surv_min,
  surv_mid = ws_surv_mid,
  surv_max = ws_surv_max,
  surv_alpha = ws_surv_alpha,
  surv_alpha2 = ws_surv_alpha2,
  surv_beta = ws_surv_beta,
  surv_beta2 = ws_surv_beta2,
  ## Size of age-1 individuals:
  recruit_mean = 112, # mean size of age-1 individuals
  recruit_sd = ws_growth_params$grow_sd, # same as grow_sd
  ## PLACEHOLDER:
  egg_viable = 0.015,
  ## Estimated from fecundity data
  egg_logslope = ws_test_logslope, #ws_egg_model$coefficients[2], # 3.1082
  egg_logintercept = ws_test_logintercept, #ws_egg_model$coefficients[1], # -9.7183
  ## Spawning Probability
  pb_max = coef(ws_matur_model)[1], # maximum probability of spawning
  pb_k = coef(ws_matur_model)[2], # rate of increase of spawning probability with size
  pb_midsize = coef(ws_matur_model)[3], # size at which 50% of individuals spawn
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

ws_g_z1z <- function(z1, z, ws_m_par) {
  mu <- ws_m_par$Linf * (1 - exp(- ws_m_par$grow_rate)) +
    exp(-ws_m_par$grow_rate) * z           # mean size next year
  sig <- ws_m_par$grow_sd                       # sd about mean
  p_den_grow <- dnorm(z1, mean = mu, sd = sig)    # pdf that you are size z1
  # given you were size z
  return(p_den_grow)
}

## Adult Survival function, 3-parameter logistic
# ws_s_z <- function(z, ws_m_par) {
#   ws_m_par$surv_min + (ws_m_par$surv_max - ws_m_par$surv_min) /
#     (1 + exp(ws_m_par$surv_beta * (log(z) - log(ws_m_par$surv_alpha))))
# }

ws_s_z <- function(z, ws_m_par) {
  ws_m_par$surv_min + ((ws_m_par$surv_mid - ws_m_par$surv_min) /
                      (1 + exp(ws_m_par$surv_beta * (log(z) - log(ws_m_par$surv_alpha)))))+
    ((ws_m_par$surv_max - ws_m_par$surv_mid) /
       (1 + exp(ws_m_par$surv_beta2 * (log(z) - log(ws_m_par$surv_alpha2)))))
}

## Reproduction, log-linear
ws_eggs_z <- function(z, ws_m_par) { # Eggs produced (note: data are in thousands)
  eggz<- exp(ws_m_par$egg_logslope*log(z) + ws_m_par$egg_logintercept)
  return(eggz)
}

## Probability of spawning, 3-parameter logistic
ws_pb_z<- function(z, ws_m_par){
  ws_m_par$pb_max / (1 + exp( -ws_m_par$pb_k * (z - ws_m_par$pb_midsize)))
}

## Recruit size pdf
ws_c_1z1 <- function(z1, ws_m_par) {
  mu <- ws_m_par$recruit_mean
  sig <- ws_m_par$recruit_sd
  p_den_recruit <- dnorm(z1, mean = mu, sd = sig)
  return(p_den_recruit)
}

#####################################################
## Section 3 - Build IPM kernels F and P
#####################################################

## Fecundity Kernel
ws_f_z1z <- function(z1, z, ws_m_par) {
  age1_dist <- ws_pb_z(z, ws_m_par) * ws_eggs_z(z, ws_m_par) *
    ws_m_par$egg_viable * ws_m_par$s0
  #returns fecundity kernel (as a matrix). Recruits= F.dot(n*delta_z)
  return(outer(ws_c_1z1(z1, ws_m_par), age1_dist))
}

## Growth and Survival Kernel
ws_p_z1z <- function(z1, z, ws_m_par) {
  N<- length(z1)
  delta_z<- z1[2]-z1[1]
  g_matrix <- matrix(0, N, N)
  for (x in 1:N) {
    g_matrix[, x] <- ws_g_z1z(z, rep(z[x], times = N), ws_m_par)
    g_matrix[, x] <- g_matrix[, x] / (sum(g_matrix[, x]) * delta_z)
  }
  return(g_matrix %*% diag(ws_s_z(z, ws_m_par)))
}

## Build the deterministic kernels ##
m <- 300 # number of meshpoints: bins for the integration
L <- 0.00   # lower size limit in mm
U <- 600.00    # upper size limit in mm - must be larger than Linf
h <- (U-L)/m # integration bin width
ws_meshpts <-  L + (1:m)*h - h/2

ws_Pmat<- h*ws_p_z1z(ws_meshpts, ws_meshpts, ws_m_par)
ws_Fmat<- h*(ws_f_z1z(ws_meshpts, ws_meshpts, ws_m_par))
# P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
# F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
ws_Kmat<- ws_Pmat+ws_Fmat
m_ws_Kmat<- ws_Kmat^0.3

## Plot the kernels to check it looks okay
# matrix.image(ws_Pmat, x=ws_meshpts, y=ws_meshpts, main='Growth+Survival')
# matrix.image(ws_Fmat, x=ws_meshpts, y=ws_meshpts, main='Reproduction')
# matrix.image(m_ws_Kmat, x=ws_meshpts, y=ws_meshpts, main='Projection Kernel^0.3')

## Calculate a few metrics to see how the model is behaving:
ws_eigz<- eigen(ws_Kmat)
# population growth rate
ws_lambda<- max(Re(ws_eigz$values))
ws_lambda

# calculate average lifespan:
#lifespan(ws_Pmat) # note that this is dependent on starting size.
# Small individuals (including YOY) live only one year on average.
# Larger individuals have an expectation of more years of life.

# plot the population size distribution at stable growth
ws_popvec<- ws_eigz$vectors[,1]
ws_popvec<- Re(ws_popvec)/sum(Re(ws_popvec))
par(mfrow=c(1,1))
barplot(ws_popvec, names.arg = ws_meshpts)

