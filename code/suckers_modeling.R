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
vbStarts<- list(Linf=500, K=0.1, t0=-3)
vbTypical<-Len~Linf*(1-exp(-K*(age-t0)))
data_forfitting<- femaleLH[,c("age", "Len")]
data_forfitting$Linf<- vbStarts$Linf
#data_forfitting$K<- vbStarts$K
fitTypical<-nls(vbTypical,data=data_forfitting, start=list(K=vbStarts$K, t0=-3))

plot(femaleLH$age, femaleLH$Len, xlim=c(1,20), ylim=c(100, 600))
exes<- 1:20
whys<- vbStarts$Linf*(1-exp(-coef(fitTypical)[1]*(exes-coef(fitTypical)[2])))
lines(exes, whys)

# sd about mean: Pierce et al. say that they use max(L_obs)-Linf, but I can't make it make sense.
grow_sd<- 25 #abs(max(femaleLH$Len)-500)
growth_params<- list(Linf=500, K=coef(fitTypical)[1], t0=coef(fitTypical)[2],
                     grow_sd = grow_sd)

###########################################################################
### Survival model
###########################################################################
# point estimates from literature
surv_points<- data.frame(len=c(10, 20, 30, 112, 200, 400),
                         surv = 1-c(0.997, 0.997, 0.997, 0.97, 0.25, 0.25))

# fit a logistic curve:
surv_model<- nls(surv~Smax/(1+exp(-k*(len-x0))), data=surv_points,
                 algorithm = "port",
                 start = list(Smax = 0.75, k = 0.1, x0 = 130),
                 lower = rep(0,3), upper = c(1, Inf, 600))
fitted_surv<- function(x){coef(surv_model)[1] / (1+exp(-coef(surv_model)[2]*(x-coef(surv_model)[3])))}

# plot:
plot(surv_points$len, surv_points$surv)
exes<- 1:600
whys<- fitted_surv(exes)
lines(exes, whys)

###########################################################################
### EGG PRODUCTION model
###########################################################################
# linear model in log-log space (Carl's results)
egg_model<- lm(formula = log(fecundity) ~ log(Len), data=femaleLH)

egg_logslope = egg_model$coefficients[2] # 2.776441
egg_logintercept =  egg_model$coefficients[1] # -7.988805

# plot it:
plot(femaleLH$Len, femaleLH$fecundity, xlim=c(0,600),
     xlab='Female length (mm)', ylab='Annual egg production')
exes<- 1:600
whys<- exp(egg_logslope*log(exes) + egg_logintercept)
lines(exes, whys)

###########################################################################
### Maturity ogive
###########################################################################
# expectation: ~5% of age 2 spawn, over 50% at age 3, 90% from age 4 onwards
# point estimates from literature
matur_points<- data.frame(len=c(70, 112, 120, 140, 170, 400),
                          p_spawn = c(0, 0.01, 0.05, 0.5, 0.9, 0.9))

# fit a logistic curve:
matur_model<- nls(p_spawn~Pmax/(1+exp(-k*(len-x0))), data=matur_points,
                  algorithm = "port",
                  start = list(Pmax = 0.9, k = 0.1, x0 = 140),
                  lower = rep(0,3), upper = c(1, Inf, 600))
fitted_matur<- function(x){coef(matur_model)[1] / (1+exp(-coef(matur_model)[2]*(x-coef(matur_model)[3])))}

# plot:
plot(matur_points$len, matur_points$p_spawn)
exes<- 1:600
whys<- fitted_matur(exes)
lines(exes, whys)


###########################################################################
### Mostly following Pierce et al. 2023, model building:
###########################################################################

# assign parameters:
m_par <- list(
  ## Growth parameters
  grow_rate = growth_params$K, # growth rate
  Linf  =   growth_params$Linf, # maximum length in mm
  grow_sd   =   growth_params$grow_sd,  # growth sd
  ## Survival parameters
  surv_max = coef(surv_model)[1], # maximum survival value
  surv_k = coef(surv_model)[2], # rate of increase of survival
  surv_midsize = coef(surv_model)[3], # size at which survival is halfway between upper and lower limit
  ## Size of age-1 individuals:
  recruit_mean = 112, # mean size of age-1 individuals
  recruit_sd = growth_params$grow_sd, # same as grow_sd
  ## PLACEHOLDER:
  egg_viable = 0.1,
  ## Estimated from fecundity data
  egg_logslope = egg_model$coefficients[2], # 2.776441
  egg_logintercept =  egg_model$coefficients[1], # -7.988805
  ## Spawning Probability
  pb_max = coef(matur_model)[1], # maximum probability of spawning
  pb_k = coef(matur_model)[2], # rate of increase of spawning probability with size
  pb_midsize = coef(matur_model)[3], # size at which 50% of individuals spawn
  ## YOY survival probability:
  s0= 0.03 # PLACEHOLDER
)

##########################
## Section 2: Model Set-up
##########################

## Growth function
# given you are size z now returns the pdf of size z1 next time
# computed from von Bertanaffy equation z(t) = L_inf(1-e^K(t-t0))
# to find z(t+1) = L_inf*(1-e^(-K)) + e^(-K)*z(t)

g_z1z <- function(z1, z, m_par) {
  mu <- m_par$Linf * (1 - exp(- m_par$grow_rate)) +
    exp(-m_par$grow_rate) * z           # mean size next year
  sig <- m_par$grow_sd                       # sd about mean
  p_den_grow <- dnorm(z1, mean = mu, sd = sig)    # pdf that you are size z1
  # given you were size z
  return(p_den_grow)
}

## Adult Survival function, 3-parameter logistic
s_z <- function(z, m_par) {
  m_par$surv_max / (1 + exp( -m_par$surv_k * (z - m_par$surv_midsize)))
}

## Reproduction, log-linear
eggs_z <- function(z, m_par) { # Eggs produced (note: data are in thousands)
  eggz<- exp(m_par$egg_logslope*log(z) + m_par$egg_logintercept)
  return(eggz)
}

## Probability of spawning, 3-parameter logistic
pb_z<- function(z, m_par){
  m_par$pb_max / (1 + exp( -m_par$pb_k * (z - m_par$pb_midsize)))
}

## Recruit size pdf
c_1z1 <- function(z1, m_par) {
  mu <- m_par$recruit_mean
  sig <- m_par$recruit_sd
  p_den_recruit <- dnorm(z1, mean = mu, sd = sig)
  return(p_den_recruit)
}

#####################################################
## Section 3 - Build IPM kernels F and P
#####################################################

## Fecundity Kernel
f_z1z <- function(z1, z, m_par) {
  age1_dist <- pb_z(z, m_par) * eggs_z(z, m_par) *
    m_par$egg_viable * m_par$s0
  #returns fecundity kernel (as a matrix). Recruits= F.dot(n*delta_z)
  return(outer(c_1z1(z1, m_par), age1_dist))
}

## Growth and Survival Kernel
p_z1z <- function(z1, z, m_par) {
  N<- length(z1)
  delta_z<- z1[2]-z1[1]
  g_matrix <- matrix(0, N, N)
  for (x in 1:N) {
    g_matrix[, x] <- g_z1z(z, rep(z[x], times = N), m_par)
    g_matrix[, x] <- g_matrix[, x] / (sum(g_matrix[, x]) * delta_z)
  }
  return(g_matrix %*% diag(s_z(z, m_par)))
}

## Build the deterministic kernels ##
m <- 300 # number of meshpoints: bins for the integration
L <- 0.00   # lower size limit in mm
U <- 600.00    # upper size limit in mm - must be larger than Linf
h <- (U-L)/m # integration bin width
meshpts <-  L + (1:m)*h - h/2

Pmat<- h*p_z1z(meshpts, meshpts, m_par)
Fmat<- h*(f_z1z(meshpts, meshpts, m_par))
# P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
# F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
Kmat<- Pmat+Fmat

## Plot the kernels to check it looks okay
matrix.image(Pmat, x=meshpts, y=meshpts, main='Growth+Survival')
matrix.image(Fmat, x=meshpts, y=meshpts, main='Reproduction')
matrix.image(Kmat, x=meshpts, y=meshpts, main='Projection Kernel')

## Calculate a few metrics to see how the model is behaving:
eigz<- eigen(Kmat)
# population growth rate
lambda<- max(Re(eigz$values))

# calculate average lifespan:
lifespan(Pmat) # note that this is dependent on starting size.
# Small individuals (including YOY) live only one year on average.
# Larger individuals have an expectation of more years of life.

# plot the population size distribution at stable growth
popvec<- eigz$vectors[,1]
popvec<- Re(popvec)/sum(Re(popvec))
par(mfrow=c(1,1))
barplot(popvec, names.arg = meshpts)
