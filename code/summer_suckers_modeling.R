###########################################################################
### Modeling life history of white suckers and summer suckers: summer sucker version
### Christina Hernandez
### July 2023
###########################################################################
# Load necessary packages:
library(readxl)
library(here)
source(here("code", "MatrixImage.R"))

# Read in data:
observations<- read_xlsx(here("data", "sucker_model_data.xlsx"))
observations_ss <- observations[observations$length < 350,]
lifehistory<- read_xlsx(here("data", "life_history_model.xlsx"))
lifehistory_ss <- lifehistory[lifehistory$genotype_cands < 0.8,]
# get rid of blank rows in both tables:
I<- which(is.na(observations_ss$year))
observations_ss<- observations_ss[-I,]
I<- which(is.na(lifehistory_ss$year))
lifehistory_ss<- lifehistory_ss[-I,]
# convert fecundity to numeric:
lifehistory_ss$fecundity<- as.numeric(lifehistory_ss$fecundity)
# fix the mislabeled male who produced eggs:
I<- which(lifehistory_ss$sex=="M" & lifehistory_ss$fecundity>0)
lifehistory_ss$sex[I]<- "F"
# pull out females only:
females<- observations_ss[observations_ss$sex=="F",]
femaleLH<- lifehistory_ss[lifehistory_ss$sex=="F",]

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
vbStarts<- list(Linf=350, K=0.1, t0=-3)
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
grow_sd<- 21 #abs(max(femaleLH$Len)-350)
growth_params<- list(Linf=350, K=coef(fitTypical)[1], t0=coef(fitTypical)[2],
                     grow_sd = grow_sd)

###########################################################################
### Survival model
###########################################################################
size_breaks<- seq(100, 350, 10)
yr1_counts<- hist(EKL_2020$length, breaks = size_breaks)

# Project these individuals forward using the growth function:
g_z1z <- function(z1, z, m_par) {
  mu <- m_par$Linf * (1 - exp(- m_par$K)) +
    exp(-m_par$K) * z           # mean size next year
  sig <- m_par$grow_sd                       # sd about mean
  p_den_grow <- dnorm(z1, mean = mu, sd = sig)    # pdf that you are size z1
  # given you were size z
  return(p_den_grow)
}

# now use that pdf to project the sizes forward:
N<- length(size_breaks)-1
L <- min(size_breaks)   # lower size limit in mm
U <- max(size_breaks)    # upper size limit in mm - must larger than Linf
h <- (U-L)/N # integration bin width
meshpts <-  L + (1:N)*h - h/2

g_matrix <- matrix(0, N, N)
for (x in 1:N) {
  g_matrix[, x] <- g_z1z(meshpts, rep(meshpts[x], times = N), growth_params)
  g_matrix[, x] <- g_matrix[, x] / (sum(g_matrix[, x]) * h)
}

yr2_predicted<- h * g_matrix %*% yr1_counts$counts

# add reproduction:
Fmat<- h*(f_z1z(meshpts, meshpts, m_par))
yr2_predicted<- yr2_predicted + Fmat %*% yr1_counts$counts


# Compare the predicted size distribution to the observed distribution:
yr2_counts<- hist(EKL_2021$length, breaks = size_breaks)

par(mfrow=c(1,1))
plot(meshpts, yr2_predicted, type='l')
lines(meshpts, yr2_counts$counts, col='red')
lines(meshpts, yr1_counts$counts, col='blue')

# okay this doesn't work, need to discuss with Steve.

###########################################################################
### EGG PRODUCTION model
###########################################################################
# linear model in log-log space (Carl's results)
egg_model<- lm(formula = log(fecundity) ~ log(Len), data=femaleLH)

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
  surv_min  =  0.03, # previous placeholder was .002 from gizzard shad Bodola (1955)
  # 0.03 from Bagley et al 2018 based on White Sucker data from McPhee 2007
  # Then et al (2015): surv_max computed from
  # 1-natural mortality = 1 - 8.872*K^.73 L^-.33
  surv_max = 1 - 8.872*growth_params$K^.73*growth_params$Linf^(-.33),
  # inflection point: will be temp dependent
  # computed for La Grange Reach
  surv_alpha = 200, # PLACEHOLDER 80.0136 - taken from Pierce et al. 2023
  # surv_alpha roughly represents the point at which predation decreases substantially.
  # This could be from gape limitation.
  surv_beta = -5, # 139.9312 slope PLACEHOLDER - taken from Pierce.
  ## How steep is the curve surrounding the inflection point?
  ## -5 from grass carp paper - Erickson et al
  ## Size of age-1 individuals:
  recruit_mean = 112, # mean size of age-1 individuals
  recruit_sd = growth_params$grow_sd, # same as grow_sd
  ## PLACEHOLDER:
  egg_viable = 0.0035, # from Bagley et al 2018 survival probability from egg to YOY
  ## Estimated from fecundity data
  egg_logslope = egg_model$coefficients[2], # 2.776441
  egg_logintercept =  egg_model$coefficients[1], # -7.988805
  ## Spawning Probability
  prob_spawn = 0.90, # PLACEHOLDER
  ## YOY survival probability:
  s0= 0.01 # PLACEHOLDER
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

## Adult Survival function, 4-parameter logistic
s_z <- function(z, m_par) {
  m_par$surv_min + (m_par$surv_max - m_par$surv_min) /
    (1 + exp(m_par$surv_beta * (log(z) - log(m_par$surv_alpha))))
}

## Reproduction, log-linear
eggs_z <- function(z, m_par) { # Eggs produced (note: data are in thousands)
  eggz<- exp(m_par$egg_logslope*log(z) + m_par$egg_logintercept)
  return(eggz)
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
  age1_dist <- m_par$prob_spawn * eggs_z(z, m_par) *
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

# Plot the kernels to check it looks okay
matrix.image(Pmat, x=meshpts, y=meshpts, main='Growth+Survival')
matrix.image(Fmat, x=meshpts, y=meshpts, main='Reproduction')
matrix.image(Kmat, x=meshpts, y=meshpts, main='Projection Kernel')

# population growth rate (currently, it shows the population as just barely growing
max(Re(eigen(Kmat, only.values=TRUE)$values))

