###########################################################################
### Modeling life history of white suckers and summer suckers
### Christina Hernandez
### July 2023
###########################################################################
# Load necessary packages:
library(readxl)
library(here)
library(exactLTRE)
library(tidyverse) # sorry Chrissy! I'll only use it for plotting, I swear...
library(RColorBrewer)
source(here("code", "MatrixImage.R"))
source('code/suckers_functions.R')

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
femaleLH_ss <- femaleLH[femaleLH$age < 10,]

###########################################################################
### Plot annual size distributions
###########################################################################
EKL_2020<- females[females$year==2020,]
EKL_2021<- females[females$year==2021,]

hist(EKL_2020$length)
hist(EKL_2021$length)

###########################################################################
### Growth model
###########################################################################

## Fit a von Bertalanffy model using non-linear least squares and a fixed Linf

# Take Linf from the literature: 500 mm, and this is reasonable for the dataset because the largest individual was 465
vbStarts<- list(Linf=450, K=0.1, t0=-3)
vbTypical<-Len~Linf*(1-exp(-K*(age-t0)))
data_forfitting<- femaleLH[,c("age", "Len")]
data_forfitting$Linf<- vbStarts$Linf
#data_forfitting$K<- vbStarts$K
fitTypical<-nls(vbTypical,data=data_forfitting, start=list(K=vbStarts$K, t0=-3))

plot(femaleLH$age, femaleLH$Len, xlim=c(1,20), ylim=c(100, 600), xaxt="n", main = "Growth Rates")
axis(1, at = c(1:20))
exes<- 1:20
whys <- vbStarts$Linf*(1-exp(-coef(fitTypical)[1]*(exes-coef(fitTypical)[2])))
lines(exes, whys, lty = 2)
# sd about mean: Pierce et al. say that they use max(L_obs)-Linf, but I can't make it make sense.
# grow_sd<- 25 #abs(max(femaleLH$Len)-500)
# growth_params<- list(Linf=500, K=coef(fitTypical)[1], t0=coef(fitTypical)[2],
#                         grow_sd = grow_sd)

###########################################################################
### Survival model
###########################################################################
# point estimates from literature
surv_points<- data.frame(len=c(130, 140, 150, 180, 250, 300),
                            surv = 1-c(0.99, 0.99, 0.99, 0.92, 0.25, 0.25))

# fit a logistic curve:
surv_model<- nls(surv~Smax/(1+exp(-k*(len-x0))), data=surv_points,
                    algorithm = "port",
                    start = list(Smax = 0.80, k = 0.05, x0 = 150),
                    lower = rep(0,3), upper = c(1, Inf, 600))
fitted_surv<- function(x){coef(surv_model)[1] / (1+exp(-coef(surv_model)[2]*(x-coef(surv_model)[3])))}

# Survival model b (4-parameter)
surv_min <-  0.003
surv_mid <- 0.62
surv_max <- 0.75 # controlled by sliding table
surv_alpha <- 112
surv_alpha2 <- 220
surv_beta <- -25
surv_beta2 <- -25
seven_fitted_surv<- function(z) {
  surv_min + ((surv_mid - surv_min) /
    (1 + exp(surv_beta * (log(z) - log(surv_alpha)))))  +
    ((surv_max - surv_mid) /
                  (1 + exp(surv_beta2 * (log(z) - log(surv_alpha2)))))
}

len <- 1:600
survl <- seven_fitted_surv(len)
surv_dat <- data.frame(len, survl)
surv_dat %>% ggplot(aes(x = len, y = survl)) +
  geom_line() +
  annotate("line",x = exes, y = ss_four_whys, lwd = 1, lty = 2) +
  labs(x = "Length (mm)", y = "Survival") +
  cowplot::theme_cowplot()

###########################################################################
### EGG PRODUCTION model
###########################################################################
# linear model in log-log space (Carl's results)
egg_model<- lm(formula = log(fecundity) ~ log(Len), data=femaleLH)

egg_logslope = egg_model$coefficients[2] # 3.1082
egg_logintercept =  egg_model$coefficients[1] # -9.7183

###########################################################################
### Maturity ogive
###########################################################################
# expectation: ~5% of age 2 spawn, over 50% at age 3, 90% from age 4 onwards
# point estimates from literature
matur_points<- data.frame(len=c(150-10, 192-10, 200-10, 220-10, 260, 490),
                             p_spawn = c(0, 0.01, 0.05, 0.5, 0.65, 0.65))

# fit a logistic curve:
matur_model<- nls(p_spawn~Pmax/(1+exp(-k*(len-x0))), data=matur_points,
                     algorithm = "port",
                     start = list(Pmax = 0.65, k = 0.1, x0 = 230),
                     lower = rep(0,3), upper = c(1, Inf, 600))
fitted_matur<- function(x){coef(matur_model)[1] / (1+exp(-coef(matur_model)[2]*(x-coef(matur_model)[3])))}

## try manually setting k and Pmax
pb_max_man <- 0.65
pb_k_man <- 0.1634

###########################################################################
### Mostly following Pierce et al. 2023, model building:
###########################################################################
## make table with values of surv_max and pb_midsize
# just biologically accurate values
# y <- rep(seq(100,300,by = 5), each = 126) # pb_midsize
# z <- rep(seq(250,500, by = 2), time = 41) # Linf

# all values but sub-sampled
#w <- rep(seq(350,500, by =6.25), each = 25) # size cut off to fit data
x <- c(rep(seq(0.58888,0.75,by = 0.005035), each = 31), rep(0.75,each=248)) # max survival
y <- rep(seq(121.6,297.6,by = 4.4), each = 31) # pb_midsize
v <- rep(seq(-1.091682,-0.3280479, by = 0.02545447), time = 41) # t0
w <- rep(seq(0.1485557,0.098245, by = -0.00167702), time = 41) # k rate
z <- rep(seq(317,482, by = 5.5), time = 41) # Linf *fit to 270, not 280
u <- c(rep(0.022, each = 341),rep(seq(0.022,0.015, by = - 0.0003043478), each = 31), rep(0.015, each = 186)) #egg viable
slide_params <- data.frame(t0 = v, k = w, pb_midsize = y, Linf = z, max_surv = x, egg_viable = u) %>% filter(pb_midsize < Linf)
# bistabe Linf-Lm relationship: pb_midsize <= (1.125*Linf - 158) + 5 & pb_midsize >= (1.125*Linf - 158) - 5
# slide_params <- data.frame(pb_midsize = y, Linf = z) %>% filter(Linf <= (2.5*pb_midsize - 100) + 20 & Linf >= (2.5*pb_midsize - 100) - 20)
## Run model for corresponding values of surv_max and mat_alpha
lambda_tab <- data.frame(pb_midsize = slide_params$pb_midsize, Linf = slide_params$Linf, max_surv = slide_params$max_surv,egg_viable = slide_params$egg_viable, lambda = NA, K = NA, t0 = NA, age_999 = NA)
surv_tab <- data.frame(len = NA, surv = NA, pb_midsize = NA)

for(i in 1:length(slide_params$pb_midsize)){
  #ifelse(slide_params$Linf[i] == 258.3, data_forfitting <- femaleLH, data_forfitting <- femaleLH_ss)
# smoothly scale vb equation for each Linf
  # sd about mean: Pierce et al. say that they use max(L_obs)-Linf, but I can't make it make sense.
  grow_sd<- 25 #abs(max(femaleLH$Len)-500)
  growth_params<- list(Linf=slide_params$Linf[i], K=slide_params$k[i], t0=slide_params$t0[i], grow_sd = grow_sd)
# assign parameters:
m_par <- list(
  ## Growth parameters
  grow_rate = slide_params$k[i], # growth rate
  Linf = slide_params$Linf[i], # scale maximum length with size at maturity in mm
  grow_sd   = growth_params$grow_sd,  # growth sd #growth_params$Linf,
  ## Survival 7-parameter, double logit
  surv_min = surv_min,
  surv_mid = ifelse(slide_params$max_surv[i]<surv_mid,slide_params$max_surv[i],surv_mid),
  surv_max = slide_params$max_surv[i],
  surv_alpha = surv_alpha,
  surv_alpha2 = surv_alpha2,
  surv_beta = surv_beta,
  surv_beta2 = surv_beta2,
  ## Size of age-1 individuals:
  recruit_mean = 112, # mean size of age-1 individuals
  recruit_sd = growth_params$grow_sd, # same as grow_sd
  ## PLACEHOLDER:
  egg_viable = slide_params$egg_viable[i],
  ## Estimated from fecundity data
  egg_logslope = egg_logslope, #egg_model$coefficients[2], # 3.1082
  egg_logintercept = egg_logintercept, #egg_model$coefficients[1], # -9.7183
  ## Spawning Probability
  pb_max = pb_max_man, #coef(matur_model)[1], # maximum probability of spawning
  pb_k = pb_k_man, #coef(matur_model)[2], # rate of increase of spawning probability with size
  pb_midsize = slide_params$pb_midsize[i], # size at which 50% of individuals spawn
  ## YOY survival probability:
  s0= 0.1 # Begley et al 2017
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

## Adult Survival function, double 4-param logistic
s_z <- function(z, m_par) {
  m_par$surv_min + ((m_par$surv_mid - m_par$surv_min) /
    (1 + exp(m_par$surv_beta * (log(z) - log(m_par$surv_alpha)))))+
    ((m_par$surv_max - m_par$surv_mid) /
       (1 + exp(m_par$surv_beta2 * (log(z) - log(m_par$surv_alpha2)))))
}

## Adult Survival function b, 4-parameter logistic
# s_z <- function(z, m_par) {
#   m_par$surv_min + (m_par$surv_max - m_par$surv_min) /
#     (1 + exp(m_par$surv_beta * (log(z) - log(m_par$surv_alpha))))
# }

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
m_Kmat<- Kmat^0.3

## Plot the kernels to check it looks okay
# matrix.image(Pmat, x=meshpts, y=meshpts, main='Growth+Survival')
# matrix.image(Fmat, x=meshpts, y=meshpts, main='Reproduction')
# matrix.image(m_Kmat, x=meshpts, y=meshpts, main='Projection Kernel^0.3')

## Calculate a few metrics to see how the model is behaving:
eigz<- eigen(Kmat)
# population growth rate
lambda<- max(Re(eigz$values))
lambda_tab$lambda[i] <- lambda
lambda_tab$K[i] <- growth_params$K
lambda_tab$t0[i] <- growth_params$t0
lambda_tab$age_999[i] <- age_most_individuals_dead(Pmat, Fmat, proportion=0.999)

## Make data.frames that have survival curve data
surv_at_size<- colSums(Pmat)
# assign(paste0("surv_out","_",lambda_tab$pb_midsize[i]), data.frame(len = meshpts, surv = surv_at_size))
tmp <- data.frame(len = meshpts, surv = surv_at_size, pb_midsize = lambda_tab$pb_midsize[i])
surv_tab <- rbind(surv_tab, tmp)
}
# calculate average lifespan:
#lifespan(Pmat) # note that this is dependent on starting size.
# Small individuals (including YOY) live only one year on average.
# Larger individuals have an expectation of more years of life.

# plot the population size distribution at stable growth
# popvec<- eigz$vectors[,1]
# popvec<- Re(popvec)/sum(Re(popvec))
# par(mfrow=c(1,1))
# barplot(popvec, names.arg = meshpts)

par(mfrow = c(1,1))
#plot(lambda_tab$lambda~lambda_tab$pb_midsize)
#plot(lambda_tab$surv_max~lambda_tab$pb_midsize)
#lines(lambda_tab$pb_midsize, lambda_tab$Linf*0.001)
#surv_tab %>% group_by(pb_midsize) %>% ggplot(aes(x = len, y = surv, group = pb_midsize)) +
#  geom_line(aes(color = pb_midsize), cex = 1) +
#  scale_color_continuous(type = "viridis") +
#  cowplot::theme_cowplot()
lambda_tab %>% ggplot(aes(x = pb_midsize, y = lambda)) +
  geom_point(aes(fill = surv_max), shape = 21, color = "black", size = 3) +
  scale_fill_continuous(type = "viridis")

## estimate Linf from length at maturity
#Linf for females
# equation from Froese and Binohlan 2000.
# original equation Lm = exp(0.9469*log(Linf)-0.1162)
#lambda_tab <- lambda_tab %>% mutate(Lm = exp(1.92*log(Linf)-6.04)) # tweaked from Froesse and Binghlan 2000
len2age <- function(Lm,Linf,k,t0){
  ifelse(Lm != Linf,age <- ((log(1-Lm/Linf))/-k)+t0, age <- NA)
  return(round(age,0))
}

# Suckers lose ~0.02561 of their length to reproduce each year
age_Lm_Linf <- data.frame(Linf = seq(317,482, by = 5.5), Lm = seq(140,290,by = 5), K = seq(0.1485557,0.098245, by = -0.00167702), t0 = seq(-1.091682,-0.3280479, by = 0.02545447))
age_Lm_Linf <- age_Lm_Linf %>% mutate(age = round(((log(1-Lm/Linf))/-K)+t0)) %>%
  mutate(Linf_est = 460/(1+0.0398)^(10-age))


lambda_tab <- lambda_tab %>% mutate(age_estm = ifelse(pb_midsize < Linf,round(((log(1-pb_midsize/Linf))/-K)+t0,0), NA)) %>%
  mutate(Linf_est = ifelse(age_estm <= 11, 460/(1+0.07)^(8-age_estm), NA))
Lm_Linf <- lambda_tab %>% filter(Linf <= Linf_est+5 & Linf > Linf_est-2)
# percent length lost per year of spawning = 0.0398

# plot fitness heatmap with Lm-Linf relationship
my_colors <- c("midnightblue","lightblue","mediumpurple","salmon", "orangered4")
# my_values <- c(0,0.32,0.36,0.40,1)
my_values <- c(0,0.61,0.65,0.69,1)
my_points <- data.frame(x = c(170,270), y = c(350,460))
#my_values <- c(0,0.38,0.40,0.42,1)

lambda_tab %>% ggplot(aes(x = pb_midsize, y = Linf)) +
  geom_tile(aes(fill = lambda)) +
  scale_fill_gradientn(colors = my_colors, values = my_values) +
  #geom_line(aes(x = Lm,y = Linf), lwd = 1, color = "red") +
  annotate("line",x = Lm_Linf$pb_midsize, y = Lm_Linf$Linf_est, lwd = 1) +
  annotate("point",x = my_points$x, y = my_points$y) +
  xlab("Lm") +
  #geom_line(aes(x = Lm, y = Linf), lty = 2) +
  cowplot::theme_cowplot()

lambda_tab_bacc <- lambda_tab %>% filter(Linf <= Linf_est+12 & Linf > Linf_est-2)# + 5 & Linf >= (1.25*pb_midsize + 100) - 5)

my_values <- c(0,0.66,0.70,0.74,1)
lambda_tab_bacc %>% ggplot(aes(x = pb_midsize, y = Linf)) +
  geom_tile(aes(fill = lambda)) +
  scale_fill_gradientn(colors = my_colors, values = my_values) +
  xlab("Lm") +
  #geom_line(aes(x = Lm2,y = Linf)) +
  cowplot::theme_cowplot()

lambda_tab_bacc %>% ggplot(aes(x = pb_midsize, y = lambda)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  geom_smooth() +
  xlab("Lm")
  cowplot::theme_cowplot()

lambda_tab_bacc %>% ggplot(aes(x = Linf, y = lambda)) +
  geom_point() +
  geom_hline(yintercept = 1)

# bacc_surv_up <- c(seq(0.71,0.80,length.out = 11), seq(0.80, 0.80, length.out = 20), seq(0.80,0.80,length.out = 10)) # max survival
# bacc_surv_dn <- c(seq(0.68,0.75,length.out = 11), seq(0.75, 0.77, length.out = 20), seq(0.77,0.77,length.out = 10))
# bacc_mat <- seq(100,300, by = 5)
# bacc_tab <- data.frame(bacc_surv_up, bacc_surv_dn, bacc_mat)
# surv_max <- NA
# pb_midsize <- NA
# Linf <- NA
# lambda <- NA
# lambda_tab_bacc <- data.frame(surv_max,pb_midsize,Linf,lambda)
#
# for(i in 1:length(lambda_tab$pb_midsize)){
#   for(k in 1:length(bacc_tab$bacc_mat)){
#     if((lambda_tab$pb_midsize[i] == bacc_tab$bacc_mat[k] & lambda_tab$surv_max[i] <= bacc_tab$bacc_surv_up[k] & lambda_tab$surv_max[i] >= bacc_tab$bacc_surv_dn[k]) == T){
#       lambda_tab_bacc <- rbind(lambda_tab_bacc,lambda_tab[i,])
#       #print(lambda_tab[i,])
#     }
#   }
# }
#
# my_values <- c(0,0.41,0.46,0.51,1)
#
# lambda_tab_bacc %>% ggplot(aes(x = pb_midsize, y = Linf)) +
#   ylim(c(0.7,0.85)) +
#   geom_tile(aes(fill = lambda)) +
#   scale_fill_gradientn(colors = my_colors, values = my_values) +
#   cowplot::theme_cowplot()
