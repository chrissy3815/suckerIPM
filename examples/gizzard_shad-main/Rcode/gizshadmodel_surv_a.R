# Use a least squares method to find the surv_alpha parameter that 
# minimizes the total square-distance between the (observed) pre-carp 
# LTRM length distribution and (predicted) model equilibrium, 
# say n(z,100). 
## Gizzard Shad (Fall 2021)
# Remember to set the working directory:
# Session -> Set Working Directory -> To Source Location
library(tidyverse)  #used to make pretty graphs
library(lubridate) # used to change date format of data
library(reshape2)

source("gizshadmodel.R")

#######################
## Section 4: Simulation Fun Time!
#######################

N <- 50 # number of size classes
l_shad <- 0.00   # lower size limit in mm
u_shad <- 500.0    # upper size limit in mm - we want this to be
# larger than L-infty
delta_z <- (u_shad - l_shad) / N
zmesh <-  l_shad + ((1:N) - 1 / 2) * (u_shad - l_shad) / N
tf <- 100 # number of years

# Initial length distribution
n <- matrix(0, length(zmesh), tf)
n0_total <- 995 # mean density of a pool in the LTRM dataset
n[, 1] <- dnorm(zmesh, mean = 0.5 * m_par$grow_max, sd = 30)
n[, 1] <- (n[, 1] / sum(n[, 1])) * n0_total / delta_z
# Note: sum(n[,1])*delta_z = n0_total

# Dynamical System
for (i in 1:(tf - 1)) {
  k_iter <- (p_z1z(zmesh, zmesh, m_par) + f_z1z(zmesh, zmesh, n[, i],
                                                m_par)) * delta_z
  n[, i + 1] <- k_iter %*% n[, i]
}

# Total Population Size
n_total <- rep(0, tf)
for (i in 1:tf) {
  n_total[i] <- sum(n[, i]) * delta_z
}

###########################
## Section 5: LTRM DATA
###########################

ltrm_gzsd <- read_csv("ltrm_fish_data.csv")
# Remove length 0 and NA
ltrm_gzsd <- ltrm_gzsd[!is.na(ltrm_gzsd$length) & (ltrm_gzsd$length > 0) &
                         !is.na(ltrm_gzsd$fdate), ]
# Convert date into new format
# Then pull year and add it as a new column
ltrm_gzsd$fdate <- as.Date(ltrm_gzsd$fdate, "%m/%d/%Y")
ltrm_gzsd <- ltrm_gzsd %>% mutate(year = year(fdate))
ltrm_gzsd <- ltrm_gzsd %>%
  filter(year != 2107)

##########################################
## Least Squares Estimation of surv_alpha in main channel
##########################################
# Filter main channel of UMR
ltrm_gzsd_main <- ltrm_gzsd %>% filter(pool %in% c("04", "08", "13", "26", "OR"))
# Filter pre-carp establishment year < 2000
# ltrm_gzsd_main <- ltrm_gzsd_main %>% filter(year < 2000)
# Round observations to nearest 10th
ltrm_gzsd_main <- ltrm_gzsd_main %>%
  mutate(length_round = round(ltrm_gzsd_main$length, -1))
# Bin data for comparison with model equilibrium
ltrm_gzsd_main <- ltrm_gzsd_main %>% 
  group_by(length_round) %>%
  summarize(count = n(), .groups = 'drop') %>%
  mutate(freq = count/sum(count)) %>%
  arrange(length_round)
# length mesh to use for model
zmesh_obs <- ltrm_gzsd_main$length_round 

## Least square function to minimize ##
least_sq <- function(x){
  # assign parameters where
  # x[1] = surv_alpha, x[2] = surv_beta
  m_par$surv_alpha <- x[1]
  m_par$surv_beta <- x[2] 
  ## Run model
  l_shad <- 0.00   # lower size limit in mm
  u_shad <- max(ltrm_gzsd_main$length_round)  # upper size limit in mm
  delta_z <- (u_shad - l_shad) / N
  zmesh <-  l_shad + ((1:N) - 1 / 2) * (u_shad - l_shad) / N
  # Initial length distribution
  n <- matrix(0, length(zmesh), tf)
  n0_total <- 995
  n[, 1] <- dnorm(zmesh, mean = 0.5 * m_par$grow_max, sd = 30)
  n[, 1] <- (n[, 1] / sum(n[, 1])) * n0_total / delta_z
  # Dynamical System
  for (i in 1:(tf - 1)) {
    k_iter <- (p_z1z(zmesh, zmesh, m_par) + f_z1z(zmesh, zmesh, n[, i],
                                                          m_par)) * delta_z
    n[, i + 1] <- k_iter %*% n[, i]
  }
  # Compute equilibrium frequencies - averaged over approximately 1 period of 8 years
  ### Since tf = 200:
    model_equil <- n[zmesh_obs/delta_z, tf]/sum(n[zmesh_obs/delta_z, tf])
  ## End model 
  sse <- sum((ltrm_gzsd_main$freq - model_equil/delta_z)^2)
  return(sse)
}

# Initial frame for model
N <- max(ltrm_gzsd_main$length_round) # number of size classes
l_shad <- 0   # lower size limit in mm
u_shad <- max(ltrm_gzsd_main$length_round)  
delta_z <- (u_shad - l_shad) / N
zmesh <-  l_shad + (1:N) * (u_shad - l_shad) / N
# Optimize step

# Initial guess surv_alpha = 90, surv_beta = -5
# opt <- optim(c(90,-50), fn = least_sq)
 
# Estimate surv_alpha and surv_beta for each year (8) of the period
# Use mean for model parameters

n_total_period <- 9

surv_params <- tibble(
  year = tf-1-n_total_period + 1:(n_total_period), 
  surv_alpha = rep(0, n_total_period),
  surv_beta = rep(0, n_total_period)
)

for(i in 91:99){
  tf <- i
  opt <- optim(c(90,-50), fn = least_sq)
  surv_params[i-90,2] <- opt$par[1] # assign alpha
  surv_params[i-90,3] <- opt$par[2] # assign beta
}

#### Now Check with La Grange
m_par$surv_alpha <- mean(surv_params$surv_alpha)
m_par$surv_beta <- mean(surv_params$surv_beta)
n <- matrix(0, length(zmesh), tf)
n0_total <- 995
n[, 1] <- dnorm(zmesh, mean = 0.5 * m_par$grow_max, sd = 30)
n[, 1] <- (n[, 1] / sum(n[, 1])) * n0_total / delta_z
# Dynamical System
for (i in 1:(tf - 1)) {
  k_iter <- (p_z1z(zmesh, zmesh, m_par) + f_z1z(zmesh, zmesh, n[, i],
                                                        m_par)) * delta_z
  n[, i + 1] <- k_iter %*% n[, i]
}

