###########################################################################
### Starting to build up the ability to generate intermediate phenotypes
### Christina Hernandez
### November 2023
###########################################################################

library(here)

## Set the original parameters -------------------------------------------------

# Run the white suckers code, then save the parameters to a new object:
source(here("code", "suckers_modeling.R"))
WS_params<- m_par

# Do the same for the summer suckers code:
source(here("code", "summer_suckers_modeling.R"))
SS_params<- m_par

## Build the WS (white sucker) and SS (summer sucker) IPM models ---------------
###### Note that this currently will not work because the two models use different survival functions.


# We'll use the same model structure for both
m <- 300 # number of meshpoints: bins for the integration
L <- 0.00   # lower size limit in mm
U <- 600.00    # upper size limit in mm - must be larger than Linf
h <- (U-L)/m # integration bin width
meshpts <-  L + (1:m)*h - h/2

# White suckers IPM:
WS_Pmat<- h*p_z1z(meshpts, meshpts, WS_params)
WS_Fmat<- h*(f_z1z(meshpts, meshpts, WS_params))
# P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
# F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
WS_Kmat<- WS_Pmat+WS_Fmat

# Summer suckers IPM:
SS_Pmat<- h*p_z1z(meshpts, meshpts, SS_params)
SS_Fmat<- h*(f_z1z(meshpts, meshpts, SS_params))
# P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
# F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
SS_Kmat<- SS_Pmat+SS_Fmat
