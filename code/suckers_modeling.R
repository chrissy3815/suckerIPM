###########################################################################
### Modeling life history of white suckers and summer suckers
### Christina Hernandez
### July 2023
###########################################################################
# Load necessary packages:
library(readxl)

# Set working directory:
setwd("/Users/chrissy/Cornell_Postdoc/suckers_modeling/")

# Read in data:
observations<- read_xlsx("sucker_model_data.xlsx")
lifehistory<- read_xlsx("life_history_model.xlsx")
# get rid of blank rows in observations table:
I<- which(is.na(observations$year))
observations<- observations[-I,]
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

# Ideally, we'd want to use a von Bertalanffy growth model with an Lmax much
# bigger than the largest observed fish, and linear growth in the observed
# range.

# estimating Linf from the Ford-Walford method:
meanLAA<- aggregate(Len~age, data=femaleLH, FUN=mean)
Lt<- meanLAA$Len[1:(length(meanLAA$Len)-1)]
Lt1<- meanLAA$Len[2:length(meanLAA$Len)]
plot(Lt, Lt1, xlim=c(100, 700), ylim=c(100, 700))
Linf_model<- lm(Lt1~Lt)
abline(61.4445, 0.8665)
abline(0,1, lty=2)
# estimates:
# Linf is the point where the linear model and the 1:1 line intersect:
Linf<- 61.4445/(1-0.8665)
# K is -ln(slope) of the linear model:
K<- -log(0.8665)
# t0 can be estimated from one of the age-meanlength pairs:
I<- which(meanLAA$age==14)
t0<- meanLAA$age[I]*log((Linf-meanLAA$Len[I])/Linf)

vbStarts<- list(Linf=Linf, K=K, t0=t0)
vbTypical<-Len~Linf*(1-exp(-K*(age-t0)))

fitTypical<-nls(vbTypical,data=femaleLH[,c("age", "Len")],start=vbStarts)

plot(femaleLH$age, femaleLH$Len, xlim=c(1,20), ylim=c(100, 600))
exes<- 1:20
whys<- Linf*(1-exp(-K*(exes-t0)))
lines(exes, whys)


m_par <- tibble(
  ## growth from Michaletz (2017) paper - Table 1
  grow_rate = 0.26, # growth rate
  grow_max  =   394.3, # maximum length in mm: L_inf = 394.30
  grow_sd   =   25,  # growth sd (Max TL - L_inf)
  surv_min  =  0.002, # min survival - Bodola (1955)
  # Then et al (2015): surv_max computed from 
  # 1-natural mortality = 1 - 8.872*K^.73 L^-.33
  surv_max = 1 - 8.872*grow_rate^.73*grow_max^(-.33), 
  # inflection point: will be temp dependent
  # computed for La Grange Reach
  surv_alpha = 80.0136, 
  surv_beta = -139.9312, # slope
  ## New recruit from Michaletz (2017)
  recruit_mean = 105,
  recruit_sd = 25, # same as grow_sd
  ## From Bodola (1955):
  egg_viable = 0.002,
  ## Estimated from Jons and Miranda (1997)
  egg_slope = coef(egg_extended_nls)[2], # -4.361915
  egg_max =  coef(egg_extended_nls)[1], # 41540.608025
  egg_infl = coef(egg_extended_nls)[3], # 935.528239
  ## Spawning Probability - Estimated from Michaletz (2009)
  prob_spawn = 0.90,
  surv0_int = coef(surv_den_exp)[1], # 0.2686
  surv0_decay = coef(surv_den_exp)[2], # 0.0030
)



