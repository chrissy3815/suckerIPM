## Old code snippets that I don't want to delete yet.

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
# estimate K is -ln(slope) of the linear model:
K<- -log(0.8665)

## Fitting a survival model
size_breaks<- seq(100, 480, 20)
yr1_counts<- hist(EKL_2020$length, breaks = size_breaks)

# Project these individuals forward using von Bertalanffy curve:
yr2_predicted <- growth_params$Linf * (1 - exp(-growth_params$K)) +
  exp(-growth_params$K) * EKL_2020$length

yr2_predcounts<-hist(yr2_predicted, breaks=size_breaks)

# Compare the predicted size distribution to the observed distribution:
yr2_obscounts<- hist(EKL_2021$length, breaks = size_breaks)

par(mfrow=c(1,1))
plot(size_breaks[2:length(size_breaks)], yr2_predcounts$counts, type='l')
lines(meshpts, yr2_obscounts$counts, col='red')
lines(meshpts, yr1_counts$counts, col='blue')

# okay this doesn't work because more individuals were collected in yr2!

