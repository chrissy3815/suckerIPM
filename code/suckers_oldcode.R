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
