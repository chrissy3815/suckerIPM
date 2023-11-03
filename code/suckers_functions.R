###########################################################################
### Functions for use in white suckers and summer suckers project
### Christina Hernandez
### Nov 2023
###########################################################################

# Function to pull out and plot survival-at-size:
surv_at_size_plot<- function(Pmat, sizes){
  # This function takes in a survival kernel, called Pmat, and a vector of the
  # centered (mean) value in each size bin of the kernel. Basically, sizes is
  # the meshpts.

  surv_at_size<- colSums(Pmat)
  plot(sizes, surv_at_size, type='l', lwd=2, col='blue',
       xlab='Length (mm)', ylab='Annual survival probability')
  abline(h=max(surv_at_size), col='grey')
  legend('bottomright', lty=1, col='grey',
         legend=paste('Max. survival =', round(max(surv_at_size),2)))

}

#Function to return the distribution of ages at death:
### Note to Carl: This has turned out to be more mathematically complicated than
#I realized. I can write a function to calculate this if you really want it, but
#maybe just knowing by what age 99% of the population has died is enough?

# Function to return the age at which X percentage of population is dead:
age_most_individuals_dead<- function(Pmat, Fmat, proportion=0.9){
  # calculate the stable size distribution:
  Kmat<- Pmat+Fmat #overall projection kernel
  eigz<- eigen(Kmat)
  stable_size_dist<- Re(eigz$vectors[,which(Re(eigz$values)==max(Re(eigz$values)))])
  # cohort of YOY:
  YOY_dist<- Fmat%*%stable_size_dist
  YOY_dist<- YOY_dist/sum(YOY_dist) #re-scale to length 1
  # initialize survival to age a:
  Pa<- diag(1, dim(Pmat)[1])
  PaC0<- Pa%*%YOY_dist
  # initialize age counter
  a<- 1
  while (sum(PaC0)>(1-proportion)){
    Pa<- Pmat %*% Pa
    PaC0<- Pa %*% YOY_dist
    a<- a+1
  }
  return(a)
}

# Function to return stable size distribution, without YOY:
stable_size_dist_noYOY<- function(Pmat, Fmat){
  # calculate the stable size distribution:
  Kmat<- Pmat+Fmat #overall projection kernel
  eigz<- eigen(Kmat)
  stable_size_dist<- Re(eigz$vectors[,which(Re(eigz$values)==max(Re(eigz$values)))])
  # rescale to sum to 1:
  stable_size_dist<- stable_size_dist/sum(stable_size_dist)

  # calculate the size distribution of YOY, as a cohort produced at stable dist:
  YOY_t1<- Fmat%*%stable_size_dist
  # project the stable distribution forward one time step:
  size_dist_t1<- Kmat%*%stable_size_dist

  # subtract YOY out:
  size_dist_noYOY<- size_dist_t1 - YOY_t1
  # rescale to sum to 1:
  size_dist_noYOY<- size_dist_noYOY/sum(size_dist_noYOY)

  return(size_dist_noYOY)
}

# Function to return stable size distribution, without individuals below a certain size:
stable_size_dist_sizeThreshold<- function(Pmat, Fmat, sizes, threshold=50){
  # Note: sizes is a vector, basically the meshpts

  # calculate the stable size distribution:
  Kmat<- Pmat+Fmat #overall projection kernel
  eigz<- eigen(Kmat)
  stable_size_dist<- Re(eigz$vectors[,which(Re(eigz$values)==max(Re(eigz$values)))])

  # subtract out the individuals below a certain size:
  I<- which(sizes<threshold)
  stable_size_dist[I]<- 0
  # rescale to sum to 1:
  stable_size_dist<- stable_size_dist/sum(stable_size_dist)

  return(stable_size_dist)

}
