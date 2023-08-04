#==========================================================
#SUPPLEMENTARY MATERIAL TO THE PAPER Fitness consequences of early life conditions and maternal size effects in a freshwater top predator
#by Vindenes et al., 2016 Journal of Animal Ecology.
#=========================================================
#R script for the IPM and corresponding elasticity/sensitivity analysis described in the main text


#MAIN STATE VARIABLES OF THE MODEL

#Dynamic trait: Body length x (cm)
#Static trait: Body length at age 1 y (cm)
#Environmental variable: Temperature  (annual mean surface temperature, in degrees Celsius)




#=========================================================
#SETUP
#=========================================================
 
library(MASS)
#--------------------------
#Define some functions
#--------------------------

#Calculate lambda, stable structure and reproductive value based on eigen-analysis.   
 
wvlambda <- function(Kmat){
 ev <- eigen(Kmat)
 tev <- eigen(t(Kmat))
 lmax <- which.max(Re(ev$values))
 W <- ev$vectors
 V <- tev$vectors
 w <- as.matrix(abs(Re(W[, lmax]))/sum(abs(Re(W[, lmax]))))
 w <- w/(sum(w))
 v <- as.matrix(abs(Re(V[, lmax])))
 v <- v/sum(w*v )
 v <- ifelse(w*v <= 0, 0, v)
 return(list("lambda"=max(Re(ev$values)),"w"=w,"v"=v))
}

#Alternative: Calculate lambda for a matrix K, using projection 
lambdafun <- function(K, tol=1e-10) {
	n <- length(K[1,])
	Nt <- rep(1,n)
	lam <- 1
	prev.lam <- .5
	while (abs(lam-prev.lam)>tol){
		prev.lam <- lam
		Nt.new <- K %*% Nt
		lam <- sum(Nt.new) / sum(Nt)
		Nt <- Nt.new
		}
	lam
	}
 
 

#--------------------------
#INITIALIZE THE MODEL 
#--------------------------


#Lower and upper limits for dynamic trait x (body length, cm)
Lx <- 10
Ux <- 130


 

#Number of "mesh points"
n <- 100
 
#Vector of dynamic trait values
x <- seq(Lx, Ux, length=n)

#Vector of static trait values (offspring length)
 

m <- max(which(x < 40))
y <- x[1:max(which(x < 40))]


#Vector of egg weights
w <- seq(0.001, 0.008, length=20)

#Vector of temperatures
Tvec <- seq(8, 13, length=100)

#Bin size
dx <- x[2]-x[1]
dy <-  y[2]-y[1]

#use midpoint rule? 
#x <- seq(Lx+dx/2,Ux-dx/2,dx)
#y <- seq(Ly+dy/2,Uy-dy/2,dy)

 
#==========================================
#VITAL RATE FUNCTIONS
#==========================================
#Four main functions: Survival, Transition, offspringnumber, Offspring distribution.

#"w" is an environmnetal effect, to be included in environmental stochasticity.

#******************************************
#offspring number (number of female 1-year olds produced per female)
#******************************************

#------------------------------------------
#FECUNDITY (EGG NUMBER)
#------------------------------------------
#Square root scale


Eggs.int <- -3.630672e+05
Eggs.x <- -8.148982e+01
Eggs.year <- 3.693226e+02
Eggs.year2 <- -9.400070e-02
Eggs.month <- 1.633226e+00 
Eggs.knsoma <- 7.726489e-01
Eggs.year.x <-  4.534849e-02
 

 
 
sd.eggnumber <- 0.06203222 #For stochastic model



mx <- function(x, Y=1982, M=11.4, K=100,  u=0 ){
	(Eggs.int + Eggs.x*x+    Eggs.year2*Y^2 +  Eggs.year*Y+ Eggs.month*M+ Eggs.knsoma*K +   Eggs.year.x *x*Y+   u)^2
}

 
 
#------------------------------------------
#Probability of maturity 
#------------------------------------------
#Values chosen so that about 50% are mature at 41.5 cm (Frost and Kipling, 1967), and 100 % at ~60 cm.
#Most (about 93%) Windermere pike are mature at age 2, very few at age 1 and some at age 3.  Maturity is largely determined by body size  (Frost and Kipling, 1967).
 

pm <- function(x, pm0=-20.75, pm1=.5 ){
	1/(1+exp(-pm0-pm1*x))
}
 
 

#------------------------------------------
#AVERAGE EGG WEIGHT
#------------------------------------------

 
EggW.int <- 6.067225  
EggW.x <-  -1.195475e-03 
EggW.temp <- 5.243521e-04
EggW.x2 <- -6.152693e-07    
EggW.year <- -6.077168e-03 
EggW.year2 <- 1.518297e-06 
EggW.month <- 4.512487e-04 
EggW.knsoma <- 9.564273e-06  
EggW.x.year <- 7.205317e-07  
EggW.x.temp <- -1.116501e-05
  
sd.eggweight <- 0.0003151809 
sdYeff.eggweight <- 8.775023e-07

wx <- function( x,y=23, Y= 1982, Tem=10.5, K =100, M =11.4,   u=0 ){
	EggW.int  + EggW.x * x + EggW.temp* Tem  +EggW.x2*x^2+EggW.year2*Y^2 +EggW.year*Y+EggW.x.year*Y*x+ EggW.month *M + EggW.knsoma *K+EggW.x.temp*Tem*x   +u
	}
 
#In appendix D:	
wx2 <- function( x, Tem=10.5, Y= 1982,  K =100, M =11.4,   u=0 ){
	#-0.005  + .00013  * x
	(-1.2e-3-EggW.x.temp*10.5*x-EggW.temp*10.5) + .00006  * x +EggW.x.temp*Tem*x + EggW.temp* Tem
	}

 
 
   

#------------------------------------------
#SURVIVAL OF EGGS TO AGE 1
#------------------------------------------


#Not known, probably high
sd.eggsurv<- .1
#With four different scenarios

egg.surv <- function(w.egg, Tem=10.5, u=0, scen="Tempsize"){
	segg.int <- 2.8e-4
	interc <- log(segg.int/(1-segg.int))
	if(scen=="Tempsize"){#Scenario 1, temperature/size negative interaction
		segg.weff <- 1500
		segg.temp <- .7
		segg.temp.w <- -130
		}
	if(scen=="Size"){#Scenario 1, Only size effect
		segg.weff <- 800 #150
		segg.temp <- 0
		segg.temp.w <- 0
		}
	if(scen=="Temp"){#Scenario 3, Only temp effect
		segg.weff <- 0
		segg.temp <- 0.6
		segg.temp.w <- 0
		}
	if(scen=="None"){#Scenario 4, No effect
		segg.weff <- 0
		segg.temp <- 0 
		segg.temp.w <- 0
		}
	1/(1+exp(-((interc-segg.temp*10.5-0.00354* segg.weff-segg.temp.w*10.5* 0.00354) +segg.weff*w.egg+ segg.temp* Tem+ segg.temp.w* Tem*w.egg+u)))
		}
	
#------------------------------------------
#PLOT FOR DIFFERENT SCENARIOS
#------------------------------------------
 

par(mfrow=c(2,2),bty="l",mar=c(3,4,2,2),las=0,family="Helvetica",cex=1)
 plot(w,egg.surv(w,Tem=10.5), xlab="", ylab="", type="l", lwd=2, ylim=c(0,.001),   main="Scenario 1: Interaction" )
 lines(w,egg.surv(w,Tem=12),lwd=2,col=2,lty=2)
  lines(w,egg.surv(w,Tem=9),lwd=2,col=4,lty=3)
    mtext(expression(paste("Offspring survival ", s[O](w,T))) ,2,2)
    mtext("Egg weight w (g)",1,2)
     mtext("A",3,at=.001,cex=1.5)
     abline(h=.00028,lty=3,lwd=2,col=8)
  legend(0.001,.001,c(expression(paste(T ," =9",degree,"C")),expression(paste(T ," =10.5",degree,"C")),expression(paste(T ," =12",degree,"C"))),lty=c(3,1,2), bty="n",lwd=2,col=c(4,1,2))
 plot(w,egg.surv(w,Tem=10.5, scen="Size"), xlab="", ylab="", type="l", lwd=2, ylim=c(0,.003),   main="Scenario 2: Egg weight" )
 lines(w,egg.surv(w,Tem=12, scen="Size"),lwd=2,col=2,lty=2)
  lines(w,egg.surv(w,Tem=9, scen="Size"),lwd=2,col=4,lty=3)
    mtext(expression(paste("Offspring survival ", s[O](w,T))) ,2,2)
    mtext("Egg weight w (g)",1,2)
     mtext("B",3,at=.001,cex=1.5)
          abline(h=.00028,lty=3,lwd=2,col=8)
    plot(w,egg.surv(w,Tem=10.5, scen="Temp"), xlab="", ylab="", type="l", lwd=2, ylim=c(0,.001),   main="Scenario 3: Temperature" )
 lines(w,egg.surv(w,Tem=12, scen="Temp"),lwd=2,col=2,lty=2)
  lines(w,egg.surv(w,Tem=9, scen="Temp"),lwd=2,col=4,lty=3)
    mtext(expression(paste("Offspring survival ", s[O](w,T))) ,2,2)
    mtext("Egg weight w (g)",1,2)
  mtext("C",3,at=.001,cex=1.5)   
       abline(h=.00028,lty=3,lwd=2,col=8) 
       plot(w,egg.surv(w,Tem=10.5, scen="None"), xlab="", ylab="", type="l", lwd=2, ylim=c(0,.001),   main="Scenario 4: Constant" )
 lines(w,egg.surv(w,Tem=12, scen="None"),lwd=2,col=2,lty=2)
  lines(w,egg.surv(w,Tem=9, scen="None"),lwd=2,col=4,lty=3)
    mtext(expression(paste("Offspring survival ", s[O](w,T))) ,2,2)
    mtext("Egg weight w (g)",1,2)
      mtext("D",3,at=.001,cex=1.5)
           abline(h=.00028,lty=3,lwd=2,col=8)
 
 
#-----------------------------------------
#TOTAL offspringnumber (1-YEAR OLD OFFSPRING)
#----------------------------------------
#Number of  (female) offspring produced per female, surviving to age 1

offspringnumber <- function(x,TP=10.5, TC=10.5, umx=0, uegg=0,uw=0,scen="Tempsize"){
	meanw <- wx(x=x,Tem=TP, u=uw)
	0.5*pm(x=x)*mx(x=x,u=umx )* egg.surv(w.egg=meanw, Tem=TC, u=uegg,scen=scen)
}
  

 
ylims<-c(0,130)
par(mfrow=c(4,3),mar=c(3,4,2,2),las=0,cex=1, bty="l")
plot(x, offspringnumber(x,TP=9,TC=10.5), type="l",lwd=2,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2),ylab="",main="T*=9" )
lines(x, offspringnumber(x,TP=9,TC=12),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP=9,TC=9),lwd=2,col=4,lty=3) 
mtext("A",3,at=10,cex=1.5)
 mtext("b(x,T,T*)",2,2)
 legend(10,130,c(expression(paste(T," =9",degree,"C")),expression(paste(T," =10.5",degree,"C")),expression(paste(T," =12",degree,"C"))),lty=c(3,1,2), bty="n",lwd=2,col=c(4,1,2))
 


plot(x, offspringnumber(x,TP=10.5,TC=10.5),type="l",xlab="Parent x", ylab="",lwd=2 ,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2),main="T*=10.5" )
lines(x, offspringnumber(x,TP= 10.5,TC=12),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP= 10.5,TC=9),lwd=2,col=4,lty=3) 
mtext("B",3,at=10,cex=1.5)
 
 
plot(x, offspringnumber(x,TP=12,TC=10.5),type="l",xlab="Parent x", ylab="", lwd=2 ,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2),main="T*=12" )
lines(x, offspringnumber(x,TP= 12,TC=12),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP= 12,TC=9),lwd=2,col=4,lty=3) 
mtext("C",3,at=10,cex=1.5)

 
plot(x, offspringnumber(x,TP=9,TC=10.5,scen="Size"),type="l",xlab="", ylab="" ,lwd=2,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2) )
lines(x, offspringnumber(x,TP=9,TC=12,scen="Size"),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP=9,TC=9,scen="Size"),lwd=2,col=4,lty=3) 
 mtext("b(x,T,T*)",2,2)
 mtext("D",3,at=10,cex=1.5)


plot(x, offspringnumber(x,TP=10.5,TC=10.5,scen="Size"),type="l",xlab="", ylab="", lwd=2 ,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2))
lines(x, offspringnumber(x,TP= 10.5,TC=12,scen="Size"),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP= 10.5,TC=9,scen="Size"),lwd=2,col=4) 
mtext("E",3,at=10,cex=1.5)

plot(x, offspringnumber(x,TP=12,TC=10.5,scen="Size"),type="l",xlab="", ylab="" ,lwd=2,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2))
lines(x, offspringnumber(x,TP= 12,TC=12,scen="Size"),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP= 12,TC=9,scen="Size"),lwd=2,col=4,lty=3) 
mtext("F",3,at=10,cex=1.5)

plot(x, offspringnumber(x,TP=9,TC=10.5,scen="Temp"),type="l",xlab="", ylab="", lwd=2,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2))
lines(x, offspringnumber(x,TP=9,TC=12,scen="Temp"),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP=9,TC=9,scen="Temp"),lwd=2,col=4,lty=3) 
 mtext("b(x,T,T*)",2,2)
 
 mtext("G",3,at=10,cex=1.5)

plot(x, offspringnumber(x,TP=10.5,TC=10.5,scen="Temp"),type="l",xlab="", ylab="" ,lwd=2 ,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2) )
lines(x, offspringnumber(x,TP= 10.5,TC=12,scen="Temp"),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP= 10.5,TC=9,scen="Temp"),lwd=2,col=4,lty=3) 
mtext("H",3,at=10,cex=1.5)

plot(x, offspringnumber(x,TP=12,TC=10.5,scen="Temp"),type="l",xlab="", ylab="", lwd=2,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2) )
lines(x, offspringnumber(x,TP= 12,TC=12,scen="Temp"),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP= 12,TC=9,scen="Temp"),lwd=2,col=4,lty=3) 
mtext("I",3,at=10,cex=1.5)

plot(x, offspringnumber(x,TP=9,TC=10.5,scen="None"),type="l",xlab="", ylab=" " ,lwd=2,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2) )
lines(x, offspringnumber(x,TP=9,TC=12,scen="None"),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP=9,TC=9,scen="None"),lwd=2,col=4,lty=3) 
mtext("J",3,at=10,cex=1.5)

 mtext("b(x,T,T*)",2,2)
 mtext(expression(paste("Length ", x, " (cm)")),1,2)
 
plot(x, offspringnumber(x,TP=10.5,TC=10.5,scen="None"),type="l",xlab=" ", ylab=" ", lwd=2 ,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2) )
lines(x, offspringnumber(x,TP= 10.5,TC=12,scen="None"),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP= 10.5,TC=9,scen="None"),lwd=2,col=4,lty=3) 
 mtext(expression(paste("Length ", x, " (cm)")),1,2)
mtext("K",3,at=10,cex=1.5)
 
 
plot(x, offspringnumber(x,TP=12,TC=10.5,scen="None"),type="l",xlab=" ", ylab=" ", lwd=2 ,ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,120,2) )
lines(x, offspringnumber(x,TP= 12,TC=12,scen="None"),lwd=2,col=2,lty=2)
lines(x, offspringnumber(x,TP= 12,TC=9,scen="None"),lwd=2,col=4,lty=3) 
 mtext(expression(paste("Length ", x, " (cm)")),1,2)
 mtext("L",3,at=10,cex=1.5)
 


 


#******************************************
#Offspring length distribution  NO HERITABILITY
#******************************************


y.int  <-  -63.43251402
y.temp <-  0.65294694
y.year <-  0.04063863 
 

sd.LAone <-  1.033711
#Res: 3.507127



xmean.off <- function(  Tem=10.5, Year=1965, u=0){
	y.int +y.temp*Tem +y.year*Year +u
}


 
 
 
#VARIANCE
xvar.off <- 13.19  
 

#DISTRIBUTION
off.dist.x <- function(xprime=x,Tem=10.5, u=0){
	mu <- xmean.off(Tem=Tem,u=u)
	var <- xvar.off 
	x2 <- dlnorm(xprime, meanlog=log(mu)-0.5*log(1+ var/(mu^2)), sdlog=sqrt(log(1+ var/(mu^2))))
	x2/sum(x2*dx)
}


#******************************************
#Offspring length distribution WITH HERITABILITY
#******************************************
herit.y <- .3 #Defines heritability (slope of parent/offspring regression)

xmean.off.inh <- function(Tem=10.5,  yparent=22, Year=1965, u=0){
	(y.int-herit.y*26.1)+yparent* herit.y+ y.temp* Tem +y.year*Year +u
}

xvar.off.inh <- 8 #Smaller because more variation in offspring length is accounted for by parental y

off.dist.x.inh <- function(xprime=x,yparent=24,Tem=10.5, u=0){
	mu <- xmean.off.inh(Tem=Tem,u=u,yparent=yparent)
	var <- xvar.off.inh
	x2 <- dlnorm(xprime, meanlog=log(mu)-0.5*log(1+ var/(mu^2)), sdlog=sqrt(log(1+ var/(mu^2))))
	x2/sum(x2*dx)
}




#Plot the mean plus some random realizations, assuming lognormal distribution: 
y2 <- seq(10,40,.1)
mu <- xmean.off()
var <- xvar.off 
mu2 <- xmean.off.inh(ypar=y2)
var2 <- xvar.off.inh

par(mfrow=c(1,2),bty="l", cex=1,las=1, mar=c(3,4,2,2))
plot(y2,rep(xmean.off(),length(y2)),type="l", lwd=2, ylim=c(5,45),xlim=c(5,45),xlab="",ylab="",main="No heritability" )
mtext("y parent",1,2,cex=1.2)
mtext("y' offspring",2,2, las=0,cex=1.2)
points(y2,rlnorm(n=length(y2), meanlog=log(mu)-0.5*log(1+ var/(mu^2)), sdlog=sqrt(log(1+ var/(mu^2)))), pch=19, cex=.2)
 
 
plot(y2,xmean.off.inh(yparent=y2),lwd=2, type="l",ylim=c(5,45),xlim=c(5,45),xlab="",ylab="",main="Heritability" )
points(y2,rlnorm(n=length(y2), meanlog=log(mu2)-0.5*log(1+ var2/(mu2^2)), sdlog=sqrt(log(1+ var2/(mu2^2)))),pch=19,cex=.2)
mtext("y parent",1,2,cex=1.2)
mtext("y' offspring",2,2, las=0,cex=1.2)
 
 

#******************************************
#SURVIVAL PROBABILITY
#****************************************** 

 
interc <- 73.258145314
xeff <- 0.489008352
teff <- 0.192934725
x2eff <- -0.003738462
yeff <- -0.043723032	
xteff <- -0.006844350



#------------------------------------------
#Effect of y 
#------------------------------------------
#Negative if cost of early growth, see Vindenes and Langangen 2015

alpha <- -0.05

 
sx1 <- function(x,y=23, Tem=10.5, u=0, Y=1972){
1/(1+exp(-(interc-alpha*23) -xeff*x- x2eff*x^2-teff* Tem-yeff*Y- xteff*x* Tem -alpha*y+u))		
	}
	
survival2 <- function(x, y=23, Tem=10.5, u=0){
 	xmax<-x[which(sx1(x=x, y=y, Tem = Tem,  u=u )==max(sx1(x=x,y=y, Tem = Tem,u=u)))]
 	ifelse(x<xmax,sx1(x=x,y=y, Tem = Tem,  u=u ),sx1(x=xmax,y=y, Tem = Tem,u=u))
 }
 
 
survival <-function(x2,y=23, Tem=10.5, u=0, b=1.8e-4){
 	 test <- survival2(x=x,y=y,Tem=Tem,u=u)
  	 test2 <- approxfun(seq(Lx, Ux, length=n), test)
  	 ifelse(x2<Ux, test2(x2), test2(x[n]))
  }
 
 
 
#******************************************
#SOMATIC GROWTH (DYNAMIC TRANSITIONS)
#****************************************** 

#Mean 

G.int <-  -1.006763e+02
G.x <- 2.789326e+00 
G.x2 <- -4.543838e-02 
G.x3 <- 4.586507e-04 
G.x4 <- -1.585666e-06 
G.L1 <- 3.710666e-01 
G.year<-  4.078231e-02
G.temp<-  1.318291e+00
G.x.temp <- -1.417239e-02
G.L1.x <- -4.098675e-03


 

sdint <- 0.5591958

sdres <- 3.352204

#Exponentially declining variance with rate 0.00632958


xmean <- function(x, Tem=10.5, y=23, u=0, Year=1966 ){
	new.x <- G.int+G.year*Year+G.x*x+G.x2*x^2+G.x3*x^3+G.x4*x^4+ G.temp* Tem +G.L1*y+G.L1.x*x*y+G.x.temp* Tem*x+u 
	ifelse(new.x < x, x, new.x)
}

 
xvar <- function(x, nu.var= -0.008106295,  sdres= 3.352204 ){
	sdres^2 *exp(nu.var*x)
}




Growth.x <- function(x, xprime, y=23, u=0, Tem=10.5){
	mu2 <- xmean(x=x, y=y, u=u, Tem=Tem)
	mu <- ifelse (mu2 < x, x, mu2)
	var <- xvar(x=x)
	z <- ifelse(xprime-x > 0, xprime, 0)
	y2 <- dlnorm(z, meanlog=log(mu)-0.5*log(1+ var/(mu^2)), sdlog=sqrt(log(1+ var/(mu^2))))
	if(sum(y2*dx)==0){
		return (c(rep(0,n-1),1/dx))
	}
	else return(y2/sum(y2*dx))
	y2/sum(y2*dx)
	} 
	
 
 
 
#******************************************
#Plot all vital rates (except offspring survival)
#******************************************
 
par(mfrow=c(2,3),mar=c(3,4,2,2),las=0,cex=1,bty="l",family="sans") 
plot(x, wx(x),type="l",lwd=2,xlab="", ylab="", col=1,ylim=c(0.001,.008),xlim=c(30,115), main="Egg weight")
lines(x,wx(x,Tem=12) , lwd=2,col=2,lty=2)
lines(x,wx(x,Tem=9) , lwd=2,col=4,lty=3)
legend(60,.0082,c(expression(paste("T* =9",degree,"C")),expression(paste("T*=10.5",degree,"C")),expression(paste( "T*=12",degree,"C"))),lty=c(3,1,2), bty="n",lwd=2,col=c(4,1,2))
mtext("Current length x (cm)",1,2,las=0,cex=1)
   mtext("Egg weight (g)",2,2)
  mtext("A",3,at=30,cex=1.2)
  
 plot(x, mx(x ),type="l",lwd=2,xlab="", ylab="", col=1, main="Fecundity",ylim=c(0,5e5), xlim=c(30,115))
mtext("Current length x (cm)",1,2,las=0,cex=1)
   mtext("Egg number",2,2)
 mtext("A",3,at=30,cex=1.2)

 plot(x, pm(x),type="l",lwd=2,xlab="", ylab="", col=1, main="Probability of maturity",ylim=c(0,1.02), xlim=c(30,60) )
mtext("Current length x (cm)",1,2,las=0,cex=1)
 mtext(expression("Maturity"),2,2)
  mtext("B",3,at=30,cex=1.2)
  

  
plot(x,off.dist.x(x) ,type="l",lwd=2,xlab="", ylab="", col=1, main="Offspring length dist",ylim=c(0,.12) ,xlim=c(10,50) )
lines(x, off.dist.x(x,Tem=12) , lwd=2,col=2,lty=2)
lines(x, off.dist.x(x,Tem=9) , lwd=2,col=4,lty=3)
legend(25,.12,c(expression(paste("T =9",degree,"C")),expression(paste("T=10.5",degree,"C")),expression(paste( "T=12",degree,"C"))),lty=c(3,1,2), bty="n",lwd=2,col=c(4,1,2))
mtext("Offspring length x' (cm)",1,2,las=0,cex=1)
   mtext("Density",2,2)
   mtext("D",3,at=10,cex=1.2)
   
plot(x, survival(x) ,type="l",lwd=2, xlab="", ylab="",cex=.5,col=1,ylim=c(-0.02,1.02),yaxp=c(0,1,2),xlim=c(10,120),xaxp=c(10,120,2),main="Survival probability")
lines(x, survival(x,Tem=12),lwd=2,col=2,lty=2)
lines(x, survival(x,Tem=9),lwd=2,col=4,lty=3)
legend(60,.5,c(expression(paste(T ," =9",degree,"C")),expression(paste(T ," =10.5",degree,"C")),expression(paste(T ," =12",degree,"C"))),lty=c(3,1,2), bty="n",lwd=2,col=c(4,1,2))
  mtext("E",3,at=10,cex=1.2)
mtext("Current length x (cm)",1, 2,las=0,cex=1)
   mtext("Survival probability",2,2)
   
plot(x, xmean(x) ,type="l",lwd=2, xlab="", ylab="",cex=.3,col=1,ylim=c(10,120),yaxp=c(10,120,2),xlim=c(10,120),xaxp=c(10,120,2),main="Mean growth")
lines(x,x,lty=3,col=8,lwd=2)
 
lines(x, xmean(x,T=12),col=2,lty=2,lwd=2)
lines(x, xmean(x,T=9),col=4,lty=3,lwd=2)
lines(0:130,0:130,lty=3)
mtext("Current length x (cm)",1, 2,las=0,cex=1)
   mtext("Next length (cm)",2,2)
 mtext("F",3,at=10,cex=1.2)


#==========================================
#PROJECTION KERNEL (full  model)
#==========================================

 

K.kernel <- function(Temp=10.5 , scen="Tempsize",inh=F){
	Sarray <- Barray <- array(0,dim=c(n,n,m,m))
	wmat <- diag(1,m,m)
	wmat2 <- diag(1,n,m)
	bxvec <- offspringnumber (x=x,TP=Temp,TC=Temp,  scen=scen )
	offvec.y <- off.dist.x(xprime=y,Tem=Temp )*dy
	for(l in 1:m){
		if(inh==T){
			offvec.y <- off.dist.x.inh(xprime=y, yparent=y[l],Tem=Temp )*dy
			}
		sxvec <- survival (x=x, y=y[l],Tem=Temp )
		for(j in 1:n){
				for(k in 1:m){
				Sarray[,j,k,l] <- Growth.x(x=x[j],xprime=x,y=y[l],Tem=Temp)* sxvec[j]*wmat[k,l]*dx
				for(i in 1:n){
					Barray[i,j,k,l] <- offvec.y[k]* bxvec[j]*wmat2[i,k]
					}
				}
			}
	}
	K.array <- Sarray + Barray
	Kmat <-  matrix(0,n*m,n*m)
	for(i in 1:m){
		for(j in 1:m){
			Kmat[(n*(i-1)+1):(i*n),(n*(j-1)+1):(j*n)] <- K.array[,,i,j]
			}
		}
	Kmat
	}
	
 
 

KS1 <- K.kernel(Temp=10.5,scen="Tempsize") 
resS1 <- wvlambda(KS1)
uS1 <- resS1$w
 vS1 <- resS1$v
resS1$lam
 
 
lambdatemp <- matrix(NA,nrow=4,ncol=7)
Tvec2 <- seq(9,12,length=7)
for(i in 1:7){
	lambdatemp[1,i]<-lambdafun(K.kernel(Temp=Tvec2[i]))
	lambdatemp[2,i]<-lambdafun(K.kernel(Temp=Tvec2[i],scen="Size") )
	lambdatemp[3,i]<-lambdafun(K.kernel(Temp=Tvec2[i],scen="Temp") )
	lambdatemp[4,i]<-lambdafun(K.kernel(Temp=Tvec2[i],scen="None") )
}
 
par(mfrow=c(1,1),bty="l",mar=c(3,4,2,2),las=0)
plot(Tvec2, lambdatemp[1,],type="l", lwd=2,lty=4,xlab="",ylab="",ylim=c(0.9,1.35),xlim=c(8.9,12.1),  main="Expected population growth rate")
abline(h=1,lty=3,col=8,lwd=2)
 lines(Tvec2, lambdatemp[2,], lwd=2,lty=2)
 lines(Tvec2, lambdatemp[3,], lwd=2,lty=3)
 lines(Tvec2, lambdatemp[4,], lwd=2,lty=1)
abline(v=9,lty=1,col=4,lwd=2)
abline(v=10.5,lty=1,col=3,lwd=2)
abline(v=12,lty=1,col=2,lwd=2)
mtext(expression(lambda),2,2.5,cex=1.2)
mtext(expression(paste("T(",degree,"C)")) ,1,2,cex=1.2)
legend(9,1.3,c("Scenario 1","Scenario 2","Scenario 3","Scenario 4"),lty=c(4,2,3,1), bty="n",lwd=2)

 

 

KS1 <- K.kernel(Temp=10.5,scen="Tempsize") 
KS2 <- K.kernel(Temp=10.5,scen="Size") 
KS3 <- K.kernel(Temp=10.5,scen="Temp") 
KS4 <- K.kernel(Temp=10.5,scen="None") 

resS1 <- wvlambda(KS1)
resS2 <- wvlambda(KS2)
resS3 <- wvlambda(KS3)
resS4 <- wvlambda(KS4)

KS19 <- K.kernel(Temp=9,scen="Tempsize") 
KS29 <- K.kernel(Temp=9,scen="Size") 
KS39 <- K.kernel(Temp=9,scen="Temp") 
KS49 <- K.kernel(Temp=9,scen="None") 
resS19 <- wvlambda(KS19)
resS29 <- wvlambda(KS29)
resS39 <- wvlambda(KS39)
resS49 <- wvlambda(KS49)

KS112 <- K.kernel(Temp=12,scen="Tempsize") 
KS212 <- K.kernel(Temp=12,scen="Size") 
KS312 <- K.kernel(Temp=12,scen="Temp") 
KS412 <- K.kernel(Temp=12,scen="None") 
resS112 <- wvlambda(KS112)
resS212 <- wvlambda(KS212)
resS312 <- wvlambda(KS312)
resS412 <- wvlambda(KS412)


 
 
lamS1 <- resS1$lam
lamS2 <- resS2$lam
lamS3 <- resS3$lam
lamS4 <- resS4$lam

uS1 <- resS1$w
uS2 <- resS2$w
uS3 <- resS3$w
uS4 <- resS4$w
vS1 <- resS1$v
vS2 <- resS2$v
vS3 <- resS3$v
vS4 <- resS4$v

#Stable structure value matrixes in different scenarios, at three temperatures
umatS1 <- umatS2 <- umatS3 <- umatS4 <- matrix(NA, nrow=n,ncol=m)
for (i in 0:(m-1)){
	umatS1[,i+1]<- uS1[(i*n+1):((i+1)*n)]
	umatS2[,i+1]<- uS2[(i*n+1):((i+1)*n)]
	umatS3[,i+1]<- uS3[(i*n+1):((i+1)*n)]
	umatS4[,i+1]<- uS4[(i*n+1):((i+1)*n)]
  
}


lamS19 <- resS19$lam
lamS29 <- resS29$lam
lamS39 <- resS39$lam
lamS49 <- resS49$lam
#uS1 <- resS1$w
vS19 <- resS19$v
vS29 <- resS29$v
vS39 <- resS39$v
vS49 <- resS49$v


lamS112 <- resS112$lam
lamS212 <- resS212$lam
lamS312 <- resS312$lam
lamS412 <- resS412$lam
#uS1 <- resS1$w
vS112 <- resS112$v
vS212 <- resS212$v
vS312 <- resS312$v
vS412 <- resS412$v

#Reproductive value matrixes in different scenarios, at three temperatures
vmatS1 <- vmatS2 <- vmatS3 <- vmatS4 <- vmatS19 <- vmatS29 <- vmatS39 <- vmatS49 <- vmatS112 <- vmatS212 <- vmatS312 <- vmatS412 <-matrix(NA, nrow=n,ncol=m)
for (i in 0:(m-1)){
	vmatS1[,i+1]<- vS1[(i*n+1):((i+1)*n)]
	vmatS19[,i+1]<- vS19[(i*n+1):((i+1)*n)]	
	vmatS112[,i+1]<- vS112[(i*n+1):((i+1)*n)]	
	vmatS2[,i+1]<- vS2[(i*n+1):((i+1)*n)]
	vmatS29[,i+1]<- vS29[(i*n+1):((i+1)*n)]	
	vmatS212[,i+1]<- vS212[(i*n+1):((i+1)*n)]	
	vmatS3[,i+1]<- vS3[(i*n+1):((i+1)*n)]
	vmatS39[,i+1]<- vS39[(i*n+1):((i+1)*n)]	
	vmatS312[,i+1]<- vS312[(i*n+1):((i+1)*n)]	
	vmatS4[,i+1]<- vS4[(i*n+1):((i+1)*n)]
	vmatS49[,i+1]<- vS49[(i*n+1):((i+1)*n)]	
	vmatS412[,i+1]<- vS412[(i*n+1):((i+1)*n)]	 
}

 
 #Reproductive value at age 1 as a function of y, in different scenarios, at three temperatures
v0S1 <- v0S19 <-  v0S112 <- v0S2 <- v0S29 <-  v0S212 <- v0S3 <- v0S39 <-  v0S312<- v0S4 <- v0S49 <-  v0S412 <- rep(NA,m)

for(i in 1:m){
	v0S1[i] <- vmatS1[i,i]
	v0S19[i] <- vmatS19[i,i]
	v0S112[i] <- vmatS112[i,i]
	v0S2[i] <- vmatS2[i,i]
	v0S29[i] <- vmatS29[i,i]
	v0S212[i] <- vmatS212[i,i]
	v0S3[i] <- vmatS3[i,i]
	v0S39[i] <- vmatS39[i,i]
	v0S312[i] <- vmatS312[i,i]
	v0S4[i] <- vmatS4[i,i]
	v0S49[i] <- vmatS49[i,i]
	v0S412[i] <- vmatS412[i,i]
}
par(mfrow=c(2,2),bty="l",mar=c(3,4,2,2),las=0,cex=1)
ylims<-c(0,2)
plot(y, v0S1,lwd=2,type="l",main="Scenario 1: Interaction", xlab="", ylab="", ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,34,2)) 
lines(y, v0S19,lwd=2,lty=3,col=4) 
lines(y, v0S112,lwd=2,lty=2,col=2)
mtext(expression(paste("Reproductive value ",v[1](y))),2,2.5,cex=1 )
mtext( "Offspring length y (cm)" ,1,2,cex=1)
legend(10,1.8,c(expression(paste(T,"=9",degree,"C")),expression(paste(T,"=10.5",degree,"C")),expression(paste(T,"=12",degree,"C"))),lwd=2,lty=c(3,1,2),col=c(4,1,2),bty="n")
mtext("A",3,at=10,cex=1.5)

plot(y, v0S2,lwd=2,type="l",main="Scenario 2: Egg weight", xlab="", ylab="", ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,34,2)) 
lines(y, v0S29,lwd=2,lty=3,col=4) 
lines(y, v0S212,lwd=2,lty=2,col=2)
 mtext("B",3,at=10,cex=1.5)
mtext(expression(paste("Reproductive value ",v[1](y))),2,2.5,cex=1 )
mtext( "Offspring length y (cm)" ,1,2,cex=1)
plot(y, v0S3,lwd=2,type="l",main="Scenario 3: Temperature", xlab="", ylab="", ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,34,2)) 
lines(y, v0S39,lwd=2,lty=3,col=4) 
lines(y, v0S312,lwd=2,lty=2,col=2)
mtext(expression(paste("Reproductive value ",v[1](y))),2,2.5,cex=1 )
mtext( "Offspring length y (cm)" ,1,2,cex=1)
mtext("C",3,at=10,cex=1.5)

plot(y, v0S4,lwd=2,type="l",main="Scenario 4: Constant", xlab="", ylab="", ylim=ylims,yaxp=c(0,ylims[2],2),xaxp=c(10,34,2)) 
lines(y, v0S49,lwd=2,lty=3,col=4) 
lines(y, v0S412,lwd=2,lty=2,col=2) 
mtext(expression(paste("Reproductive value ",v[1](y))),2,2.5,cex=1 )
mtext( "Offspring length y (cm)" ,1,2,cex=1)
  mtext("D",3,at=10,cex=1.5)

par(mfrow=c(1,1),bty="l",mar=c(3,4,2,2),las=0,cex=1)
plot(x,apply(vmatS1,1,mean),type="l",xlab="x",ylab="v(x)",lwd=2) 
lines(x,apply(vmatS2,1,mean),lty=2,lwd=2) 
lines(x,apply(vmatS3,1,mean),lty=3,lwd=2)  
lines(x,apply(vmatS4,1,mean),lty=4,lwd=2) 

#***************************
#Sensitivity/Elasticity kernels
#***************************

SensmatS1 <- vS1 %*% t(uS1)
SensmatS2 <- vS2 %*% t(uS2)
SensmatS3 <- vS3 %*% t(uS3)
SensmatS4 <- vS4 %*% t(uS4)

EmatS1 <- SensmatS1*KS1/lamS1
EmatS2 <- SensmatS2*KS2/lamS2
EmatS3 <- SensmatS3*KS3/lamS3
EmatS4 <- SensmatS4*KS4/lamS4

#Integrate over y
SvecS1<- apply(SensmatS1, 1,sum)
KvecS1<- apply(KS1, 1,sum)
SvecS2<- apply(SensmatS2, 1,sum)
SvecS3<- apply(SensmatS3, 1,sum)
SvecS4<- apply(SensmatS4, 1,sum)

 
EvecS1<- apply(EmatS1, 1,sum)
EvecS2<- apply(EmatS2, 1,sum)
EvecS3<- apply(EmatS3, 1,sum)
EvecS4<- apply(EmatS4, 1,sum)

#Elasticity matrix across x and y 

EmatSumS1 <-EmatSumS2 <-EmatSumS3 <-EmatSumS4 <- matrix(NA, nrow=n,ncol=m)
for (i in 0:(m-1)){
	EmatSumS1[,i+1]<- EvecS1[(i*n+1):((i+1)*n)]
	EmatSumS2[,i+1]<- EvecS2[(i*n+1):((i+1)*n)]
	EmatSumS3[,i+1]<- EvecS3[(i*n+1):((i+1)*n)]
	EmatSumS4[,i+1]<- EvecS4[(i*n+1):((i+1)*n)]
}

 
SmatSumS1 <-SmatSumS2 <-SmatSumS3 <-SmatSumS4 <- Kmatsum <- matrix(NA, nrow=n,ncol=m)
for (i in 0:(m-1)){
	Kmatsum[,i+1]<- KvecS1[(i*n+1):((i+1)*n)]
	SmatSumS1[,i+1]<- SvecS1[(i*n+1):((i+1)*n)]
	SmatSumS2[,i+1]<- SvecS2[(i*n+1):((i+1)*n)]
	SmatSumS3[,i+1]<- SvecS3[(i*n+1):((i+1)*n)]
	SmatSumS4[,i+1]<- SvecS4[(i*n+1):((i+1)*n)]
}
 
par(mfrow=c(1,2),bty="l",mar=c(3,4,2,2),las=0,cex=1)
plot(x,apply(vmatS1,1,mean),type="l",xlab="x",ylab="v(x)",lwd=2) 
lines(x,apply(vmatS2,1,mean),lty=2,lwd=2) 
lines(x,apply(vmatS3,1,mean),lty=3,lwd=2)  
lines(x,apply(vmatS4,1,mean),lty=4,lwd=2) 

plot(x,apply(SmatSumS1,1,mean),type="l",xlab="x",ylab="v(x)",lwd=2) 
lines(x,apply(SmatSumS2,1,mean),lty=2,lwd=2) 
lines(x,apply(SmatSumS3,1,mean),lty=3,lwd=2)  
lines(x,apply(SmatSumS4,1,mean),lty=4,lwd=2) 


 
par(mfrow=c(1,2))
persp(x,y,SmatSumS1,theta=-60,phi=10)
persp(x,y,EmatSumS1,theta=30,phi=10)

 par(mfrow=c(1,2), mar=c(4,4,2,4),las=0,cex=1,bty="l",xaxs="i",yaxs="i")
plot(x,apply(SmatSumS1,1,sum),type="l",bty="l",lwd=2,ylab="",xlab=" ",main="Sensitivity",xlim=c(10,125))  
  plot(x,apply(EmatSumS1,1,sum),type="l",bty="l",lwd=2,ylab="",xlab=" ",main="Elasticity",xlim=c(10,125),ylim=c(0,.025))  
    
par(mfrow=c(1,1), mar=c(4,4,2,4),las=0,cex=1,bty="l",xaxs="i",yaxs="i")
plot(x,apply(EmatSumS1,1,sum),type="l",bty="l",lwd=2,ylab="",xlab=" ",main="",xlim=c(10,125),ylim=c(0,.025)) 
lines(x,apply(EmatSumS2,1,sum),lty=2,lwd=2)
lines(x,apply(EmatSumS3,1,sum),lty=2,lwd=2)
lines(x,apply(EmatSumS4,1,sum),lty=4,lwd=2)
par(new=T)
plot(x,log(apply(umatS1,1,sum)),lwd=2 , lty=5,col=3,yaxt="n",xlab="",ylab="",type="l",ylim=c(-13,-2),xlim=c(10,125)) 
axis(4,at=c(-13,-11,-9,-7,-5,-3,-1),tick=T,col=3,col.axis=3 )
mtext(expression(paste("Elasticity of ", lambda)),2,2.5,cex=1.2)
mtext("Log (u(x))",4,2.5,cex=1.2,col=3)
mtext( "Length x (cm)" ,1,2,cex=1.2)
legend(45, y=-8,c("Scenario 1","Scenario 2","Scenario 3","Scenario 4", "Stable distribution of x"),lty=1:5, col=c(rep(1,4),3),text.col=c(rep(1,4),3),bty="n",lwd=2)
 
 
par(mfrow=c(1,1), mar=c(4,4,2,4),las=0,cex=1,bty="l",xaxs="i",yaxs="i")
plot(x,apply(EmatSumS1,1,sum),type="l",bty="l",lwd=2,ylab="",xlab=" ",main="",xlim=c(10,125),ylim=c(0,.025)) 
lines(x,apply(EmatSumS2,1,sum),lty=2,lwd=2)
lines(x,apply(EmatSumS3,1,sum),lty=2,lwd=2)
lines(x,apply(EmatSumS4,1,sum),lty=4,lwd=2)
par(new=T)
plot(x,(apply(umatS1,1,sum)),lwd=2 , lty=5,col=3,yaxt="n",xlab="",ylab="",type="l",ylim=c(0,.13),xlim=c(10,125)) 
axis(4,tick=T,col=3,col.axis=3 )
mtext(expression(paste("Elasticity of ", lambda)),2,2.5,cex=1.2)
mtext("u(x)",4,2.5,cex=1.2,col=3)
mtext( "Length x (cm)" ,1,2,cex=1.2)
legend(45, y=.06,c("Scenario 1","Scenario 2","Scenario 3","Scenario 4", "Stable distribution of x"),lty=1:5, col=c(rep(1,4),3),text.col=c(rep(1,4),3),bty="n",lwd=2)
 
#===============================================================
#SENSITIVITY/ELASTICITY ANALYSIS OF LAMBDA TO UNDERLYING PARAMETERS
#===============================================================




#-----------------
#TO LENGTH x
#-----------------
 
 

Sensit.x <- function(x, Tem=10.5, scen="Tempsize", perturb=1e-5,res=resS1 ){
umat <-  vmat <- matrix(NA, nrow=n,ncol=m)
for (i in 0:(m-1)){
	umat[,i+1]<- res$w[(i*n+1):((i+1)*n)]
	vmat[,i+1]<- res$v[(i*n+1):((i+1)*n)]
	}
uy <- apply(umat,2,sum)
ux <- apply(umat,1,sum)
vx<-apply(vmat,1,sum)
vy<-apply(vmat,2,sum)
dw.dx <- (wx(x+ perturb,Tem=Tem)-wx(x,Tem=Tem))/perturb
wvec <- wx(x,Tem=Tem)
dm.dx <- (mx(x+perturb)-mx(x))/perturb
mxvec <- mx(x)
dpm.dx <- (pm(x+ perturb)-pm(x))/perturb
pmvec<- pm(x)
dse.dx <- (egg.surv(wx(x+ perturb,Tem=Tem),scen=scen,Tem=Tem)-egg.surv(wx(x,Tem=Tem),scen=scen,Tem=Tem))/perturb
dse.dw <- (egg.surv(wx(x,Tem=Tem)+ perturb,scen=scen,Tem=Tem)-egg.surv(wx(x,Tem=Tem),scen=scen,Tem=Tem))/perturb
sevec<- egg.surv(wx(x,Tem=Tem),scen=scen,Tem=Tem)
fxvec<- off.dist.x(xprime=x,Tem=Tem)
dfx.dx<- (off.dist.x(xprime=x+ perturb,Tem=Tem )-off.dist.x(xprime=x,Tem=Tem ))/perturb
wmat <- diag(1,m,m)
wmat2 <- diag(1,n,m)
dSxy.dx  <- sxymat<- matrix(0,n,m)#x,y
for(j in 1:m) {#y
	dSxy.dx[,j]<- (survival(x+ perturb,y=y[j],Tem=Tem)-survival(x,y=y[j],Tem=Tem))/perturb
	sxymat[,j]<-survival(x,y=y[j],Tem=Tem)
  }
dGxy.dx  <- gxymat<- array(NA,dim=c(n,n,m))#x,x',y
for(j in 1:m) {#y
	for(i in 1:n){#x
		gxymat[i,,j] <- Growth.x(x=x[i],xprime=x,y=y[j],Tem=Tem)#x,x',y
		dGxy.dx[i,,j] <- (Growth.x(x=x[i]+ perturb,xprime=x,y=y[j],Tem=Tem)-Growth.x(x=x[i],xprime=x,y=y[j],Tem=Tem))/perturb
		}
  }
MatM<- MatPM<- MatW <- MatG<- MatS<- array(NA,dim=c(n,n,m,m))#x,x',y,y'
for(i in 1:m){#y
	for(j in 1:m){#y'
		for (l in 1:n){#x'
			MatS[,l,i,j] <- vmat[l,j]*umat[,i]* dSxy.dx[,i]* gxymat[,l,i]*wmat[i,j]*dx
			MatG[,l,i,j] <- vmat[l,j]*umat[,i]* sxymat[,i]* dGxy.dx[,l,i]*wmat[i,j]*dx
			MatM[,l,i,j] <- vmat[l,j]*umat[,i]* .5*dm.dx* sevec* pmvec* fxvec[l]*wmat2[l,j]*dx 
			MatPM[,l,i,j] <- vmat[l,j]*umat[,i]* .5*dpm.dx* sevec* mxvec* fxvec[l]*wmat2[l,j]*dx 
			MatW[,l,i,j] <- vmat[l,j]*umat[,i]* .5*dw.dx*dse.dw* pmvec* mxvec* fxvec[l]*wmat2[l,j]*dx 
			}
		}
	}
list("M"=MatM, "S"=MatS,"G"= MatG,"W"=MatW,"PM"=MatPM)
}	
sensdec1 <- Sensit.x(x)
sensdec2 <- Sensit.x(x,scen="Size",res=resS2)
sensdec3 <- Sensit.x(x,scen="Temp",res=resS3)
sensdec4 <- Sensit.x(x,scen="None",res=resS4)


#SUMs
Sensit.x.Sum <- function(sensdec= sensdec1, Tem=10.5, res=resS1){
	Sensum<- Elastsum<- Elastsum2<- rep(NA,5)
	for(i in 1:5){
		Sensum[i]<- sum(apply(sensdec[[i]],c(1,2),sum))
		Elastsum[i]<- sum(apply(sensdec[[i]],c(1,2),sum)*x/res$lam)
	}
	names(Sensum) <- names(Elastsum) <- c("M", "S","G","W", "PM" )
  list("Sensitivities"=Sensum, "Elasticities"=Elastsum)
	}
	
sensdecsum1 <- Sensit.x.Sum()
sensdecsum2 <- Sensit.x.Sum(sensdec=sensdec2,Tem=10.5,res=resS2)
sensdecsum3 <- Sensit.x.Sum(sensdec=sensdec3,Tem=10.5,res=resS3)
sensdecsum4 <- Sensit.x.Sum(sensdec=sensdec4,Tem=10.5,res=resS4)

	


#-----------------
#TO TEMPERATURE T 
#----------------- 

Sensit.Temp <- function(x,Tem=10.5,scen="Tempsize",perturb=1e-5,res=resS1){
umat <-  vmat <- matrix(NA, nrow=n,ncol=m)
for (i in 0:(m-1)){
	umat[,i+1]<- res$w[(i*n+1):((i+1)*n)]
	vmat[,i+1]<- res$v[(i*n+1):((i+1)*n)]
	}
dw.dT <- (wx(x,Tem=Tem+ perturb)-wx(x,Tem=Tem))/perturb
wvec <- wx(x,Tem=Tem)
mxvec <- mx(x)
pmvec<- pm(x)
dse.dw <- (egg.surv(wx(x,Tem=Tem)+ perturb,scen=scen,Tem=Tem)-egg.surv(wx(x,Tem=Tem),scen=scen,Tem=Tem))/perturb
dse.dT <- (egg.surv(wx(x,Tem=Tem+ perturb),scen=scen,Tem=Tem+ perturb)-egg.surv(wx(x,Tem=Tem),scen=scen,Tem=Tem))/perturb
dse.dT2 <- (egg.surv(wx(x,Tem=Tem),scen=scen,Tem=Tem+ perturb)-egg.surv(wx(x,Tem=Tem),scen=scen,Tem=Tem))/perturb
sevec<- egg.surv(wx(x),scen=scen,Tem=Tem)
fxvec<- off.dist.x(xprime=x,Tem=Tem)
dfx.dT<- (off.dist.x(xprime=x,Tem=Tem+ perturb )-off.dist.x(xprime=x,Tem=Tem ))/perturb
wmat <- diag(1,m,m)
wmat2 <- diag(1,n,m)
dSxy.dT  <- sxymat<- matrix(0,n,m)#x,y
for(j in 1:m) {#y
	dSxy.dT[,j]<- (survival(x,y=y[j],Tem=Tem+ perturb)-survival(x,y=y[j],Tem=Tem))/perturb
	sxymat[,j]<-survival(x,y=y[j],Tem=Tem)
  }
dGxy.dT  <- gxymat<- array(NA,dim=c(n,n,m))#x,x',y
for(j in 1:m) {#y
	for(i in 1:n){#x
		gxymat[i,,j] <- Growth.x(x=x[i],xprime=x,y=y[j],Tem=Tem)#x,x',y
		dGxy.dT[i,,j] <- (Growth.x(x=x[i],xprime=x,y=y[j],Tem=Tem+ perturb)-Growth.x(x=x[i],xprime=x,y=y[j],Tem=Tem))/perturb
		}
  }
MatFX<- MatSE <- MatG<- MatS<-MatSE2<-MatW<- array(NA,dim=c(n,n,m,m))#x,x',y,y'
for(i in 1:m){#y
	for(j in 1:m){#y'
		for (l in 1:n){#x'
			MatS[,l,i,j] <- vmat[l,j]*umat[,i]* dSxy.dT[,i]* gxymat[,l,i]*wmat[i,j]*dx 
			MatG[,l,i,j] <- vmat[l,j]*umat[,i]* sxymat[,i]* dGxy.dT[,l,i]*wmat[i,j]*dx 
			MatSE[,l,i,j] <- vmat[l,j]*umat[,i]* .5*dse.dT* pmvec* mxvec* fxvec[l]*wmat2[l,j]*dx 
			MatSE2[,l,i,j] <- vmat[l,j]*umat[,i]* .5*dse.dT2* pmvec* mxvec* fxvec[l]*wmat2[l,j]*dx 
			MatW[,l,i,j] <- vmat[l,j]*umat[,i]* .5* dse.dw*dw.dT* pmvec* mxvec* fxvec[l]*wmat2[l,j]*dx 
			MatFX[,l,i,j] <- vmat[l,j]*umat[,i]* sevec* pmvec* mxvec* dfx.dT[l]*wmat2[l,j]*dx 
			}
		}
	}
list("S"=MatS,"G"= MatG,"SE"=MatSE,"FX"=MatFX,"W"=MatW,"SE2"=MatSE2)
}
senstemp1 <- Sensit.Temp(x)
senstemp2 <- Sensit.Temp(x,scen="Size",res=resS2)
senstemp3 <- Sensit.Temp(x,scen="Temp",res=resS3)
senstemp4 <- Sensit.Temp(x,scen="None",res=resS4)


Sensit.Temp.Sum <- function(senstemp=senstemp1,Tem=10.5,res=resS1){
	Sensum<- Elastsum<- rep(NA,6)
	for(i in 1:6){
		Sensum[i]<- sum(apply(senstemp[[i]],c(1,2),sum))
		Elastsum[i]<- Sensum[i]*Tem/res$lam
	}
	names(Sensum) <- names(Elastsum) <- c("S", "G","SE","FX", "W","SE2")
  list("Sensitivities"=Sensum, "Elasticities"=Elastsum)
	}
	
tempsum1 <- Sensit.Temp.Sum()
tempsum2 <- Sensit.Temp.Sum(senstemp=senstemp2,Tem=10.5,res=resS2)
tempsum3 <- Sensit.Temp.Sum(senstemp=senstemp3,Tem=10.5,res=resS3)
tempsum4 <- Sensit.Temp.Sum(senstemp=senstemp4,Tem=10.5,res=resS4)

#-----------------
#TO OFFSPRING LENGTH x
#-----------------

Sensit.y <- function(x,Tem=10.5,scen="Tempsize", perturb=1e-5, res=resS1){
umat <-  vmat <- matrix(NA, nrow=n,ncol=m)
for (i in 0:(m-1)){
	umat[,i+1]<- res$w[(i*n+1):((i+1)*n)]
	vmat[,i+1]<- res$v[(i*n+1):((i+1)*n)]
	}
wmat <- diag(1,m,m)
dSxy.dy  <- sxymat<- matrix(0,n,m)#x,y
for(j in 1:m) {#y
	dSxy.dy[,j]<- (survival(x,y=y[j]+perturb,Tem=Tem)-survival(x,y=y[j],Tem=Tem))/perturb
	sxymat[,j]<-survival(x,y=y[j],Tem=Tem)
  }
dGxy.dy  <- gxymat<- array(NA,dim=c(n,n,m))#x,x',y
for(j in 1:m) {#y
	for(i in 1:n){#x
		gxymat[i,,j] <- Growth.x(x=x[i],xprime=x,y=y[j],Tem=Tem)#x,x',y
		dGxy.dy[i,,j] <- (Growth.x(x=x[i],xprime=x,y=y[j]+ perturb,Tem=Tem)-Growth.x(x=x[i],xprime=x,y=y[j],Tem=Tem))/perturb
		}
  }
MatG<- MatS<- array(NA,dim=c(n,n,m,m))#x,x',y,y'
for(i in 1:m){#y
	for(j in 1:m){#y'
		for (l in 1:n){#x'
			MatS[,l,i,j] <- vmat[l,j]*umat[,i]* dSxy.dy[,i]* gxymat[,l,i]*wmat[i,j]*dx 
			MatG[,l,i,j] <- vmat[l,j]*umat[,i]* sxymat[,i]* dGxy.dy[,l,i]*wmat[i,j]*dx 
			}
	}
}

list("S"=MatS,"G"= MatG)
}	


sensy1 <- Sensit.y(x)
sensy2 <- Sensit.y(x,scen="Size",res=resS2)
sensy3 <- Sensit.y(x,scen="Temp",res=resS3)
sensy4 <- Sensit.y(x,scen="None",res=resS4)


Elast.y <- function(x,Tem=10.5,scen="Tempsize", perturb=1e-5, res=resS1){
lam <- res$lam 
umat <-  vmat <- matrix(NA, nrow=n,ncol=m)
for (i in 0:(m-1)){
	umat[,i+1]<- res$w[(i*n+1):((i+1)*n)]
	vmat[,i+1]<- res$v[(i*n+1):((i+1)*n)]
	}
uy <- apply(umat,2,sum)
ux <- apply(umat,1,sum)
vx<-apply(vmat,1,sum)
vy<-apply(vmat,2,sum)
wmat <- diag(1,m,m)
dSxy.dy  <- sxymat<- matrix(0,n,m)#x,y
for(j in 1:m) {#y
	dSxy.dy[,j]<- (survival(x,y=y[j]+perturb,Tem=Tem)-survival(x,y=y[j],Tem=Tem))/perturb
	sxymat[,j]<-survival(x,y=y[j],Tem=Tem)
  }
dGxy.dy  <- gxymat<- array(NA,dim=c(n,n,m))#x,x',y
for(j in 1:m) {#y
	for(i in 1:n){#x
		gxymat[i,,j] <- Growth.x(x=x[i],xprime=x,y=y[j],Tem=Tem)#x,x',y
		dGxy.dy[i,,j] <- (Growth.x(x=x[i],xprime=x,y=y[j]+ perturb,Tem=Tem)-Growth.x(x=x[i],xprime=x,y=y[j],Tem=Tem))/perturb
		}
  }
MatG<- MatS<- array(NA,dim=c(n,n,m,m))#x,x',y,y'
for(i in 1:m){#y
	for(j in 1:m){#y'
		for (l in 1:n){#x'
			MatS[,l,i,j] <- vmat[l,j]*umat[,i]* dSxy.dy[,i]* gxymat[,l,i]*wmat[i,j]*y[i]/lam*dx 
			MatG[,l,i,j] <- vmat[l,j]*umat[,i]* sxymat[,i]* dGxy.dy[,l,i]*wmat[i,j]*y[i]/lam*dx
			}
	}
}

list("S"=MatS,"G"= MatG)
}	


elasty1 <- Elast.y(x)
elasty2 <- Elast.y(x,scen="Size",res=resS2)
elasty3 <- Elast.y(x,scen="Temp",res=resS3)
elasty4 <- Elast.y(x,scen="None",res=resS4)

#SUMS
Sensit.y.Sum <- function(sensy= sensy1, elasty=elasty1,res=resS1){
	Sensum<- Elastsum<-  rep(NA,2)
	for(i in 1:2){
		Sensum[i]<- sum(apply(sensy[[i]],c(1,2),sum))
		Elastsum[i]<- sum(apply(elasty[[i]],c(1,2),sum))
	}
	names(Sensum) <- names(Elastsum) <- c( "S","G" )
  list("Sensitivities"=Sensum, "Elasticities"=Elastsum)
	}
	
sensysum1 <- Sensit.y.Sum()
sensysum2 <- Sensit.y.Sum(sensy= sensy2, elasty=elasty2,res=resS2)
sensysum3 <- Sensit.y.Sum(sensy= sensy3, elasty=elasty3,res=resS3)
sensysum4 <- Sensit.y.Sum(sensy= sensy4, elasty=elasty4,res=resS4)

	


#==================================================
#Plot sensitivities, T=10.5 
#==================================================

par(mfrow=c(4,3), mar=c(4,4,2,4),las=0,cex=1,bty="l",xaxs="i",yaxs="i",cex=1)
ylims <-c(-0.003,.008)
xlims <- c(10,120)
plot(x,apply(apply(senstemp1[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="Sensitivity contribution",xlab=" ",ylim= ylims,xlim=xlims,main="To temperature T")
mtext("A",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp1[[2]],c(1,2),sum),1,sum) , lwd=2,col=3)
lines(x,apply(apply(senstemp1[[4]],c(1,2),sum),1,sum) , lwd=2,col=4)
lines(x,apply(apply(senstemp1[[5]],c(1,2),sum),1,sum) , lwd=2,col="Orange")
lines(x,apply(apply(senstemp1[[6]],c(1,2),sum),1,sum) , lwd=2,col="Purple")
legend(10,.008,c("Offspring survival (0.051)", "Offspring length (0.048)","Survival (-0.040)","Growth (0.019)", "Egg weight (-0.009)"),col=c("Purple",4,2,3,"Orange"),lty=1,bty="n",lwd=2,text.col=c("Purple",4,2,3,"Orange"))

plot(x,apply(apply(sensdec1[[2]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To current length x")
mtext("B",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec1[[3]],c(1,2),sum),1,sum) , lwd=2,col=3)
lines(x,apply(apply(sensdec1[[1]],c(1,2),sum),1,sum) , lwd=2,col=1)
lines(x,apply(apply(sensdec1[[5]],c(1,2),sum),1,sum) , lwd=2,col=6)
lines(x,apply(apply(sensdec1[[4]],c(1,2),sum),1,sum) , lwd=2,col="Orange")
legend(30,.008,c("Survival (0.051)","Growth (0.025)",  "Fecundity (0.009)", "Maturity (0.001)","Egg weigth (0.001)"),col=c(2,3,1,6,"Orange"),lty=1,bty="n",lwd=2,text.col=c(2,3,1,6,"Orange"))

plot(x,apply(apply(sensy1[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To offspring length y")
mtext("C",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensy1[[2]],c(1,2),sum),1,sum), lwd=2,col=3)
legend(10,.008,c("Survival (-0.022)","Growth (0.005)"),col=c(2,3),lty=1,bty="n",lwd=2,text.col=c(2,3))
plot(x,apply(apply(senstemp2[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="Sensitivity contribution",xlab=" ",ylim= ylims,xlim=xlims,main="To temperature T")
mtext("D",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp2[[2]],c(1,2),sum),1,sum) , lwd=2,col=3)
lines(x,apply(apply(senstemp2[[4]],c(1,2),sum),1,sum) , lwd=2,col=4)
lines(x,apply(apply(senstemp2[[5]],c(1,2),sum),1,sum) , lwd=2,col="Orange")
legend(10,.008,c( "Egg weight (-0.051)", "Offspring length (0.042)","Survival (-0.042)","Growth (0.019)"),col=c("Orange",4,2,3 ),lty=1,bty="n",lwd=2,text.col=c( "Orange",4,2,3 ))

plot(x,apply(apply(sensdec2[[2]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To current length x")
mtext("E",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec2[[3]],c(1,2),sum),1,sum) , lwd=2,col=3)
lines(x,apply(apply(sensdec2[[1]],c(1,2),sum),1,sum) , lwd=2,col=1)
lines(x,apply(apply(sensdec2[[5]],c(1,2),sum),1,sum) , lwd=2,col=6)
lines(x,apply(apply(sensdec2[[4]],c(1,2),sum),1,sum) , lwd=2,col="Orange")
legend(30,.008,c("Survival (0.046)","Growth (0.025)",  "Fecundity (0.007)","Egg weigth (0.002)", "Maturity (0.000)"),col=c(2,3,1, "Orange",6),lty=1,bty="n",lwd=2,text.col=c(2,3,1, "Orange", 6))

plot(x,apply(apply(sensy2[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To offspring length y")
mtext("F",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensy2[[2]],c(1,2),sum),1,sum), lwd=2,col=3)
legend(10,.008,c("Survival (-0.021)","Growth (0.005)"),col=c(2,3),lty=1,bty="n",lwd=2,text.col=c(2,3))
plot(x,apply(apply(senstemp3[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="Sensitivity contribution",xlab=" ",ylim= ylims,xlim=xlims,main="To temperature T")
mtext("G",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp3[[2]],c(1,2),sum),1,sum) , lwd=2,col=3)
lines(x,apply(apply(senstemp3[[4]],c(1,2),sum),1,sum) , lwd=2,col=4)
lines(x,apply(apply(senstemp3[[6]],c(1,2),sum),1,sum) , lwd=2,col="Purple")
legend(10,.008,c("Offspring survival (0.119)", "Offspring length (0.050)","Survival (-0.039)","Growth (0.020)" ),col=c("Purple",4,2,3 ),lty=1,bty="n",lwd=2,text.col=c("Purple",4,2,3 ))

plot(x,apply(apply(sensdec3[[2]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To current length x")
mtext("H",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec3[[3]],c(1,2),sum),1,sum) , lwd=2,col=3)
lines(x,apply(apply(sensdec3[[1]],c(1,2),sum),1,sum) , lwd=2,col=1)
lines(x,apply(apply(sensdec3[[5]],c(1,2),sum),1,sum) , lwd=2,col=6)

legend(30,.008,c("Survival (0.053)","Growth (0.026)",  "Fecundity (0.010)", "Maturity (0.001)" ),col=c(2,3,1,6 ),lty=1,bty="n",lwd=2,text.col=c(2,3,1,6 ))

plot(x,apply(apply(sensy3[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To offspring length y")
mtext("I",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensy3[[2]],c(1,2),sum),1,sum), lwd=2,col=3)
legend(10,.008,c("Survival (-0.022)","Growth (0.005)"),col=c(2,3),lty=1,bty="n",lwd=2,text.col=c(2,3))

plot(x,apply(apply(senstemp4[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="Sensitivity contribution",xlab=" ",ylim= ylims,xlim=xlims,main="To temperature T")
mtext("J",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp4[[2]],c(1,2),sum),1,sum) , lwd=2,col=3)
lines(x,apply(apply(senstemp4[[4]],c(1,2),sum),1,sum) , lwd=2,col=4)

legend(10,.008,c("Offspring length (0.050)","Survival (-0.039)","Growth (0.020)" ),col=c( 4,2,3 ),lty=1,bty="n",lwd=2,text.col=c(4,2,3 ))

plot(x,apply(apply(sensdec4[[2]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To current length x")
mtext("K",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec4[[3]],c(1,2),sum),1,sum) , lwd=2,col=3)
lines(x,apply(apply(sensdec4[[1]],c(1,2),sum),1,sum) , lwd=2,col=1)
lines(x,apply(apply(sensdec4[[5]],c(1,2),sum),1,sum) , lwd=2,col=6)
legend(30,.008,c("Survival (0.053)","Growth (0.026)",  "Fecundity (0.010)", "Maturity (0.001)" ),col=c(2,3,1,6 ),lty=1,bty="n",lwd=2,text.col=c(2,3,1,6 ))

plot(x,apply(apply(sensy4[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To offspring length y")
mtext("L",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensy4[[2]],c(1,2),sum),1,sum), lwd=2,col=3)
legend(10,.008,c("Survival (-0.022)","Growth (0.005)"),col=c(2,3),lty=1,bty="n",lwd=2,text.col=c(2,3))
 

#==================================================
#Plot sensitivities, T=10.5, S1 only
#==================================================

par(mfrow=c(4,3), mar=c(4,4,2,4),las=0,cex=1,bty="l",xaxs="i",yaxs="i",cex=1)
ylims <-c(-0.003,.007)
xlims <- c(10,120)
plot(x,apply(apply(senstemp1[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="Sensitivity contribution",xlab=" ",ylim= ylims,xlim=xlims,main="To temperature T")
mtext("A",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp1[[2]],c(1,2),sum),1,sum) , lwd=2,col=3)
lines(x,apply(apply(senstemp1[[4]],c(1,2),sum),1,sum) , lwd=2,col=4)
legend(10,.007,c( "Offspring length (0.048)","Survival (-0.040)","Growth (0.019)"),col=c(4,2,3),lty=1,bty="n",lwd=2,text.col=c(4,2,3))

plot(x,apply(apply(sensdec1[[2]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To current length x")
mtext("B",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec1[[3]],c(1,2),sum),1,sum) , lwd=2,col=3)
lines(x,apply(apply(sensdec1[[1]],c(1,2),sum),1,sum) , lwd=2,col=1)
lines(x,apply(apply(sensdec1[[5]],c(1,2),sum),1,sum) , lwd=2,col=6)
legend(30,.007,c("Survival (0.051)","Growth (0.025)",  "Fecundity (0.009)", "Maturity (0.001)"),col=c(2,3,1,6),lty=1,bty="n",lwd=2,text.col=c(2,3,1,6))

plot(x,apply(apply(sensy1[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To offspring length y")
mtext("C",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensy1[[2]],c(1,2),sum),1,sum), lwd=2,col=3)
legend(10,.007,c("Survival (-0.022)","Growth (0.005)"),col=c(2,3),lty=1,bty="n",lwd=2,text.col=c(2,3))

plot(x,apply(apply(senstemp1[[5]],c(1,2),sum),1,sum) ,ylim= ylims,xlim=xlims, lty=1,lwd=2,col=1,type="l", ylab="Sensitivity contribution",xlab="Female length x (cm)",main="To T via egg weight")
mtext("D",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp2[[5]],c(1,2),sum),1,sum) , lwd=2,col=1,lty=2)
legend(10,.007,c("Scenario 1 Interaction (-0.009)","Scenario 2 Size (-0.051)"),col=1,lty=c(1,2),bty="n",lwd=2)

plot(x,apply(apply(senstemp1[[6]],c(1,2),sum),1,sum),ylim= ylims,xlim=xlims,lwd=2,col=1,type="l", ylab="",xlab="Female length x (cm)",main="To T via offspring survival")
mtext("E",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp3[[6]],c(1,2),sum),1,sum) , lwd=2,col=1,lty=3)
legend(10,.007,c("Scenario 1 Interaction (0.051)",  "Scenario 3 Temp (0.119)"),col=1,lty=c(1,3),bty="n",lwd=2)

plot(x,apply(apply(sensdec1[[4]],c(1,2),sum),1,sum) ,ylim= ylims,xlim=xlims, lwd=2,col=1,type="l", ylab="",xlab="Female length x (cm)",main="To x via egg weight")
mtext("F",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec2[[4]],c(1,2),sum),1,sum)  , lwd=2,col=1,lty=2)
legend(10,.007,c("Scenario 1 Interaction (0.001)","Scenario 2 Size (0.002)"),col=1,lty=c(1,2),bty="n",lwd=2)
 

#==================================================
#Plot elasticities, T=10.5, S1 only (A-C)
#==================================================
par(mfrow=c(2,3), mar=c(4,4,2,4),las=0,cex=1,bty="l",xaxs="i",yaxs="i",cex=1)
ylims <- c(-0.05,.15)
xlims <- c(10,120)
plot(x,apply(apply(senstemp1[[1]],c(1,2),sum),1,sum)*10.5/resS1$lam,xlim=xlims, lwd=2,col=2,type="l", ylab="Elasticity contribution",xlab=" ",ylim= ylims,main="To temperature T")
mtext("A",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp1[[2]],c(1,2),sum),1,sum)*10.5/resS1$lam, lwd=2,col=3)
lines(x,apply(apply(senstemp1[[4]],c(1,2),sum),1,sum)*10.5/resS1$lam, lwd=2,col=4)
legend(10,ylims[2],c( "Offspring length (0.483)","Survival (-0.399)","Growth (0.195)"),col=c(4,2,3),lty=1,bty="n",lwd=2,text.col=c(4,2,3))

plot(x,apply(apply(sensdec1[[2]],c(1,2),sum),1,sum)*x/resS1$lam,xlim=xlims,lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,main="To current length x")
mtext("B",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec1[[3]],c(1,2),sum),1,sum)*x/resS1$lam, lwd=2,col=3)
lines(x,apply(apply(sensdec1[[1]],c(1,2),sum),1,sum)*x/resS1$lam, lwd=2,col=1)
lines(x,apply(apply(sensdec1[[5]],c(1,2),sum),1,sum)*x/resS1$lam, lwd=2,col=6)
legend(30,ylims[2],c("Survival (1.350)","Growth (1.062)",  "Fecundity (0.616)", "Maturity (0.028)"),col=c(2,3,1,6),lty=1,bty="n",lwd=2,text.col=c(2,3,1,6))

plot(x,apply(apply(elasty1[[1]],c(1,2),sum),1,sum),xlim=xlims , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,main="To offspring length y")
mtext("C",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(elasty1[[2]],c(1,2),sum),1,sum), lwd=2,col=3)
legend(10,ylims[2],c("Survival (-0.537)","Growth (0.127)"),col=c(2,3),lty=1,bty="n",lwd=2,text.col=c(2,3))

plot(x,apply(apply(senstemp1[[5]],c(1,2),sum),1,sum)*10.5/resS1$lam,xlim=xlims, ylim= ylims, lty=1,lwd=2,col=1,type="l", ylab="Elasticity contribution",xlab="Female length x (cm)",main="To T via egg weight")
mtext("D",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp2[[5]],c(1,2),sum),1,sum)*10.5/resS2$lam, lwd=2,col=1,lty=2)
legend(10,ylims[2],c("Scenario 2 Size (-0.523)","Scenario 1 Interaction (-0.086)"),col=1,lty=c(2,1),bty="n",lwd=2)

plot(x,apply(apply(senstemp1[[6]],c(1,2),sum),1,sum)*10.5/resS1$lam, ylim= ylims,lwd=2,col=1,type="l", ylab=" ",xlab="Female length x (cm)",main="To T via offspring survival",xlim=xlims)
mtext("E",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp3[[6]],c(1,2),sum),1,sum)*10.5/resS3$lam, lwd=2,col=1,lty=3)

legend(10,ylims[2],c("Scenario 3 Temp (1.187)","Scenario 1 Interaction (0.515)"),col=1,lty=c(3,1),bty="n",lwd=2)

plot(x,apply(apply(sensdec1[[4]],c(1,2),sum),1,sum)*x/resS1$lam , ylim=ylims,xlim=xlims , lwd=2,col=1,type="l", ylab=" ",xlab="Female length x (cm)",main="To x via egg weigth")
mtext("F",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec2[[4]],c(1,2),sum),1,sum)*x/resS2$lam , lwd=2,col=1,lty=2)
legend(10,ylims[2],c("Scenario 2 Size (0.147)","Scenario 1 Interaction (0.033)"),col=1,lty=c(2,1),bty="n",lwd=2)
 

#==================================================
#Plot elasticities, T=10.5 
#==================================================

par(mfrow=c(4,3), mar=c(4,4,2,4),las=0,cex=1,bty="l",xaxs="i",yaxs="i",cex=1)
ylims <- c(-0.05,.15)
xlims <- c(10,120)
plot(x,apply(apply(senstemp1[[1]],c(1,2),sum),1,sum)*10.5/resS1$lam , lwd=2,col=2,type="l", ylab="Elasticity contribution",xlab=" ",ylim= ylims,xlim=xlims,main="To temperature T")
mtext("A",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp1[[2]],c(1,2),sum),1,sum)*10.5/resS1$lam , lwd=2,col=3)
lines(x,apply(apply(senstemp1[[4]],c(1,2),sum),1,sum)*10.5/resS1$lam , lwd=2,col=4)
lines(x,apply(apply(senstemp1[[5]],c(1,2),sum),1,sum)*10.5/resS1$lam , lwd=2,col="Orange")
lines(x,apply(apply(senstemp1[[6]],c(1,2),sum),1,sum)*10.5/resS1$lam , lwd=2,col="Purple")
legend(10,ylims[2],c("Offspring survival (0.515)", "Offspring length (0.483)","Survival (-0.399)","Growth (0.195)", "Egg weight (-0.086)"),col=c("Purple",4,2,3,"Orange"),lty=1,bty="n",lwd=2,text.col=c("Purple",4,2,3,"Orange"))

plot(x,apply(apply(sensdec1[[2]],c(1,2),sum),1,sum)*x/resS1$lam  , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To current length x")
mtext("B",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec1[[3]],c(1,2),sum),1,sum)*x/resS1$lam  , lwd=2,col=3)
lines(x,apply(apply(sensdec1[[1]],c(1,2),sum),1,sum)*x/resS1$lam  , lwd=2,col=1)
lines(x,apply(apply(sensdec1[[5]],c(1,2),sum),1,sum)*x/resS1$lam  , lwd=2,col=6)
lines(x,apply(apply(sensdec1[[4]],c(1,2),sum),1,sum)*x/resS1$lam  , lwd=2,col="Orange")
legend(30,ylims[2],c("Survival (1.350)","Growth (1.062)",  "Fecundity (0.616)", "Egg weigth (0.033)", "Maturity (0.028)"),col=c(2,3,1, "Orange",6),lty=1,bty="n",lwd=2,text.col=c(2,3,1, "Orange",6))

plot(x,apply(apply(elasty1[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims,main="To offspring length y")
mtext("C",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(elasty1[[2]],c(1,2),sum),1,sum), lwd=2,col=3)
legend(10,ylims[2],c("Survival (-0.537)","Growth (0.127)"),col=c(2,3),lty=1,bty="n",lwd=2,text.col=c(2,3))

plot(x,apply(apply(senstemp2[[1]],c(1,2),sum),1,sum)*10.5/resS2$lam , lwd=2,col=2,type="l", ylab="Elasticity contribution",xlab=" ",ylim= ylims,xlim=xlims )
mtext("D",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp2[[2]],c(1,2),sum),1,sum)*10.5/resS2$lam , lwd=2,col=3)
lines(x,apply(apply(senstemp2[[4]],c(1,2),sum),1,sum)*10.5/resS2$lam , lwd=2,col=4)
lines(x,apply(apply(senstemp2[[5]],c(1,2),sum),1,sum)*10.5/resS2$lam , lwd=2,col="Orange")
legend(10,ylims[2],c( "Egg weight (-0.523)", "Offspring length (0.426)","Survival (-0.424)","Growth (0.193)"),col=c("Orange",4,2,3 ),lty=1,bty="n",lwd=2,text.col=c( "Orange",4,2,3 ))

plot(x,apply(apply(sensdec2[[2]],c(1,2),sum),1,sum)*x/resS2$lam  , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims )
mtext("E",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec2[[3]],c(1,2),sum),1,sum)*x/resS2$lam  , lwd=2,col=3)
lines(x,apply(apply(sensdec2[[1]],c(1,2),sum),1,sum)*x/resS2$lam  , lwd=2,col=1)
lines(x,apply(apply(sensdec2[[5]],c(1,2),sum),1,sum)*x/resS2$lam  , lwd=2,col=6)
lines(x,apply(apply(sensdec2[[4]],c(1,2),sum),1,sum)*x/resS2$lam  , lwd=2,col="Orange")
legend(30,ylims[2],c("Survival (1.245)","Growth (1.122)",  "Fecundity (0.540)","Egg weigth (0.147)", "Maturity (0.010)"),col=c(2,3,1, "Orange",6),lty=1,bty="n",lwd=2,text.col=c(2,3,1, "Orange", 6))

plot(x,apply(apply(elasty2[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims )
mtext("F",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(elasty2[[2]],c(1,2),sum),1,sum), lwd=2,col=3)
legend(10,ylims[2],c("Survival (-0.530)","Growth (0.126)"),col=c(2,3),lty=1,bty="n",lwd=2,text.col=c(2,3))

plot(x,apply(apply(senstemp3[[1]],c(1,2),sum),1,sum)*10.5/resS3$lam , lwd=2,col=2,type="l", ylab="Elasticity contribution",xlab=" ",ylim= ylims,xlim=xlims )
mtext("G",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp3[[2]],c(1,2),sum),1,sum)*10.5/resS3$lam , lwd=2,col=3)
lines(x,apply(apply(senstemp3[[4]],c(1,2),sum),1,sum)*10.5/resS3$lam , lwd=2,col=4)
lines(x,apply(apply(senstemp3[[6]],c(1,2),sum),1,sum)*10.5/resS3$lam , lwd=2,col="Purple")
legend(10,ylims[2],c("Offspring survival (1.187)", "Offspring length (0.499)","Survival (-0.392)","Growth (0.195)" ),col=c("Purple",4,2,3 ),lty=1,bty="n",lwd=2,text.col=c("Purple",4,2,3 ))

plot(x,apply(apply(sensdec3[[2]],c(1,2),sum),1,sum)*x/resS3$lam  , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims )
mtext("H",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec3[[3]],c(1,2),sum),1,sum)*x/resS3$lam , lwd=2,col=3)
lines(x,apply(apply(sensdec3[[1]],c(1,2),sum),1,sum)*x/resS3$lam , lwd=2,col=1)
lines(x,apply(apply(sensdec3[[5]],c(1,2),sum),1,sum)*x/resS3$lam , lwd=2,col=6)

legend(30,ylims[2],c("Survival (1.379)","Growth (1.045)",  "Fecundity (0.638)", "Maturity (0.034)" ),col=c(2,3,1,6 ),lty=1,bty="n",lwd=2,text.col=c(2,3,1,6 ))

plot(x,apply(apply(elasty3[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab=" ",ylim= ylims,xlim=xlims )
mtext("I",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(elasty3[[2]],c(1,2),sum),1,sum), lwd=2,col=3)
legend(10,ylims[2],c("Survival (-0.539)","Growth (0.127)"),col=c(2,3),lty=1,bty="n",lwd=2,text.col=c(2,3))

plot(x,apply(apply(senstemp4[[1]],c(1,2),sum),1,sum)*10.5/resS4$lam , lwd=2,col=2,type="l", ylab="Elasticity contribution",xlab="Current length x (cm)",ylim= ylims,xlim=xlims )
mtext("J",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(senstemp4[[2]],c(1,2),sum),1,sum)*10.5/resS4$lam , lwd=2,col=3)
lines(x,apply(apply(senstemp4[[4]],c(1,2),sum),1,sum)*10.5/resS4$lam , lwd=2,col=4)

legend(10,ylims[2],c("Offspring length (0.499)","Survival (-0.392)","Growth (0.195)" ),col=c( 4,2,3 ),lty=1,bty="n",lwd=2,text.col=c(4,2,3 ))

plot(x,apply(apply(sensdec4[[2]],c(1,2),sum),1,sum)*x/resS4$lam , lwd=2,col=2,type="l", ylab="",xlab="Current length x (cm)",ylim= ylims,xlim=xlims)
mtext("K",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(sensdec4[[3]],c(1,2),sum),1,sum)*x/resS4$lam , lwd=2,col=3)
lines(x,apply(apply(sensdec4[[1]],c(1,2),sum),1,sum)*x/resS4$lam , lwd=2,col=1)
lines(x,apply(apply(sensdec4[[5]],c(1,2),sum),1,sum)*x/resS4$lam , lwd=2,col=6)
legend(30,ylims[2],c("Survival (1.379)","Growth (1.045)",  "Fecundity (0.638)", "Maturity (0.034)" ),col=c(2,3,1,6 ),lty=1,bty="n",lwd=2,text.col=c(2,3,1,6 ))

plot(x,apply(apply(elasty4[[1]],c(1,2),sum),1,sum) , lwd=2,col=2,type="l", ylab="",xlab="Current length x (cm)",ylim= ylims,xlim=xlims,main=" ")
mtext("L",3,at=13,cex=1.5)
abline(h=0,lty=3,lwd=2,col=8)
lines(x,apply(apply(elasty4[[2]],c(1,2),sum),1,sum), lwd=2,col=3)
legend(10,ylims[2],c("Survival (-0.539)","Growth (0.127)"),col=c(2,3),lty=1,bty="n",lwd=2,text.col=c(2,3))
 
  