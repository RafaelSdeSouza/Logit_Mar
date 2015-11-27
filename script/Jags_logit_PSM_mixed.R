# JAGS Code 
#  Required libraries
library(rjags);library(ggmcmc);library(ggplot2);library(ggthemes);library(pander);library(Cairo);library(MASS);library(parallel) 
library(scales);library(plyr);require(gdata);require(runjags);require(gdata);require(caret);require(pROC);require(plyr);library(grid)
cl       <- makeCluster(3) # For parallel computing

source("facet_wrap_labeller.R")
# Read and format data
data     <- read.table("..//data/matched.txt",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","logM200_L","RprojLW_Rvir","zoo")]


X          <- model.matrix( ~  logM200_L + RprojLW_Rvir, data = data_cut) # Predictors
K          <- ncol(X)                   # Number of Predictors including the intercept 
y          <- data_cut$bpt # Response variable (0/1)
n          <- length(y)                 # Sample size
gal        <- as.numeric(data_cut$zoo)

# Grid of values for prediction 

Mx <- seq(2*min(X[,2]),2*max(X[,2]),length.out=300)
Rx <- seq(2*min(X[,3]),2*max(X[,3]),length.out=300)


jags.data  <- list(Y= y,N = n,X=X,b0 = rep(0,K),B0=diag(1e-4,K),Mx=Mx,
                   Rx=Rx,gal=gal)
# b0 and B0 are the mean vector and the inverse of the covariance matrix of prior distribution of the regression coefficients

### Prior Specification and Likelihood function
model<-"model{
 
#0. Hyperprior
tau ~ dgamma(0.001,0.001)
mu  ~ dnorm(0,1e-3)

#1. Priors

#a.Normal 
for(j in 1:2){
for(k in 1:3){
beta[k,j]~dnorm(mu,tau)                                    
}
}

#2. Likelihood
for(i in 1:N){

Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
eta[i] <- beta[1,gal[i]]*X[i,1]+beta[2,gal[i]]*X[i,2]+beta[3,gal[i]]*X[i,3]

}
#3.Probability for each variable 
# 
#for(l in 1:300){

#a. Ellipticals

logit(pMxE[l])<-beta[1,1]+beta[2,1]*Mx[l]
logit(pRxE[l])<-beta[1,1]+beta[3,1]*Rx[l]
logit(pxE[l])<-beta[1,1]+beta[2,1]*Mx[l]+beta[3,1]*Rx[l]

#b. Spirals 

logit(pMxS[l])<-beta[1,2]+beta[2,2]*Mx[l]
logit(pRxS[l])<-beta[1,2]+beta[3,2]*Rx[l]
logit(pxS[l])<-beta[1,2]+beta[2,2]*Mx[l]+beta[3,2]*Rx[l]
             }
             }"
params <- c("beta","pi","pMx","pRx","px")        # Monitor these parameters.
inits0  <- function () {list(beta.0 = rnorm(2,0, 0.01),beta.1 = rnorm(2,0, 0.01),beta.2 = rnorm(2,0, 0.01))} # A function to generat initial values for mcmc
inits1=inits0();inits2=inits0();inits3=inits0()         # Generate initial values for three chains

# Run mcmc
bin    = 10^4   # burn-in samples
ad     = 10^4   # Number of adaptive samples
s      = 3*10^4 # Number of samples for each chain
nc     = 3      # Number of mcmc
th     = 10     # Thinning value
jags.logit  <- run.jags(method="rjparallel",data = jags.data,inits = list(inits1,inits2,inits3),model=model,
                       n.chains = nc,adapt=ad,monitor=c(params),burnin=bin,thin=th,sample=s,summarise=FALSE,plots=FALSE)




#-----------##-----------#-----------##-----------
pi_Rx<-summary(as.mcmc.list(jags.logit, vars="pRx"))
pi_Rx<-pi_Rx$quantiles
gRx<-data.frame(x=Rx,mean=pi_Rx[,3],lwr1=pi_Rx[,2],lwr2=pi_Rx[,1],upr1=pi_Rx[,4],upr2=pi_Rx[,5])

#-----------##
pi_Rx<-summary(as.mcmc.list(jags.logit, vars="pRx"))
pi_Rx<-pi_Rx$quantiles
gRx<-data.frame(x=Rx,mean=pi_Rx[,3],lwr1=pi_Rx[,2],lwr2=pi_Rx[,1],upr1=pi_Rx[,4],upr2=pi_Rx[,5])
#-----------##-----------#-----------##-----------




#-----------##-----------#-----------##-----------
pi_Mx<-summary(as.mcmc.list(jags.logit, vars="pMx"))
pi_Mx<-pi_Mx$quantiles
gMx<-data.frame(x=Mx,mean=pi_Mx[,3],lwr1=pi_Mx[,2],lwr2=pi_Mx[,1],upr1=pi_Mx[,4],upr2=pi_Mx[,5])

#-----------##
pi_Mx<-summary(as.mcmc.list(jags.logit, vars="pMx"))
pi_Mx<-pi_Mx$quantiles
gMx<-data.frame(x=Mx,mean=pi_Mx[,3],lwr1=pi_Mx[,2],lwr2=pi_Mx[,1],upr1=pi_Mx[,4],upr2=pi_Mx[,5])
#-----------##-----------#-----------##-----------



write.table(gRx,"gRx_S.dat",row.names = F)
write.table(gMx,"gMx_S.dat",row.names = F)


# Plots with ggmcmc
jagssamples <- as.mcmc.list(jags.logit)

G1<-ggs(jagssamples,family="beta")
gplot<-G1[,3:4]
gplot$gal<-rep(c("E","S"),each=nrow(gplot)/2)
write.table(gplot,"gplot.dat",row.names = F)


ggs_density(G1)


# Trace plots and diagnostic analysis to investigate convregence
#plot(beta) 
#gelman.diag(beta)
#gelman.plot(beta)
#geweke.diag(beta)
#geweke.plot(beta)
#autocorr.plot(beta)

# Summary Statistics for parameters of interest
#print(summary(window(beta)))
#print(summary(window(pi)))




# Plot all probabilities 
#pi_AGN<-summary(as.mcmc.list(jags.logit, vars="pi"))
#pi_AGN<-pi_AGN$quantiles

#gdata<-data.frame(x=data_n2$sfr_tot_p50,mean=pi_AGN[,3],lwr1=pi_AGN[,2],lwr2=pi_AGN[,1],upr1=pi_AGN[,4],upr2=pi_AGN[,5])




