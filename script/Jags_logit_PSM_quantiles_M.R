# JAGS Code 
#  Required libraries
library(rjags);library(ggmcmc);library(ggplot2);library(ggthemes);library(pander);library(Cairo);library(MASS);library(parallel) 
library(scales);library(plyr);require(gdata);require(runjags);require(gdata);require(caret);require(pROC);require(plyr);library(grid)
cl       <- makeCluster(3) # For parallel computing

source("facet_wrap_labeller.R")
# Read and format data
dataE     <- read.table("..//data/matched_E.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S.txt",header=TRUE,na.strings="")
data<-rbind(dataE,dataS)
data_cut <- data[,c("bpt","logM200_L","RprojLW_Rvir","zoo")]



X          <- model.matrix( ~  logM200_L + RprojLW_Rvir, data = data_cut) # Predictors
K          <- ncol(X)                   # Number of Predictors including the intercept 
y          <- data_cut$bpt # Response variable (0/1)
n          <- length(y)                 # Sample size
gal        <- as.numeric(data_cut$zoo)

# Grid of values for prediction 

Mx <- seq(1.05*min(X[,2]),1.05*max(X[,2]),length.out=300)
Rx <- seq(1.05*min(X[,3]),1.05*max(X[,3]),length.out=300)


#jags.data  <- list(Y= y,N = n,X=X,b0 = rep(0,K),B0=diag(1e-4,K),Mx=Mx,
#                   Rx=Rx,gal=gal)

jags.data  <- list(Y= y,N = n,X=X,Mx=Mx,Rx=Rx,gal=gal)
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
for(l in 1:300){

#a. Ellipticals

#logit(pMxE[l])<-beta[1,1]+beta[2,1]*Mx[l]

logit(pRxEminus[l])<-beta[1,1]+beta[3,1]*Rx[l]+ beta[2,1]*(mean(X[,2])-sd(X[,2])) # at mean -1sigma
logit(pRxE[l])<-beta[1,1]+beta[3,1]*Rx[l]+ beta[2,1]*mean(X[,2]) # at mean Mass
logit(pRxEplus[l])<-beta[1,1]+beta[3,1]*Rx[l]+ beta[2,1]*(mean(X[,2])+sd(X[,2])) # at mean +1sigma

#logit(pxE[l])<-beta[1,1]+beta[2,1]*Mx[l]+beta[3,1]*Rx[l]

#b. Spirals 

#logit(pMxS[l])<-beta[1,2]+beta[2,2]*Mx[l]

logit(pRxSminus[l])<-beta[1,2]+beta[3,2]*Rx[l]+beta[2,2]*(mean(X[,2])-sd(X[,2]))
logit(pRxS[l])<-beta[1,2]+beta[3,2]*Rx[l]+ beta[2,2]*mean(X[,2])
logit(pRxSplus[l])<-beta[1,2]+beta[3,2]*Rx[l]+beta[2,2]*(mean(X[,2])+sd(X[,2]))

#logit(pxS[l])<-beta[1,2]+beta[2,2]*Mx[l]+beta[3,2]*Rx[l]
             }
             }"
params <- c("beta","pi","pRxS","pRxSminus","pRxSplus","pRxE","pRxEminus","pRxEplus")        # Monitor these parameters.
inits0  <- function () {list(beta = matrix(rnorm(6,0, 0.01),ncol=2))} # A function to generat initial values for mcmc
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
pi_RxE<-summary(as.mcmc.list(jags.logit, vars="pRxE"))
pi_RxE<-pi_RxE$quantiles
gRxE<-data.frame(x=Rx,mean=pi_RxE[,3],lwr1=pi_RxE[,2],lwr2=pi_RxE[,1],upr1=pi_RxE[,4],upr2=pi_RxE[,5])


#-----------##
pi_RxE_minus<-summary(as.mcmc.list(jags.logit, vars="pRxEminus"))
pi_RxE_minus<-pi_RxE_minus$quantiles
gRxE_minus<-data.frame(x=Rx,mean=pi_RxE_minus[,3],lwr1=pi_RxE_minus[,2],lwr2=pi_RxE_minus[,1],upr1=pi_RxE_minus[,4],upr2=pi_RxE_minus[,5])

#-----------##
pi_RxE_plus<-summary(as.mcmc.list(jags.logit, vars="pRxEplus"))
pi_RxE_plus<-pi_RxE_plus$quantiles
gRxE_plus<-data.frame(x=Rx,mean=pi_RxE_plus[,3],lwr1=pi_RxE_plus[,2],lwr2=pi_RxE_plus[,1],upr1=pi_RxE_plus[,4],upr2=pi_RxE_plus[,5])


#-----------##
write.table(gRxE[1:300,],"..//data/gRxE.dat",row.names = F)
write.table(gRxE_minus,"..//data/gRxE_minus.dat",row.names = F)
write.table(gRxE_plus,"..//data/gRxE_plus.dat",row.names = F)
#-----------##-----------#-----------##-----------


#-----------###-----------##-----------#-----------##-----------
pi_RxS<-summary(as.mcmc.list(jags.logit, vars="pRxS"))
pi_RxS<-pi_RxS$quantiles
gRxS<-data.frame(x=Rx,mean=pi_RxS[,3],lwr1=pi_RxS[,2],lwr2=pi_RxS[,1],upr1=pi_RxS[,4],upr2=pi_RxS[,5])



#-----------##
pi_RxS_minus<-summary(as.mcmc.list(jags.logit, vars="pRxSminus"))
pi_RxS_minus<-pi_RxS_minus$quantiles
gRxS_minus<-data.frame(x=Rx,mean=pi_RxS_minus[,3],lwr1=pi_RxS_minus[,2],lwr2=pi_RxS_minus[,1],upr1=pi_RxS_minus[,4],upr2=pi_RxS_minus[,5])


#-----------##
pi_RxS_plus<-summary(as.mcmc.list(jags.logit, vars="pRxSplus"))
pi_RxS_plus<-pi_RxS_plus$quantiles
gRxS_plus<-data.frame(x=Rx,mean=pi_RxS_plus[,3],lwr1=pi_RxS_plus[,2],lwr2=pi_RxS_plus[,1],upr1=pi_RxS_plus[,4],upr2=pi_RxS_plus[,5])

#-----------##-----------#-----------##-----------


#-----------##
write.table(gRxS[1:300,],"..//data/gRxS.dat",row.names = F)
write.table(gRxS_minus,"..//data/gRxS_minus.dat",row.names = F)
write.table(gRxS_plus,"..//data/gRxS_plus.dat",row.names = F)
#-----------##-----------#-----------##-----------













#-----------##-----------#-----------##-----------
pi_MxE<-summary(as.mcmc.list(jags.logit, vars="pMxE"))
pi_MxE<-pi_MxE$quantiles
gMxE<-data.frame(x=Mx,mean=pi_MxE[,3],lwr1=pi_MxE[,2],lwr2=pi_MxE[,1],upr1=pi_MxE[,4],upr2=pi_MxE[,5])

#-----------##
pi_MxS<-summary(as.mcmc.list(jags.logit, vars="pMxS"))
pi_MxS<-pi_MxS$quantiles
gMxS<-data.frame(x=Mx,mean=pi_MxS[,3],lwr1=pi_MxS[,2],lwr2=pi_MxS[,1],upr1=pi_MxS[,4],upr2=pi_MxS[,5])
#-----------##-----------#-----------##-----------



write.table(gRxS,"..//data/gRx_S.dat",row.names = F)
write.table(gRxE,"..//data/gRx_E.dat",row.names = F)
write.table(gMxS,"..//data/gMx_S.dat",row.names = F)
write.table(gMxE,"..//data/gMx_E.dat",row.names = F)

# Plots with ggmcmc
jagssamples <- as.mcmc.list(jags.logit)

G1<-ggs(jagssamples,family="beta")
gplot<-G1[,3:4]
gplot$gal<-rep(c("E","S"),each=nrow(gplot)/2)
write.table(gplot,"gplot.dat",row.names = F)




#ggs_density(G1)


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




