# JAGS Code with Adaptive Shrinkage
#  Required libraries
library(rjags)
library(ggmcmc)
library(ggplot2)
library(ggthemes)
library(pander)
library(Cairo)
library(MASS)
library(scales)
library(plyr)
require(gdata)
require(runjags)
require(gdata)
require(caret)
require(pROC)
require(plyr)


# Read and format data
data<-read.csv("sample_agn.csv",header=TRUE,na.strings="")
data2<-na.omit(data)
data2<-data2[data2$logMstar>0,]
data2<-data2[which(data2$vlos_sigma!=Inf),]

#data2$bpt<-as.factor(data2$bpt)
data2$bpt <- revalue(data2$bpt,c("SF"="0","Composite"="0",
                                 "LINER"="1","Seyfert/LINER"="1",
                                 "Seyfert"="1","BLANK"="0"))


# Prepare data for JAGS
data2$logMstar<-(data2$logMstar-mean(data2$logMstar))/sd(data2$logMstar)
data2$logMhalo<-(data2$logMhalo-mean(data2$logMhalo))/sd(data2$logMhalo)
data2$vlos_sigma<-(data2$vlos_sigma-mean(data2$vlos_sigma))/sd(data2$vlos_sigma)
data2$r_rvir<-(data2$r_rvir-mean(data2$r_rvir))/sd(data2$r_rvir)


X<-model.matrix(~logMstar+logMhalo+vlos_sigma+r_rvir,data=data2)
K<-ncol(X)

jags.data <- list(Y= data2$bpt,
                  N = nrow(data2),
                  X=X,
                  b0 = rep(0,K),
                  B0=diag(1e-4,K),
                  Npred = K
)
model<-"model{
#1. Priors 
#beta~dmnorm(b0[],B0[,]) # Normal Priors
# Jefreys priors for sparseness 
for(j in 1:Npred)   {
lnTau[j] ~ dunif(-50, 50)   
TauM[j] <- exp(lnTau[j])
beta[j] ~ dnorm(0, TauM[j]) 
Ind[j] <- step(abs(beta[j]) - 0.05)
}

#2. Likelihood

for(i in 1:N){
Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
eta[i] <- inprod(beta[], X[i,])


#3. Prediction
NewPred[i]~dbern(pi[i])
}

}"

params <- c("beta","pi","Ind","NewPred")

inits0  <- function () {
  list(beta = rnorm(K, 0, 0.1))}



inits1=inits0()
inits2=inits0()
inits3=inits0()

library(parallel)
cl <- makeCluster(3)
jags.logit <- run.jags(method="rjparallel", 
                       data = jags.data, 
                       inits = list(inits1,inits2,inits3),
                       model=model,
                       n.chains = 3,
                       adapt=1000,
                       monitor=c(params),
                       burnin=5000,
                       sample=15000,
                       summarise=FALSE,
                       plots=FALSE
)

jagssamples <- as.mcmc.list(jags.logit)



