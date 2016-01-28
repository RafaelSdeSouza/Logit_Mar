# JAGS Code 
#  Required libraries
library(rjags);library(ggmcmc);library(ggplot2);library(ggthemes);library(pander);library(Cairo);library(MASS);library(parallel) 
library(scales);library(plyr);require(gdata);require(runjags);require(gdata);require(caret);require(pROC);require(plyr);library(grid)
cl       <- makeCluster(3) # For parallel computing

source("facet_wrap_labeller.R")
# Read and format data
dataE     <- read.table("..//data/matched_E_original.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S_original.txt",header=TRUE,na.strings="")
data<-rbind(dataE,dataS)
data_cut <- data[,c("bpt","logM200_L","RprojLW_Rvir","zoo","groupID")]



X          <- model.matrix( ~ RprojLW_Rvir, data = data_cut) # Predictors
K          <- ncol(X)                   # Number of Predictors including the intercept 
y          <- data_cut$bpt # Response variable (0/1)
n          <- length(y)                 # Sample size
gal        <- as.numeric(data_cut$zoo)
GroupID    <- match(data_cut$groupID,unique(data_cut$groupID))
NGroup     <- length(unique(data_cut$groupID))


# Grid of values for prediction 

#Mx <- seq(1.05*min(X[,2]),1.05*max(X[,2]),length.out=300)
Rx <- seq(0.95*min(X[,2]),1.05*max(X[,2]),length.out=300)


#jags.data  <- list(Y= y,N = n,X=X,b0 = rep(0,K),B0=diag(1e-4,K),Mx=Mx,
#                   Rx=Rx,gal=gal)

jags.data  <- list(Y= y,N = n,X=X,Rx=Rx,gal=gal,GroupID=GroupID,NGroup=NGroup)
# b0 and B0 are the mean vector and the inverse of the covariance matrix of prior distribution of the regression coefficients

### Prior Specification and Likelihood function
model<-"model{
 
#0. Hyperprior
tau ~ dgamma(0.001,0.001)
mu  ~ dnorm(0,1e-3)

#tau.R<-pow(sdBeta,-1)
#sdBeta ~ dgamma(0.01,0.01)
#1. Priors

#a.Normal 
for(j in 1:2){
for(k in 1:2){
beta[k,j]~dnorm(mu,tau)                                    
}
}

#.b Random Intercept for Halo ID
#for (t in 1:NGroup){
#ranef[t]~ddexp(0,tau.R)
#}


#2. Likelihood
for(i in 1:N){

Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
#eta[i] <- beta[1,gal[i]]*X[i,1]+beta[2,gal[i]]*X[i,2]+ranef[GroupID[i]]
eta[i] <- beta[1,gal[i]]*X[i,1]+beta[2,gal[i]]*X[i,2]

}
#3.Probability for each variable 
# 
for(l in 1:300){

#a. Ellipticals


logit(pRxE[l])<-beta[1,1]+beta[2,1]*Rx[l]
logit(pxE[l])<-beta[1,1]+beta[2,1]*Rx[l]

#b. Spirals 


logit(pRxS[l])<-beta[1,2]+beta[1,2]*Rx[l]
logit(pxS[l])<-beta[1,2]+beta[2,2]*Rx[l]
             }
             }"
params <- c("beta","pi","pRxS","pxS","pRxE","pxE")        # Monitor these parameters.
inits0  <- function () {list(beta = matrix(rnorm(4,0, 0.01),ncol=2))} # A function to generat initial values for mcmc
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
pi_RxS<-summary(as.mcmc.list(jags.logit, vars="pRxS"))
pi_RxS<-pi_RxS$quantiles
gRxS<-data.frame(x=Rx,mean=pi_RxS[,3],lwr1=pi_RxS[,2],lwr2=pi_RxS[,1],upr1=pi_RxS[,4],upr2=pi_RxS[,5])
#-----------##-----------#-----------##-----------



write.table(gRxS,"..//data/gRx_S_N.dat",row.names = F)
write.table(gRxE,"..//data/gRx_E_N.dat",row.names = F)


# Plots with ggmcmc
jagssamples <- as.mcmc.list(jags.logit)

G1<-ggs(jagssamples,family="beta")
gplot<-G1[,3:4]
gplot$gal<-rep(c("E","S"),each=nrow(gplot)/2)
write.table(gplot,"gplot.dat",row.names = F)




S.full <- ggs(jagssamples,family=c("ranef"))
library(RColorBrewer)
blues_fun <- colorRampPalette(brewer.pal(9,"Blues")[4:9])
blues=blues_fun(762)


gplot2<-S.full[,3:4]
gplot2$Parameter<-as.factor(gplot2$Parameter)
ggplot(data=gplot2,aes(x=value),group=Parameter)+geom_density(aes(colour=Parameter,alpha=0.4,linetype="dashed"))+
  theme_hc()+ 
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(size=25,vjust=0.75),
        axis.title.x=element_text(size=25,vjust=-0.25),axis.text.x =element_text(size=25),
        text = element_text(size=17))+scale_color_economist()


ggs_caterpillar(S.full)
ggs_caterpillar(S.full)+geom_vline(aes(yintercept=0),color="gray80",size=1,linetype="dashed")+
  theme_hc()+ 
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(size=25,vjust=0.75),
        axis.title.x=element_text(size=25,vjust=-0.25),axis.text.x =element_text(size=25),
        text = element_text(size=17))+aes(color=Parameter)+
  scale_color_manual(guide="none",values = blues)+ylab("Group ID")+
  xlab(expression(paste(zeta[j]," Highest Posterior Density"," ")))


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




