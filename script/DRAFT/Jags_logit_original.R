# JAGS Code with Adaptive Shrinkage
#  Required libraries
library(rjags);library(ggmcmc);library(ggplot2);library(ggthemes);library(pander);library(Cairo);library(MASS);library(parallel) 
library(scales);library(plyr);require(gdata);require(runjags);require(gdata);require(caret);require(pROC);require(plyr)
cl       <- makeCluster(3) # For parallel computing


# Read and format data
data     <- read.csv("..//data/sample_CRP02_sub.csv",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","lgm_tot_p50","logM200_L","RprojLW_Rvir","sfr_tot_p50","color_gr","zoo")]
data2    <- na.omit(data_cut)
data2    <- data2[data2$lgm_tot_p50>0,]
data2    <- data2[which(data2$logM200_L>0),]
data2    <- data2[which(data2$RprojLW_Rvir>=0),]
data2    <- data2[which(data2$sfr_tot_p50>=-100),]

# Standardized variables

data2<-data.frame(bpt=data2$bpt,as.data.frame(scale(data2[,2:6])),zoo=data2$zoo)

trainIndex <- sample(1:nrow(data2),2000)
data3      <- data2[trainIndex,]
#data3    <- subset(data3, bpt!="LINER")    # remove Liners
data3$bpt  <- revalue(data3$bpt,c("Star Forming"="0","Composite"="0",
                                 "LINER"="1","Seyfert/LINER"="1","Star Fo"="0",
                                 "Seyfert"="1","BLANK"="0"))



data_n     <- data3
data_n2    <- subset(data_n, zoo=="E" | zoo == "S")
galtype    <- match(data_n2$zoo,c("E", "S")) 
galtype    <- galtype-1                 # 0 = E and 1 = S
X          <- model.matrix( ~ lgm_tot_p50 + logM200_L + RprojLW_Rvir + 
                              sfr_tot_p50 + color_gr + galtype, data = data_n2) # Predictors
K          <- ncol(X)                   # Number of Predictors including the intercept 
y          <- as.numeric(data_n2$bpt)-1 # Response variable (0/1)
n          <- length(y)                 # Sample size


# Grid of values for prediction 
gx <- seq(1.75*min(X[,2]),1.75*max(X[,2]),length.out=250)
Mx <- seq(1.75*min(X[,2]),1.75*max(X[,3]),length.out=250)
Rx <- seq(1.75*min(X[,2]),1.75*max(X[,4]),length.out=250)
sfrx <- seq(1.75*min(X[,2]),1.75*max(X[,5]),length.out=250)
grx <-  seq(1.75*min(X[,2]),1.75*max(X[,6]),length.out=250) 
 

jags.data  <- list(Y= y,N = n,X=X,b0 = rep(0,K),B0=diag(1e-4,K),gx=gx,Mx=Mx,
                   Rx=Rx,sfrx=sfrx, grx = grx )
# b0 and B0 are the mean vector and the inverse of the covariance matrix of prior distribution of the regression coefficients

### Prior Specification and Likelihood function
model<-"model{
#1. Priors

#a.Normal 
beta~dmnorm(b0,B0)                                    

#b.Jefreys priors for sparseness
#for(j in 1:K)   {
#lnTau[j] ~ dunif(-50, 50)
#TauM[j] <- exp(lnTau[j])
#beta[j] ~ dnorm(0, TauM[j])
#Ind[j] <- step(abs(beta[j]) - 0.05)
#}

#c.LASSO

#for(j in 1:K){
#beta[j]~ddexp(0,tau)
#}
#c1.prior for tau
#tau <-pow(sdBeta,-1)
#sdBeta ~ dgamma(0.01,0.01)

#2. Likelihood
for(i in 1:N){
Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
eta[i] <- inprod(beta[], X[i,])
}

#3.Probability for each variable 
# E galaxies
for(l in 1:250){
logit(pgx[l])<-beta[1]+beta[2]*gx[l]
logit(pMx[l])<-beta[1]+beta[3]*Mx[l]
logit(pRx[l])<-beta[1]+beta[4]*Rx[l]
logit(psfrx[l])<-beta[1]+beta[5]*sfrx[l]
logit(pgrx[l])<-beta[1]+beta[6]*grx[l]

# S galaxies
logit(pgxS[l])<-beta[1]+beta[2]*gx[l]+beta[7]
logit(pMxS[l])<-beta[1]+beta[3]*Mx[l]+beta[7]
logit(pRxS[l])<-beta[1]+beta[4]*Rx[l]+beta[7]
logit(psfrS[l])<-beta[1]+beta[5]*sfrx[l]+beta[7]
logit(pgrS[l])<-beta[1]+beta[6]*grx[l]+beta[7]

             }
             }"
params <- c("beta","pi","pgx","pMx","pRx","psfrx","pgrx",
            "pgxS","pMxS","pRxS","psfrS","pgrS")        # Monitor these parameters.
inits0  <- function () {list(beta = rnorm(K, 0, 0.01))} # A function to generat initial values for mcmc
inits1=inits0();inits2=inits0();inits3=inits0()         # Generate initial values for three chains

# Run mcmc
bin    = 10^4   # burn-in samples
ad     = 10^4   # Number of adaptive samples
s      = 3*10^4 # Number of samples for each chain
nc     = 3      # Number of mcmc
th     = 10     # Thinning value
jags.logit  <- run.jags(method="rjparallel",data = jags.data,inits = list(inits1,inits2,inits3),model=model,
                       n.chains = nc,adapt=ad,monitor=c(params),burnin=bin,thin=th,sample=s,summarise=FALSE,plots=FALSE)

beta     <- as.mcmc.list(jags.logit,vars="beta")
pi       <- as.mcmc.list(jags.logit,vars="pi")
pgx     <- as.mcmc.list(jags.logit,vars="pgx")

# Trace plots and diagnostic analysis to investigate convregence
plot(beta) 
gelman.diag(beta)
gelman.plot(beta)
geweke.diag(beta)
geweke.plot(beta)
autocorr.plot(beta)

# Summary Statistics for parameters of interest
print(summary(window(beta)))
print(summary(window(pi)))


# Plots with ggmcmc

jagssamples <- as.mcmc.list(jags.logit)

labels<-c(expression(beta[0]),expression(beta[1]),expression(beta[2]),
          expression(beta[3]),expression(beta[4]),expression(beta[5]),expression(beta[6]))
L.radon.intercepts <- data.frame(
  Parameter=paste("beta[", seq(1:7), "]", sep=""),
  Label=labels)
head(L.radon.intercepts)

#"lgm_tot_p50","logM200_L","RprojLW_Rvir","sfr_tot_p50","color_gr","zoo"

G1<-ggs(jagssamples,family="beta")
plotbeta<-ggs_caterpillar(G1)+theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
   geom_vline(xintercept=0,linetype="dashed",colour=c("#034e7b")) +
  aes(color="#034e7b")+ylab("")+
  scale_y_discrete(breaks=c("beta[1]", "beta[2]", "beta[3]","beta[4]","beta[5]","beta[6]","beta[7]"),
                   labels=c(expression(beta[1]),expression(beta[2]),expression(beta[3]),
                   expression(beta[4]),expression(beta[5]),expression(beta[6]),expression(beta[7])))

CairoPDF("betas.pdf",width = 5, height = 6)
plotbeta
dev.off()



# Plot all probabilities 
pi_AGN<-summary(as.mcmc.list(jags.logit, vars="pi"))
pi_AGN<-pi_AGN$quantiles

gdata<-data.frame(x=data_n2$sfr_tot_p50,mean=pi_AGN[,3],lwr1=pi_AGN[,2],lwr2=pi_AGN[,1],upr1=pi_AGN[,4],upr2=pi_AGN[,5])



pi_AGN<-summary(as.mcmc.list(jags.logit, vars="pi"))
pi_AGN<-pi_AGN$quantiles
gdata1<-data.frame(x=data_n2$sfr_tot_p50,mean=pi_AGN[,3],lwr1=pi_AGN[,2],lwr2=pi_AGN[,1],upr1=pi_AGN[,4],upr2=pi_AGN[,5])

pi_gx<-summary(as.mcmc.list(jags.logit, vars="pgx"))
pi_gx<-pi_gx$quantiles
ggx<-data.frame(x=sfrx,mean=pi_sfrS[,3],lwr1=pi_sfrS[,2],lwr2=pi_sfrS[,1],upr1=pi_sfrS[,4],upr2=pi_sfrS[,5])

pi_gx<-summary(as.mcmc.list(jags.logit, vars="pgx"))
pi_gx<-pi_gx$quantiles
ggx<-data.frame(x=sfrx,mean=pi_sfrS[,3],lwr1=pi_sfrS[,2],lwr2=pi_sfrS[,1],upr1=pi_sfrS[,4],upr2=pi_sfrS[,5])


#-----------##-----------#-----------##-----------
pi_Mx<-summary(as.mcmc.list(jags.logit, vars="pMx"))
pi_Mx<-pi_Mx$quantiles
gMx<-data.frame(x=Mx,mean=pi_Mx[,3],lwr1=pi_Mx[,2],lwr2=pi_Mx[,1],upr1=pi_Mx[,4],upr2=pi_Mx[,5])

pi_MxS<-summary(as.mcmc.list(jags.logit, vars="pMxS"))
pi_MxS<-pi_MxS$quantiles
gMxS<-data.frame(x=Mx,mean=pi_MxS[,3],lwr1=pi_MxS[,2],lwr2=pi_MxS[,1],upr1=pi_MxS[,4],upr2=pi_MxS[,5])
#-----------##-----------#-----------##-----------


pi_Rx<-summary(as.mcmc.list(jags.logit, vars="pRx"))
pi_Rx<-pi_Rx$quantiles
gRx<-data.frame(x=Rx,mean=pi_Rx[,3],lwr1=pi_Rx[,2],lwr2=pi_Rx[,1],upr1=pi_Rx[,4],upr2=pi_Rx[,5])


pi_RxS<-summary(as.mcmc.list(jags.logit, vars="pRxS"))
pi_RxS<-pi_RxS$quantiles
gRxS<-data.frame(x=Rx,mean=pi_RxS[,3],lwr1=pi_RxS[,2],lwr2=pi_RxS[,1],upr1=pi_RxS[,4],upr2=pi_RxS[,5])

#-----------##-----------#-----------##-----------
pi_sfrx<-summary(as.mcmc.list(jags.logit, vars="psfrx"))
pi_sfrx<-pi_sfrx$quantiles
gsfr<-data.frame(x=sfrx,mean=pi_sfrx[,3],lwr1=pi_sfrx[,2],lwr2=pi_sfrx[,1],upr1=pi_sfrx[,4],upr2=pi_sfrx[,5])

pi_sfrS<-summary(as.mcmc.list(jags.logit, vars="psfrS"))
pi_sfrS<-pi_sfrS$quantiles
gsfrS<-data.frame(x=sfrx,mean=pi_sfrS[,3],lwr1=pi_sfrS[,2],lwr2=pi_sfrS[,1],upr1=pi_sfrS[,4],upr2=pi_sfrS[,5])
#-----------##-----------#-----------##-----------


pi_grx<-summary(as.mcmc.list(jags.logit, vars="pgrx"))
pi_grx<-pi_grx$quantiles

#p_all<-rbind(pi_gx,pi_Mx,pi_Rx)
#p_all2<-cbind(p_all,data.frame(rep(c("gal","Mhalo","Rproj"),each=100)))
#colnames(p_all2)<-c("2.5","25","50","75","97.5","probs")

# Plots for each parameter

#a) SFR
Psfr<-ggplot(aes(x=x,y=mean),data=gsfr)+geom_line()+
  geom_ribbon(data=gsfr,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.50, fill=c("#d7301f")) +
  geom_ribbon(data=gsfr,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.40, fill=c("#feb24c")) +
  geom_line(aes(x=x,y=mean),data=gsfrS)+
  geom_ribbon(data=gsfrS,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.50, fill=c("#034e7b")) +
  geom_ribbon(data=gsfrS,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.40, fill=c("#a6bddb")) +
  theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(log~SFR))+ylab(expression(P[AGN]))
cairo_pdf("P_sfr.pdf",width = 8, height = 7)
Psfr
dev.off()

#b) M_gal
PMx<-ggplot(aes(x=x,y=mean),data=gMx)+geom_line()+
  geom_ribbon(data=gMx,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.50, fill=c("#d7301f")) +
  geom_ribbon(data=gMx,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.40, fill=c("#feb24c")) +
  geom_line(aes(x=x,y=mean),data=gMxS)+
  geom_ribbon(data=gMxS,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.50, fill=c("#034e7b")) +
  geom_ribbon(data=gMxS,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.40, fill=c("#a6bddb")) +
  theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(M[gal]))+ylab(expression(P[AGN]))
cairo_pdf("P_Mx.pdf",width = 8, height = 7)
PMx
dev.off()

#c) R_proj
PRx<-ggplot(aes(x=x,y=mean),data=gRx)+geom_line()+
  geom_ribbon(data=gRx,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.50, fill=c("#d7301f")) +
  geom_ribbon(data=gRx,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.40, fill=c("#feb24c")) +
#  geom_line(aes(x=x,y=mean),data=gRxS)+
#  geom_ribbon(data=gRxS,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.50, fill=c("#034e7b")) +
#  geom_ribbon(data=gRxS,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.40, fill=c("#a6bddb")) +
  theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(R/R[vir]))+ylab(expression(P[AGN]))
cairo_pdf("P_Rx.pdf",width = 8, height = 7)
PRx
dev.off()




#gdata_pall<-data.frame(x=data_n2$sfr_tot_p50,mean=pi_AGN[,3],lwr1=pi_AGN[,2],lwr2=pi_AGN[,1],upr1=pi_AGN[,4],upr2=pi_AGN[,5])
# Plot fit 

P1<-ggplot(aes(x=x,y=mean),data=gdata)+geom_line()+
  geom_ribbon(data=gdata,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.45, fill=c("#005502")) +
  geom_ribbon(data=gdata,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.35, fill=c("#3A5F0B")) +
  theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(SFR))+ylab(expression(P[AGN]))

cairo_pdf("logit_AGN.pdf",width = 8, height = 7)
P1
dev.off()
