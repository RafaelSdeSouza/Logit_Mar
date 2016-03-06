# One dimensional 

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
data<-read.csv("..//data/sample_CRP02_sub.csv",header=TRUE,na.strings="")
#data$tempbar<-data$t03_bar_a06_bar_flag*2+data$t03_bar_a07_no_bar_flag
#data_cut<-data[,c("bpt","lgm_tot_p50","logM200_L","RprojLW_Rvir","sfr_tot_p50","zoo","tempbar")]
data_cut<-data[,c("bpt","lgm_tot_p50","logM200_L","RprojLW_Rvir","sfr_tot_p50","zoo")]
#data_cut<-data_cut[which(data$tempbar==1 | data$tempbar==2),]



data2<-na.omit(data_cut)
data2<-data2[data2$lgm_tot_p50>0,]
data2<-data2[which(data2$logM200_L>0),]
data2<-data2[which(data2$RprojLW_Rvir>=0),]
data2<-data2[which(data2$sfr_tot_p50>=-100),]

# bar 


#


trainIndex <- sample(1:nrow(data2),10000)
data3<-data2[trainIndex]

#data2$bpt<-as.factor(data2$bpt)
data3$bpt <- revalue(data3$bpt,c("Star Forming"="0","Composite"="0",
                                 "LINER"="1","Seyfert/LINER"="1","Star Fo"="0",
                                 "Seyfert"="1","BLANK"="0"))




# Standardized variables

data_n<-data3
#data_n<-data.frame(data3$bpt,as.data.frame(scale(data3[,2:5])),zoo=data3$zoo)


#data3$R_rvir = data3$Rproj_L/data3$rvir
# Prepare data for JAGS
#data3$lgm_tot_p50<-(data3$lgm_tot_p50-mean(data3$lgm_tot_p50))/sd(data3$lgm_tot_p50)
#data3$logMhalo<-(data3$logMhalo-mean(data3$logMhalo))/sd(data3$logMhalo)
#data2$vlos_sigma<-(data2$vlos_sigma-mean(data2$vlos_sigma))/sd(data2$vlos_sigma)
#data3$R_rvir<-(data3$R_rvir-mean(data3$R_rvir))/sd(data3$R_rvir)

data_n2<-subset(data_n, zoo=="E" | zoo == "S")
galtype<-match(data_n2$zoo,c("E", "S"))
#galtype<-as.numeric(data_n2$zoo)-1
Ntype<-length(unique(data_n2$zoo))


X<-model.matrix(~sfr_tot_p50,data=data_n2)
K<-ncol(X)

data_n2$zoo<-droplevels(data_n2$zoo)

jags.data <- list(Y= as.numeric(data_n2$bpt)-1,
                  N = nrow(data_n2),
                  X=X,
                                    b0 = rep(0,K),
                                   B0=diag(1e-4,K),
                  galtype = galtype,
                  #                 bar=bar,
                  Ntype=Ntype,
                  Npred = K
)
model<-"model{
#1. Priors
beta~dmnorm(b0[],B0[,]) # Normal Priors
# Jefreys priors for sparseness
#for(j in 1:Npred)   {
#lnTau[j] ~ dunif(-50, 50)
#TauM[j] <- exp(lnTau[j])
#beta[j] ~ dnorm(0, TauM[j])
#Ind[j] <- step(abs(beta[j]) - 0.05)
#}

#LASSO

#for(j in 1:Npred){
#beta[j]~ddexp(0,tau)
#}
#prior for tau
#tau <-pow(sdBeta,-1)
#sdBeta ~ dgamma(0.01,0.01)


# Random intercept 
tau2~dgamma(1e-3,1e-3)
ranef[1]~dbern(0)
ranef[2]~dnorm(0,1/tau2)

#for (j in 2:Ntype){
#ranef[j]~dnorm(0,1/tau2)
#}


#2. Likelihood

for(i in 1:N){
Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
eta[i] <- inprod(beta[], X[i,])+ranef[galtype[i]]


#3. Prediction
NewPred[i]~dbern(pi[i])
}

}"

#params <- c("beta","ranef")
params <- c("beta","pi","Ind","NewPred")

inits0  <- function () {
  list(beta = rnorm(K, 0, 0.01))}



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
                       adapt=5000,
                       monitor=c(params),
                       burnin=1000,
                       sample=30000,
                       summarise=FALSE,
                       plots=FALSE
)

jagssamples <- as.mcmc.list(jags.logit)

G1<-ggs(jagssamples,family="beta")
ggs_density(G1)+theme_few()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
  scale_color_stata()+
  scale_fill_stata()



pi_AGN<-summary(as.mcmc.list(jags.logit, vars="pi"))
pi_AGN<-pi_AGN$quantiles

gdata<-data.frame(x=data_n2$sfr_tot_p50,mean=pi_AGN[,3],lwr1=pi_AGN[,2],lwr2=pi_AGN[,1],upr1=pi_AGN[,4],upr2=pi_AGN[,5])



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