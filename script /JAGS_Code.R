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
data<-read.csv("..//data/Logit_sample-2.csv",header=TRUE,na.strings="")
data2<-na.omit(data)
data2<-data2[data2$logMstar_p50>0,]
data2<-data2[which(data2$sigma>0),]
data2<-data2[which(data2$rvir>0),]

trainIndex <- sample(1:nrow(data2),5000)
data3<-data2[trainIndex,]

#data2$bpt<-as.factor(data2$bpt)
data3$bpt <- revalue(data3$bpt,c("Star Forming"="0","Composite"="0","Composi"="0",
                                 "LINER"="1","Seyfert/LINER"="1","Star Fo"="0",
                                 "Seyfert"="1","BLANK"="0"))

data3$R_rvir = data3$Rproj_L/data3$rvir
# Prepare data for JAGS
data3$logMstar_p50<-(data3$logMstar_p50-mean(data3$logMstar_p50))/sd(data3$logMstar_p50)
data3$logMhalo<-(data3$logMhalo-mean(data3$logMhalo))/sd(data3$logMhalo)
#data2$vlos_sigma<-(data2$vlos_sigma-mean(data2$vlos_sigma))/sd(data2$vlos_sigma)
data3$R_rvir<-(data3$R_rvir-mean(data3$R_rvir))/sd(data3$R_rvir)


X<-model.matrix(~logMstar_p50+logMhalo+R_rvir,data=data3)
K<-ncol(X)

jags.data <- list(Y= as.numeric(data3$bpt)-1,
                  N = nrow(data3),
                  X=X,
#                  b0 = rep(0,K),
#                  B0=diag(1e-4,K),
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
                       adapt=1000,
                       monitor=c(params),
                       burnin=5000,
                       sample=15000,
                       summarise=FALSE,
                       plots=FALSE
)

jagssamples <- as.mcmc.list(jags.logit)

G1<-ggs(jagssamples)
ggs_density(G1)+theme_few()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
  scale_colour_manual(values=c("#00526D", "#00A3DB", "#7A2713", "#939598", "#6CCFF6"))+
  scale_fill_manual(values=c("#00526D", "#00A3DB", "#7A2713", "#939598", "#6CCFF6"))

# Plot fit 

P1<-ggplot(aes(x=x,y=y),data=Msigma2)+geom_point()+
  theme_few()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
  scale_colour_manual(values=c("#00526D", "#00A3DB", "#7A2713", "#939598", "#6CCFF6"))+
  scale_fill_manual(values=c("#00526D", "#00A3DB", "#7A2713", "#939598", "#6CCFF6"))


s<-summary(posterior.normal)
capture.output(s, file = "myfile.txt")


