# JAGS Code with Adaptive Shrinkage
#  Required libraries
library(rjags);library(ggmcmc);library(ggplot2);library(ggthemes);library(pander);library(Cairo);library(MASS);library(parallel) 
library(scales);library(plyr);require(gdata);require(runjags);require(gdata);require(caret);require(pROC);require(plyr);library(grid)
cl       <- makeCluster(3) # For parallel computing

source("facet_wrap_labeller.R")
# Read and format data
data     <- read.table("..//data/matched.txt",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","logM200_L","RprojLW_Rvir","zoo")]


#data_n2    <- subset(data_n, zoo=="E" | zoo == "S")

X          <- model.matrix( ~  logM200_L + RprojLW_Rvir, data = data_cut) # Predictors
K          <- ncol(X)                   # Number of Predictors including the intercept 
y          <- data_cut$bpt # Response variable (0/1)
n          <- length(y)                 # Sample size


# Grid of values for prediction 

Mx <- seq(1.75*min(X[,2]),1.75*max(X[,2]),length.out=250)
Rx <- seq(1.75*min(X[,3]),1.75*max(X[,3]),length.out=250)


jags.data  <- list(Y= y,N = n,X=X,b0 = rep(0,K),B0=diag(1e-4,K),Mx=Mx,
                   Rx=Rx)
# b0 and B0 are the mean vector and the inverse of the covariance matrix of prior distribution of the regression coefficients

### Prior Specification and Likelihood function
model<-"model{
#1. Priors

#a.Normal 
beta~dmnorm(b0,B0)                                    

#2. Likelihood
for(i in 1:N){
Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
eta[i] <- inprod(beta[], X[i,])
}

#3.Probability for each variable 
# 
for(l in 1:250){

logit(pMx[l])<-beta[1]+beta[2]*Mx[l]
logit(pRx[l])<-beta[1]+beta[3]*Rx[l]
logit(px[l])<-beta[1]+beta[2]*Mx[l]+beta[3]*Rx[l]
             }
             }"
params <- c("beta","pi","pMx","pRx","px")        # Monitor these parameters.
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
pMx    <- as.mcmc.list(jags.logit,vars="pMx")

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

gplot<-G1[,3:4]
gplot$Parameter<-as.factor(gplot$Parameter)
pL<-ggplot(data=gplot,aes(x=value,group=Parameter,fill=Parameter))+
  geom_density(colour="white",size=0.01,alpha=0.8)+facet_grid(Parameter~.,labeller = label_parsed)+
 theme_hc()+
theme(legend.position="none",panel.background = element_rect(fill = "white"),plot.background = element_rect(
    fill = "white"),plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=22),
        axis.text.y=element_text(size=22),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=22),axis.title.x=element_text(size=rel(1)),strip.background=element_rect(
          fill = "white"),strip.text=element_text(
            size = 25))+
#  geom_vline(xintercept=0,linetype="dashed",colour=c("grey80")) +
  scale_fill_manual(values=c("#E0FFFF","#00CED1","cyan4"))+
ylab("Density")+xlab("Parameter value")



G1<-ggs(jagssamples,family="beta")
plotbeta<-ggs_density(G1)+theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
   geom_vline(xintercept=0,linetype="dashed",colour=c("grey80")) + scale_fill_manual(values=c("#E0FFFF","#00CED1","cyan4"))
  aes(fill="#034e7b")+
  ylab("")+
  scale_y_discrete(breaks=c("beta[1]", "beta[2]", "beta[3]"),
                   labels=c(expression(beta[1]),expression(beta[2]),expression(beta[3])))

  
CairoPDF("betas.pdf",width = 5, height = 6)
plotbeta
dev.off()



# Plot all probabilities 
pi_AGN<-summary(as.mcmc.list(jags.logit, vars="pi"))
pi_AGN<-pi_AGN$quantiles

gdata<-data.frame(x=data_n2$sfr_tot_p50,mean=pi_AGN[,3],lwr1=pi_AGN[,2],lwr2=pi_AGN[,1],upr1=pi_AGN[,4],upr2=pi_AGN[,5])



pi_AGN<-summary(as.mcmc.list(jags.logit, vars="px"))
pi_AGN<-pi_AGN$quantiles
gdata1<-data.frame(x=data_n2$sfr_tot_p50,mean=pi_AGN[,3],lwr1=pi_AGN[,2],lwr2=pi_AGN[,1],upr1=pi_AGN[,4],upr2=pi_AGN[,5])



#-----------##-----------#-----------##-----------
pi_Mx<-summary(as.mcmc.list(jags.logit, vars="pMx"))
pi_Mx<-pi_Mx$quantiles
gMx<-data.frame(x=Mx,mean=pi_Mx[,3],lwr1=pi_Mx[,2],lwr2=pi_Mx[,1],upr1=pi_Mx[,4],upr2=pi_Mx[,5])
#-----------##-----------#-----------##-----------



#-----------##-----------#-----------##-----------
pi_Rx<-summary(as.mcmc.list(jags.logit, vars="pRx"))
pi_Rx<-pi_Rx$quantiles
gRx<-data.frame(x=Rx,mean=pi_Rx[,3],lwr1=pi_Rx[,2],lwr2=pi_Rx[,1],upr1=pi_Rx[,4],upr2=pi_Rx[,5])
#-----------##-----------#-----------##-----------






# Plots for each parameter


#b) M_star
PMx<-ggplot(aes(x=x,y=mean),data=gMx)+
  geom_ribbon(data=gMx,aes(x=x,y=mean,ymin=lwr2, ymax=upr2),alpha=0.95,  fill=c("#E0FFFF")) +
  geom_ribbon(data=gMx,aes(x=x,y=mean,ymin=lwr1, ymax=upr1),alpha=0.85,  fill=c("#00CED1")) +
  geom_line(size=1,linetype="dashed")+
#  geom_line(aes(x=x,y=mean),data=gMxS)+
#  geom_ribbon(data=gMxS,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.50, fill=c("#034e7b")) +
#  geom_ribbon(data=gMxS,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.40, fill=c("#a6bddb")) +
  theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(log~M[halo]))+ylab(expression(P[AGN]))
cairo_pdf("P_Mx.pdf",width = 8, height = 7)
PMx
dev.off()

#c) R_proj
PRx<-ggplot(aes(x=x,y=mean),data=gRx)+
  geom_ribbon(data=gRx,aes(x=x,y=mean,ymin=lwr2, ymax=upr2),alpha=0.95,  fill=c("#E0FFFF")) +
  geom_ribbon(data=gRx,aes(x=x,y=mean,ymin=lwr1, ymax=upr1),alpha=0.85,  fill=c("#00CED1")) +
  geom_line(size=1,linetype="dashed")+
  # geom_line(aes(x=x,y=mean),data=gRxS)+
#  geom_ribbon(data=gRxS,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.50, fill=c("#034e7b")) +
#  geom_ribbon(data=gRxS,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.40, fill=c("#a6bddb")) +
  theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(R/R[vir]))+ylab(expression(P[AGN]))
cairo_pdf("P_Rx.pdf",width = 8, height = 7)
PRx
dev.off()


# Plot 3D

z<-as.numeric(pi_AGN[,3])
x<-Mx
y<-Rx

library(rsm)
library(lattice)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
#p<-persp(x,y,z, theta=150, phi=20, 
#         expand = 0.5,shade = 0.1,
#         xlab="Z", ylab=expression(NII.Ha), zlab=expression(log10.EW.Ha),ticktype='detailed',
#         col = YlOrBr,border=NA,xlog=T,ylog=T)
cor = topo.colors(200)

cairo_pdf("logit3D.pdf")
trellis.par.set("axis.line",list(axis.text=list(cex=20),col=NA,lty=1,lwd=2))
par(mar=c(1,1,1,1))
wireframe(z~x+y,data=data.frame(x=x, y=rep(y, each=length(x)), z=z),
          par.settings = list(regions=list(alpha=0.4)),
          col.regions =cor,drape=T,light.source = c(5,5,5),colorkey = FALSE,
          xlab=list(label=expression(log10.NII.Ha.),cex=1.25), ylab=list(label=expression(log10.EW.Ha..),cex=1.25), 
          zlab=list(rot=90,label=expression(pi),cex=1.25,dist=-1,rot=0),
          scale=list(tck=0.75,arrows=FALSE,distance =c(0.75, 0.75, 0.75)))


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
