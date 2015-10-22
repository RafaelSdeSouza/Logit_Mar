# JAGS Code with Adaptive Shrinkage
#  Required libraries
library(rjags);library(ggmcmc);library(ggplot2);library(ggthemes);library(pander);library(Cairo);library(MASS);library(parallel) 
library(scales);library(plyr);require(gdata);require(runjags);require(gdata);require(caret);require(pROC);require(plyr)
cl       <- makeCluster(3) # For parallel computing


# Read and format data
data     <- read.csv("..//data/sample_CRP02_sub.csv",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","lgm_tot_p50","logM200_L","RprojLW_Rvir","sfr_tot_p50","zoo")]
data2    <- na.omit(data_cut)
data2    <- data2[data2$lgm_tot_p50>0,]
data2    <- data2[which(data2$logM200_L>0),]
data2    <- data2[which(data2$RprojLW_Rvir>=0),]
data2    <- data2[which(data2$sfr_tot_p50>=-100),]

trainIndex <- sample(1:nrow(data2),1000)
data3      <- data2[trainIndex,]
data3    <- subset(data3, bpt!="LINER")    # remove Liners
data3$bpt  <- revalue(data3$bpt,c("Star Forming"="0","Composite"="0",
                                 "LINER"="1","Seyfert/LINER"="1","Star Fo"="0",
                                 "Seyfert"="1","BLANK"="0"))



data_n     <- data3
data_n2    <- subset(data_n, zoo=="E" | zoo == "S")
galtype    <- match(data_n2$zoo,c("E", "S")) 
galtype    <- galtype-1                 # 0 = E and 1 = S
X          <- model.matrix(~lgm_tot_p50+logM200_L+RprojLW_Rvir+sfr_tot_p50+galtype,data=data_n2) # Predictors
K          <- ncol(X)                   # Number of Predictors including the intercept 
y          <- as.numeric(data_n2$bpt)-1 # Response variable (0/1)
n          <- length(y)                 # Sample size

jags.data  <- list(Y= y,N = n,X=X,b0 = rep(0,K),B0=diag(1e-4,K))
# b0 and B0 are the mean vector and the inverse of the covariance matrix of prior distribution of the regression coefficients

### Prior Specification and Likelihood function
model<-"model{
#1. Priors
beta~dmnorm(b0,B0) 
#2. Likelihood
for(i in 1:N){
Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
eta[i] <- inprod(beta, X[i,])
             }
             }"
params <- c("beta","pi")                                # Monitor these parameters.
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

G1<-ggs(jagssamples,family="beta")
ggs_caterpillar(G1)+theme_few()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
  scale_color_fivethirtyeight()+
  scale_fill_fivethirtyeight()


