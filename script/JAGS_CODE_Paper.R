# JAGS code paper

#  Required libraries
library(rjags);library(MASS);
library(parallel);library(rjags)
cl       <- makeCluster(3) # For parallel computing using 3 cores

############### Data
data <-read.csv("https://raw.githubusercontent.com/RafaelSdeSouza/Logit_Mar/master/data/Seyfert.csv",header=T)

X          <- model.matrix( ~  logM200 + r_r200, data = data) # Predictors
K          <- ncol(X)                   # Number of Predictors including the intercept 
y          <- data$bpt # Response variable (0/1)
n          <- length(y)                 # Sample size
gal        <- as.numeric(data$zoo) # Galaxy Type (0 = E, 1 = S)



############### JAGS Data
jags_data  <- list(Y= y,
                   N = n,
                   X=X,
                   gal=gal)

############### Model
jags_model<-"model{
# Shared hyperpriors for beta coefficients
tau ~ dgamma(0.001,0.001) # Precision
mu  ~ dnorm(0,1e-3) # mean

# Diffuse prior for beta coefficients
for(j in 1:2){
for(k in 1:3){
beta[k,j]~dnorm(mu,tau)                                    
}
}

# Likelihood
for(i in 1:N){
Y[i] ~ dbern(pi[i])
logit(pi[i]) <-  eta[i]
eta[i] <- beta[1,gal[i]]*X[i,1]+
beta[2,gal[i]]*X[i,2]+beta[3,gal[i]]*X[i,3]
}
}"
      
params <- c("beta") # Monitor these parameter  
inits0  <- function () {list(beta = matrix(rnorm(6,0, 0.01),ncol=2))} # A function to generat initial values for mcmc
inits1=inits0();inits2=inits0();inits3=inits0()         # Generate initial values for three chains

# Run mcmc
bin    = 10^4   # burn-in samples
ad     = 10^4   # Number of adaptive samples
s      = 3*10^4 # Number of samples for each chain
nc     = 3      # Number of mcmc
th     = 10     # Thinning value
jags_fit  <- run.jags(method="rjparallel",data = jags_data,inits = list(inits1,inits2,inits3),model=jags_model,
                        n.chains = nc,adapt=ad,monitor=c(params),burnin=bin,thin=th,sample=s,summarise=FALSE,plots=FALSE)

############### Output
print(jags_fit,intervals=c(0.025, 0.975), digits=3)





