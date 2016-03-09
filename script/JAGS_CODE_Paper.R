# JAGS Script

#  Required libraries
library(R2jags)

############### Data
data <-read.csv("https://goo.gl/gvFSH8",header=T)

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
tau ~ dgamma(1e-3,1e-3) # Precision
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
      
params <- c("beta") # Monitor this parameter  
inits  <- function () {list(beta = matrix(rnorm(6,0, 0.01),ncol=2))} # A function to generate initial values for mcmc


# Run mcmc
burn    = 2*10^4   # burn-in samples
s       = 5*10^4   # Number of samples for each chain
nc      = 3        # Number of mcmc
th      = 10        # Thinning value

# JAGS MCMC
jags_fit <- jags(data = jags_data,
                 inits      = inits,
                 parameters = params,
                 model.file = textConnection(jags_model),
                 n.chains   = nc,
                 n.thin     = th,
                 n.iter     = s,
                 n.burnin   = burn)

############### Output
print(jags_fit,intervals=c(0.025, 0.975), digits=3)





