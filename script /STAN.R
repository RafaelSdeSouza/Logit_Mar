# Stan version speed-up
require(rstan)



# Set Seed
set.seed(1234)

# Number of Households
H=100

# True parameters
true.b.mu = c(10,-5)
true.b.sig = c(3,2)
true.y.sig = 3

# Storage 
b =  matrix(0,100,2)
ctr = 0
# Each "household" can have between 5 and 10 reps
Ti = sample(5:10,H,replace=T)
id = x = y = matrix(0,sum(Ti),1)

# Simulate Data
for(i in 1:100) {
  b[i,1] = rnorm(1,true.b.mu[1],true.b.sig[1])
  b[i,2] = rnorm(1,true.b.mu[2],true.b.sig[2])	
  for(j in 1:Ti[i]) {
    ctr = ctr + 1
    x[ctr] = runif(1)*3 + 1
    y[ctr] = b[i,1]+b[i,2]*x[ctr]+rnorm(1,0,true.y.sig)
    id[ctr] = i 
  }
}



hcode = "
data {
int<lower=0> N; // Observations
int<lower=0> H; // persons
int<lower=0> id[N]; // ID variable
real y[N]; // Y Vector
real x[N]; // X Vector
}

parameters {
vector[2] beta[H];	// [H,2] dim matrix for b
real beta_mu[2]; // means for b
real<lower=0> beta_sig[2];   // var parameters for b
real<lower=0> y_sig;   // regression variance
}

model {
//Priors
y_sig~uniform(0,1000); 

for(j in 1:2) {
beta_mu[j] ~ normal(0,100);
beta_sig[j] ~ uniform(0,1000);
for(h in 1:H) 
beta[h,j] ~ normal(beta_mu[j], beta_sig[j]); 
}

// Likelihood
for(n in 1:N) 
y[n] ~ normal(beta[id[n],1]+beta[id[n],2] .* x[n], y_sig);

} 
"


data ="{
  int<lower=0> N; // Observations
  int<lower=0> H; // persons
  int<lower=0> id[N]; // ID variable
  real y[N]; // Y Vector
  real x[N]; // X Vector
}"



parameters=" {
  vector[2] beta[H];	// [H,2] dim matrix for b
  real beta_mu[2]; // means for b
  real<lower=0> beta_sig[2];   // var parameters for b
  real<lower=0> y_sig;   // regression variance
}"


model = "{
#Priors
  y_sig~uniform(0,1000); 
  
  for(j in 1:2) {
    beta_mu[j] ~ normal(0,100);
    beta_sig[j] ~ uniform(0,1000);
    for(h in 1:H) 
      beta[h,j] ~ normal(beta_mu[j], beta_sig[j]); 
  }
# Likelihood
  for(n in 1:N) 
    y[n] ~ normal(beta[id[n],1]+beta[id[n],2] .* x[n], y_sig);
  
}"


dat = list(N=length(y),H=dim(b)[1],y=y,x=x,id=id)
hlm = stan(model_name="Hierarchical Linear Model", model_code = hcode, data=dat , iter = 2000, n_chains = 2, verbose = FALSE)

