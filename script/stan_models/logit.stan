data {
    int<lower=0> n; // number of data points
    int<lower=1> k; // number of variables
    matrix[n, k] x; // predictor variables
    int<lower=0,upper=1> y[N]; // dependent variable
}
parameters {
    real alpha;     // intercept
    vector[k] beta; // coefficients for predictors    
}
model {
    y ~ bernoulli_logit(alpha + beta * x);
}
