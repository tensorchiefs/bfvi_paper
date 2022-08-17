 data {
    int N;   //Number of Datapoints
    int M;   //Number of features (here 1)
    real X[N, M];
    int<lower=0, upper=1> y[N];
}

parameters {
    real beta[M];
    real beta0;
}

model {
    for (i in 1:N)
        y[i] ~ bernoulli(inv_logit (beta0 + dot_product(X[i] , beta)));
    beta[M] ~ normal(0, 1);
    beta0 ~ normal(0, 1);
}

