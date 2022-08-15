data{
  int<lower=0> N;
  vector[N] y;
}
parameters{
  real xi;
}
model{
    y ~ cauchy(xi, 0.5);
    xi ~ normal(0, 1);
}