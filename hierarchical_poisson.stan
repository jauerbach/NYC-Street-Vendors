data {
  int<lower=1> K;            // number of subregions
  array[K] int<lower=0> n0;  // number of vendors without permits in each subregion
  array[K] int<lower=0> n1;  // number of vendors with permits in each subregion
  int<lower=0> N1;           // known number of permits/licenses citywide
}

parameters {
  real mu_p;
  real<lower=0> sigma_p;
  vector[K] alpha_raw;
  real mu_0;
  real<lower=0> sigma_0;
  vector[K] eta0_raw;
  real<lower=0> sigma_1;
  vector[K-1] eta1_raw;
}

transformed parameters {
  vector[K] alpha = mu_p + sigma_p * alpha_raw;
  vector<lower=0, upper=1>[K] p = inv_logit(alpha);
  vector[K] eta0 = mu_0 + sigma_0 * eta0_raw;
  vector<lower=0>[K] lambda0 = exp(eta0);
  vector[K] eta1;
  eta1[1:(K-1)] = sigma_1 * eta1_raw;
  eta1[K]       = -sum(eta1[1:(K-1)]);
  simplex[K] r = softmax(eta1);                      
}

model {
  alpha_raw ~ normal(0,1);
  eta0_raw ~ normal(0,1);
  eta1_raw ~ normal(0, 1);
  
  for (i in 1:K)
    n0[i] ~ poisson(p[i] * lambda0[i]);
  
  vector[K+1] theta;
  theta[1:K] = p .* r;
  theta[K+1] = 1 - sum(p .* r);

  array[K+1] int n1_aug;
  n1_aug[1:K] = n1;
  n1_aug[K+1] = N1 - sum(n1);

  target += multinomial_lpmf(n1_aug | theta);
}
