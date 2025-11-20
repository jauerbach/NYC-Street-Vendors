functions {
  real nb_lpmf(int n, real mu, real rho) {
    real r = mu / rho;
    real q = 1.0 - 1.0 / rho;
    if (r <= 0 || q <= 0 || q >= 1) return negative_infinity();
    if (n - r + 1 <= 0)             return negative_infinity();
    return lgamma(n) - lgamma(r) - lgamma(n - r + 1)
         + (n - r) * log(q) + r * log1m(q);
  }
  real mnh_lpmf(array[] int n_aug, vector r) {
    int J = size(n_aug);
    int N = 0; for (j in 1:J) N += n_aug[j];
    real R = sum(r);
    if (R <= 0 || N - R + 1 <= 0) return negative_infinity();
    real lp = 0;
    for (j in 1:J) {
      if (n_aug[j] <= 0 || r[j] <= 0) return negative_infinity();
      if (n_aug[j] - r[j] + 1 <= 0)   return negative_infinity();
      lp += lgamma(n_aug[j]) - lgamma(r[j]) - lgamma(n_aug[j] - r[j] + 1);
    }
    lp -= (lgamma(N) - lgamma(R) - lgamma(N - R + 1));
    return lp;
  }
}

data {
  int<lower=1> K;            // number of subregions
  array[K] int<lower=0> n0;  // number of vendors without permits in each subregion
  array[K] int<lower=0> n1;  // number of vendors with permits in each subregion
  int<lower=0> N1;           // known number of permits/licenses citywide
  vector<lower=1>[K] rho;    // expected proportion of respondents per market
}

parameters {
  real mu_p;
  real<lower=0> sigma_p;
  vector[K] alpha_raw;
  real mu_0;
  real<lower=0> sigma_0;
  vector[K] eta0_raw;
  real<lower=0> sigma_1;
  vector[K] eta1_raw;
}

transformed parameters {
  vector[K] alpha = mu_p + sigma_p * alpha_raw;
  vector<lower=0, upper=1>[K] p = inv_logit(alpha);
  vector[K] eta0 = mu_0 + sigma_0 * eta0_raw;
  vector<lower=0>[K] lambda0 = exp(eta0);
  vector<lower=0>[K] mu0 = p .* lambda0;
  vector<lower=0>[K] M = mu0 ./ rho;
  vector[K] eta1 = sigma_1 * eta1_raw;
  vector<lower=0>[K] lambda1 = exp(eta1);
  vector<lower=0>[K+1] s;
  s[1:K] = sum(lambda0 ./ rho) * p .* (lambda1 / sum(lambda1));
  s[K+1] = sum(lambda0 ./ rho) * (1 - sum(p .* (lambda1 / sum(lambda1))));
}

model {
  alpha_raw ~ normal(0, 1);
  eta0_raw  ~ normal(0, 1);
  eta1_raw  ~ normal(0, 1);

  for (i in 1:K)
    target += nb_lpmf(n0[i] | mu0[i], rho[i]);

  array[K+1] int n1_aug;
  n1_aug[1:K] = n1;
  n1_aug[K+1] = N1 - sum(n1);
  target += mnh_lpmf(n1_aug | s);
}
