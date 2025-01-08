data {
  int<lower=1> N;
  vector[N] effect;
  vector[N] SE;
  vector[N] study_id;
}
parameters {
  real zeta;
  vector[N] zeta_n;
  real<lower = 0> tau;
}
model {
  target += normal_lpdf(zeta_n | zeta, tau);
  target += normal_lpdf(effect| zeta + zeta_n, sqrt(square(SE) + square(tau)));
  target += normal_lpdf(zeta | 0, 100);
  target += normal_lpdf(tau | 0, 100)
    - normal_lccdf(0 | 0, 100);
}
