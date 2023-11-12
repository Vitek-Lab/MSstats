data {
  int<lower=0> N_obs;
  int<lower=0> N_missing;
  int<lower=0> R;
  int<lower=0> Feat;
  int<lower=0> P;

  vector[N_obs] obs;
  int run_id[N_obs];
  int feature_id[N_obs];
  int protein_id[N_obs]; 
  
  int run_id_missing[N_missing];
  int feature_id_missing[N_missing];
  int protein_id_missing[N_missing];

  int zeros[N_obs] ;
  int ones[N_missing];
  real<lower=0, upper=40> run_priors[R];
  vector<lower=0>[P] sigma;
  vector<lower=0>[R] sigma_run;
  vector<lower=0>[Feat] sigma_feat;
  
  real beta0;
  real beta1;
}
parameters {
  // real beta0;
  // real beta1;
  // real mar;
  vector<lower=0, upper=50>[R] run_mu;
  vector<lower=-10, upper=10>[Feat] feature_mu;
  // vector<lower=0, upper=5>[P] sigma;
  vector<lower=0, upper=50>[N_missing] obs_mis;
  // vector<lower=0, upper=5>[N_missing] missing_sigma;
}
model {
  // beta0 ~ normal(8., 2.);
  // beta1 ~ normal(.8, .2);
  // mar ~ beta(1, 20);
  
  sigma_run ~ lognormal(1,1);
  run_priors ~ normal(run_priors, 2.);
  run_mu ~ normal(run_priors, sigma_run);
  feature_mu ~ normal(0, sigma_feat);
  // sigma ~ exponential(sigma_list);
  
  obs ~ normal(run_mu[run_id] + feature_mu[feature_id], sigma[protein_id]);
  obs_mis ~ normal(run_mu[run_id_missing] + feature_mu[feature_id_missing] - (beta1*sigma[protein_id_missing]), sigma[protein_id_missing]);

  // zeros ~ bernoulli(.05 + ((1-.05)*(1.0 ./ (1.0 + exp(-beta0 + (beta1*obs))))));
  // ones ~ bernoulli(.05 + ((1-.05)*(1.0 ./ (1.0 + exp(-beta0 + (beta1*obs_mis)))))); //mar + ((1 - mar) *
}
