data {
  int<lower=0> N_obs;
  // int<lower=0> N_missing;
  int<lower=0> R;
  int<lower=0> Feat;
  int<lower=0> P;

  vector[N_obs] obs;
  int run_id[N_obs];
  int feature_id[N_obs];
  int protein_id[N_obs]; 
  
  // int run_id_missing[N_missing];
  // int feature_id_missing[N_missing];
  // int protein_id_missing[N_missing];

  // int zeros[N_obs];
  // int ones[N_missing];
  
  vector<lower=0, upper=40>[R] run_mu_prior;
  // vector<lower=0>[R] sigma_prior;
  
  // vector<lower=0>[Feat] feat_mu_prior;
  // vector<lower=0>[Feat] sigma_feat_prior;
  
  // vector<lower=0>[R] sigma_prior;
  // real sigma_prior;
  
  // real beta0;
  // real beta1;
}
parameters {
  // real beta0;
  // real beta1;
  // real mar;
  vector<lower=0, upper=40>[R] run_mu;
  // vector<lower=0> run_std;
  real<lower=0> run_std;
  vector<lower=0, upper=40>[R] run;
  
  vector<lower=-10, upper=10>[Feat] feature_mu;
  // vector<lower=0> feature_std;
  real<lower=0> feature_std;
  vector<lower=-20, upper=20>[Feat] feature;
  // vector<lower=0, upper=5>[P] sigma;
  // vector<lower=0, upper=50>[N_missing] obs_mis;
  // vector<lower=0, upper=5>[N_missing] missing_sigma;
  
  vector<lower=0>[R] sigma;
  // real sigma;
}
model {
  // beta0 ~ normal(8., 2.);
  // beta1 ~ normal(.8, .2);
  // mar ~ beta(1, 20);
  
  // Run effect
  run_mu ~ normal(run_mu_prior, 2.);
  run_std ~ lognormal(0, 1);
  run ~ normal(run_mu, run_std);
  
  // Feature effect
  feature_mu ~ normal(0., 2.5);
  feature_std ~ lognormal(0, 1);
  feature ~ normal(feature_mu, feature_std);
  
  sigma ~ lognormal(0, 1);
  
  obs ~ normal(run[run_id] + feature[feature_id], sigma[run_id]);
  // obs_mis ~ normal(run_mu[run_id_missing] + feature_mu[feature_id_missing], sigma[run_id_missing]);
   // - (beta1*sigma[protein_id_missing])

  // zeros ~ bernoulli(.05 + ((1-.05)*(1.0 ./ (1.0 + exp(-beta0 + (beta1*obs))))));
  // ones ~ bernoulli(.05 + ((1-.05)*(1.0 ./ (1.0 + exp(-beta0 + (beta1*obs_mis)))))); //mar + ((1 - mar) *
}
