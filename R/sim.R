#' 
#' library(tidyverse)
#' library(fastDummies)
#' 
#' # Questions: 
#' # - Sample variances from prior dist?
#' # - Independent variance terms?
#' 
#' #' Sim data
#' #' 
#' #' @param N_proteins number of proteins in simulation
#' #' @param N_runs number of runs in sim
#' #' @param N_features number of features in sim
#' simulate_data = function(N_proteins=1, N_runs=4, N_features=8) {
#'     
#'     run_mu = rnorm(4, 20, 5)
#'     run_std = rlnorm(4, 0, .5)
#'     run = rnorm(N_runs, run_mu, run_std)
#'     
#'     # feature_mu = rnorm(1, 0., 2.5)
#'     # feature_std = rlnorm(1, 0, .25)
#'     feature = rnorm(N_features, 0, 2.5)
#'     
#'     sigma_mu = rnorm(4, 0, .5)
#'     sigma_std = rlnorm(4, 0, .25)
#'     sigma = rlnorm(N_runs, 0., .5)
#'     
#'     run_id = rep(seq_len(N_runs), N_features)
#'     feature_id = rep(seq_len(N_features), each=N_runs)
#'     
#'     obs = rnorm(n=N_runs*N_features, 
#'                 run[run_id] + feature[feature_id], sigma[run_id])
#'     
#'     data = data.frame(run=run_id, feature=feature_id, intensity=obs)
#'     # samples = list(run=run, feature=feature,sigma=sigma)
#'     samples = list(run_mu=run_mu, run_std=run_std, run=run,
#'                    # feature_mu=feature_mu, feature_std=feature_std,feature=feature, sigma_mu=sigma_mu, sigma_std=sigma_std,
#'                    sigma=sigma)
#'     return(list(data, samples))
#' }
#' 
#' sim = simulate_data(N_features=10)
#' test_data=sim[[1]]
#' samples=sim[[2]]
#' samples
#' test_data %>% ggplot() + 
#'     geom_point(aes(x=run, y=intensity, color=as.character(feature))) + 
#'     geom_line(aes(x=run, y=intensity, color=as.character(feature)))
#' 
#' stan_file="/Users/kohler.d/Documents/Northeastern/msstats/MSstats/inst/stan/bayes_model.stan"
#' stan_input = list(N_obs=nrow(test_data), 
#'                   R=4,
#'                   Feat=10,
#'                   P=1,
#'                   obs=test_data[,3],
#'                   run_id=test_data[,1],
#'                   feature_id=test_data[,2],
#'                   protein_id=rep(1, nrow(test_data)))
#' 
#' fit = stan(file = stan_file, data = stan_input,
#'            chains = 4, iter = 20000, cores = 4,
#'            seed=100, model_name="MSstats_model")
#' 
#' 
#' ## TODO: Automate this to test different numbers of runs/features
#' simulate_data = function(N_proteins=1, N_runs=3, N_features=4) {
#'     
#'     ## Generate data
#'     run_mu = runif(N_runs, 10, 20)
#'     run_std = runif(N_runs, 1, 3)
#'     run = rnorm(N_features*N_runs*N_runs, run_mu, run_std)
#'     
#'     run = matrix(run, nrow=N_features*N_runs, byrow=TRUE)
#'     
#'     feature_mu = runif(N_features, -3, 3)
#'     feature_std = runif(N_features, .25, 1.5)
#'     feature = rnorm(N_features*N_runs*N_features, feature_mu, feature_std)
#'     feature = matrix(feature, nrow=N_features*N_runs, byrow=TRUE)
#'     
#'     error = rnorm(N_features*N_runs, 0, 1.)
#'     
#'     ## Make dummy matrix
#'     X = dummy_cols(
#'         data.frame("Run" = rep(paste0("run_", seq_len(N_runs)), each=N_features),
#'                "Feature" = rep(paste0("feature_", seq_len(N_features)), N_runs))) %>% 
#'         select(-Run, -Feature) %>% as.matrix()
#'     
#'     B = cbind(run, feature)
#'     obs = rowSums(X*B) + error
#'     one_hot_X = X[,c(2:N_runs, (N_runs+2):(N_runs+N_features))]
#'     one_hot_X = cbind(rep(1, nrow(one_hot_X)),one_hot_X)
#'     
#'     data = list(observations = obs, one_hot_X=one_hot_X,
#'                    sample_info=list(
#'         run_mu=run_mu, run_std=run_std,
#'         feature_mu=feature_mu, feature_std=feature_std, error=error))
#'     return(data)
#' }
#' 
#' input_data = simulate_data(N_proteins=1, N_runs=3, N_features=4)
#' 
#' stan_file="/Users/kohler.d/Documents/Northeastern/msstats/MSstats/inst/stan/test_model.stan"
#' 
#' stan_input = list(N = length(input_data$observations),
#'                   P = ncol(input_data$one_hot_X),
#'                   X=input_data$one_hot_X,
#'                   y=input_data$observations
#' )
#' 
#' fit = stan(file = stan_file, data = stan_input,
#'            chains = 4, iter = 200000, cores = 4,
#'            seed=100, model_name="MSstats_model")
#' 
#' fit
#' 
#' lm(obs~test)

