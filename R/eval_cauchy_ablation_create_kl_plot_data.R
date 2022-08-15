library(readr)
library(ggplot2)
library(data.table)
library(keras)
library(tensorflow)
library(tfprobability)
source('R/bern_utils.R')
source('R/eval_utils.R')
library(readr)
library(dplyr)
library(rstan)

############### Creation of the MCMC Data ########
if (FALSE) {
  y_data = c(1.2083935, -2.7329216,  4.1769943,  1.9710574, -4.2004027, -2.384988)
  m = stan_model(file='mcmc/cauchy/mcmc_cauchy.stan')
  data = list(
    N = 6,
    y = y_data
  )
  fit <- sampling(m, data=data, iter = 8000, chains=8)
  posts_mcmc = extract(fit)
  save(posts_mcmc, file = 'mcmc/cauchy/mcmc_samples_cauchy.rda')
  #fit = pystan.stan(model_code=stan_code, data=coin_dat, iter=1000, chains=1)
}
load('mcmc/cauchy/mcmc_samples_cauchy.rda')
mcmc_df = as.data.frame(posts_mcmc) 
names(mcmc_df) = c("w","lp__")
mcmc_df = mcmc_df %>%  arrange(w) #Sorting

#Normalization
n_const = 1/sum(diff(mcmc_df$w) * exp(mcmc_df$lp__)[-1])
mcmc_df$lp__ = n_const*exp(mcmc_df$lp__) 

N_samples = 5000
mcmc_df = sample_n(mcmc_df,N_samples)

Ms = c(2,3,6,10,20,30,50,60)#,100,200,300)

methods = c('F1F2', 'SigmoidF2', 'TruncF2')
for (method in methods){
dir = "bfvi/runs/cauchy_ablation/"
#Checking if all files exists 
for (M in Ms){
  param_fn = paste0(dir,"cauchy_eval_ablation_M_", M, "_", method, "_params.csv")
  print(file.exists(param_fn))
}

########## Main Loop #######
#We set M=0 to indicate MCMC
df_plot = data.frame(M=0, method='MCMC', x = mcmc_df$w, density = mcmc_df$lp__, seed=1)
df_kl = NULL
for (M in Ms){
  param_fn = paste0(dir,"cauchy_eval_ablation_M_", M, "_", method, "_params.csv")
  print(paste0("Starting with M ", M, "file", param_fn))
  params <- read_csv(param_fn, col_names = TRUE)
  params = params[,-1] 
  M1 = ncol(params) - 5 + 1 #M+1
  colnames(params) = c('seed','epoch','az','bz', paste0('theta_p', 1:M1))
  seeds = unique(params$seed)
  # theta = to_theta(tf$reshape(tf$Variable(as.numeric(params[line,5:ncol(params)])),c(1L,-1L)))
  beta_dists_for_h = init_beta_dist_for_h(M1)
  
  for (seed in seeds){
    line = which(params$seed == seed & params$epoch!='initial')
    theta = to_theta(tf$reshape(tf$Variable(as.numeric(params[line,5:ncol(params)])),c(1L,-1L)))
    num_samples = N_samples
    df_w_logqw = get_w_logqw(theta = theta, num_samples = num_samples,
                             a_np = params$az[line], 
                             b_np = params$bz[line], 
                             trafo = method) %>% 
                    arrange(w)  #Sorting
    
    
    df_plot = rbind(df_plot, data.frame(M = M1 - 1,
      method=method, x = df_w_logqw$w, density = exp(df_w_logqw$log_q_w), seed=seed))
  }
  
  #Calculation of KL Divergence (might take some time)
  if(TRUE){
    for (seed in seeds){
      line = which(params$seed == seed & params$epoch!='initial')
      num_samples = 1e5L
      df_w_logqw = get_w_logqw(theta = theta, num_samples = num_samples,
                               a_np = params$az[line], 
                               b_np = params$bz[line], 
                               trafo = method)
      kl = eval_kl_cauchy(df_w_logqw$w, df_w_logqw$log_q_w)
      df_kl_c = data.frame(seed, kl, M=M1-1, method, num_samples)
      if (is.null(df_kl)) {
        df_kl = df_kl_c
      } else{
        df_kl = rbind(df_kl, df_kl_c)
      }
    }
  }
  
}#M in Ms

if (method=='SigmoidF2'){
  df_kl_SigmoidF2 = df_kl 
  write_csv(df_kl_SigmoidF2, file='R/runs/Cauchy_1D/df_kl_SigmoidF2.csv')
  write_csv(df_plot, file='R/runs/Cauchy_1D/df_plot_SigmoidF2.csv')
}
if (method=='F1F2'){
  df_kl_F1F2 = df_kl 
  write_csv(df_kl_F1F2, file='R/runs/Cauchy_1D/df_kl_F1F2.csv')
  write_csv(df_plot, file='R/runs/Cauchy_1D/df_plot_F1F2.csv')
}
if (method=='TruncF2'){
  df_kl_TruncF2 = df_kl 
  write_csv(df_kl_TruncF2, file='R/runs/Cauchy_1D/df_kl_TruncF2.csv')
  write_csv(df_plot, file='R/runs/Cauchy_1D/df_plot_TruncF2.csv')
}
}

