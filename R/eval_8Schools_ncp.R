#############################################################
# The non-central parametrisation (easier for MCMC)
# --- EASIER FOR MCMC -----
library(keras)
library(tensorflow)
library(tfprobability)
library(tidyverse)
library(rstan)
library(loo)
source("R/eval_utils_multi.R")

#Loading of cached data
dir_name = 'R/runs/gpu_8SCHOOLS_F1F2_Epo_100000_M_50_T_10/'
reps = 5
df = make_plots_and_stats(dir_name, reps)
df$k_bar #0.3633735
df$k_rubin_lower #0.17110
df$k_rubin_upper #0.5556

load('mcmc/8Schools_CP/mcmc_samples.rda')
mcmc_samples = data.frame(posts_mcmc$mu, posts_mcmc$tau, posts_mcmc$theta)
if (FALSE){
  #For the MCMC we use the easier parametrisation
  m = stan_model(file='eight_schools_ncp.stan')
  fit <- sampling(m, data=data, iter = 10000, seed=194838,control=list(adapt_delta=0.9))
  posts_mcmc = extract(fit)
  save(posts_mcmc, file = '~/Dropbox/__ZHAW/__Projekte_Post_ZHAH/shared_Oliver_Beate/tmvi/runs/8SCHOOLS_NCP_Run1/mcmc_samples.rda')

  #Using Stan AUVI
  fit_vb=vb(m,  data=data, iter=1e6,output_samples=1e5,tol_rel_obj=0.001)
  vb_sample=extract(fit_vb)
  vb_samples = data.frame(vb_sample$mu, vb_sample$tau, vb_sample$theta)
  error_mean_sigma(mcmc_samples, vb_samples)
  #L2_mu   L2_sig
  #1.892609 1.207052

  #Second run
  #L2_mu   L2_sig
  #1.985089 1.227350
  psis(vb_sample)
}

load('R/runs/gpu_8SCHOOLS_F1F2_Epo_100000_M_50_T_10/samples_1.rda')
wws = samples$w
theta_tilde_vi = wws[,3:10]
mu = wws[,2]
tau = exp(0.1 * wws[,1])
theta = mu + tau*theta_tilde_vi

pdf('figures/8schools_ncp.pdf')
par(mfrow=(c(1,2)))
boxplot(as.matrix(posts_mcmc$theta_tilde), main='MCMC',
        #sub=dir_name,ylim=c(-4.5,4.5),
        #ylab='theta_tilde'
        )
boxplot(as.matrix(theta_tilde_vi), main='VI',
        #sub=dir_name,
        ylim=c(-4.5,4.5)
        #,ylab='theta_tilde'
        )
par(mfrow=(c(1,1)))
dev.off()



