#############################################################
# The center parametrisation (harder for MCMC)
# --- HARDER FOR MCMC -----

library(keras)
library(tensorflow)
library(tfprobability)
library(ggplot2)
library(rstan)
library(loo)
source("R/eval_utils_multi.R")

#Loading of cached data
dir_name = 'R/runs/gpu_8CPprior_fix_8SCHOOLS_CP_F1F2_Epo_100000_M_50_T_10/'
reps = 5
df = make_plots_and_stats(dir_name, reps)
df$k_bar #0.5344975
df$k_rubin_lower #0.1107464
df$k_rubin_upper #0.9582486

load('mcmc/8Schools_CP/mcmc_samples.rda')
mcmc_samples = data.frame(posts_mcmc$mu, posts_mcmc$tau, posts_mcmc$theta)
load('R/runs/gpu_8CPprior_fix_8SCHOOLS_CP_F1F2_Epo_100000_M_50_T_10/samples_1.rda')
ww = samples$w
theta_vi = ww[,3:10]
pdf('figures/8schools_cp.pdf')
par(mfrow=(c(1,2)))
boxplot(posts_mcmc$theta, ylim=c(-20.5,30), main='MCMC', cex.sub=0.5, 
        #sub=dir_name,
        #xlab=expression(theta)
        )
boxplot(as.matrix(theta_vi), main='VI', ylim=c(-20.5,30), 
        #xlab=expression(theta)
        )
par(mfrow=(c(1,1)))
dev.off()


