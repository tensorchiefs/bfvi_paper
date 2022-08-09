#############################################################
library(keras)
library(tensorflow)
library(tfprobability)
library(ggplot2)
library(rstan)
library(loo)
library(readr)
source('R/eval_utils_multi.R')

dir_name = 'R/runs/gpu_2D_F1F2_Epo_100000_M_50_T_10/'
reps = 5
#Make Standard Plots and save to dir_name
df2 = make_plots_and_stats(dir_name = dir_name, reps=reps)
#0.68 (0.51, 0.86) 90$\perc$ CI

#### Loading of Data and MCMC Samples #######
load('data/2D.rdata')
x_r
y_r
N = 6L # number data points
P = 2L #number of predictors
D = P + 2L  #Dimension (number of weights) predictors + intercept + sigma
#Loading MCMC SAMPLES --> posts_mcmc
load('mcmc/2D/2D_mcmc_samples.rda')
if (FALSE){ #Time consuming MCMC
  dl = list(N=N,x=x_r,y=y_r, P = P)
  samples = stan(file="mcmc/2D/tmvi_multiple_lr_sigma.stan", data=dl, iter = 24000)
  samples
  #traceplot(samples)
  posts_mcmc = extract(samples)
  save(posts_mcmc, file = 'mcmc/2D/2D_mcmc_samples.rda')
}

#### Loading of samples from variational posterior
load('R/runs/gpu_2D_F1F2_Epo_100000_M_50_T_10/samples_1.rda')
w = samples$w



####### 
# For paper
sigma = tf$math$softplus(w[,P+2])$numpy()
TT = 1000 #Subsampling to avoid too large plots
df = data.frame(beta1 = posts_mcmc$w[1:TT,1], 
                beta2=posts_mcmc$w[1:TT,2], 
                sigma=posts_mcmc$sigma[1:TT],
                intercept = posts_mcmc$b[1:TT],
                type='MCMC')
df = rbind(df, data.frame(
  beta1 = w[1:TT,2], 
  beta2=w[1:TT,3], 
  sigma=sigma[1:TT],
  intercept  =w[1:TT,1],
  type='BF-VI'))
library(ggplot2)
library(GGally)
p = ggpairs(df, 
            aes(color = type, alpha=0.5), 
            columns = c(1,2,3,4),
            columnLabels = c('beta1', "beta2", "sigma", 'intercept')
) + 
  scale_color_manual(values = c('blue', 'red')) + 
  scale_fill_manual(values = c('blue', 'red')) 
  #labs(title = dir_name)
p = p + theme_bw() + theme(
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank()
)
p
ggsave(p, filename = 'figures/2D_pairs.pdf', width = 7, height = 7)
#Attention Manual relabeling has been done
#ggsave(p, filename = '~/Dropbox/Apps/Overleaf/bernvi/images/2D_pairs.pdf',width = 7, height=7)

###### Pairs Plot for the appendix ########
#########  Mean Field Gaussian ########
# Mean field gaussian
# See python NB 2D_MF_Gaussian_VI.ipynb

##### Calculation of k_hat
df3 = NULL
for (i in 1:5){
  log_qs = read_csv(paste0("R/runs/2D_MF_Gauss/2D_MFG_qs_samples_",i,".csv.gz"))
  log_joints = read_csv(paste0("R/runs/2D_MF_Gauss/2D_MFG_log_joints_samples_",i,".csv.gz"))
  log_ratios = as.matrix(log_joints - log_qs, ncols=1)
  df3t = get_losses_metrics(dir_name_r = dir_name, reps=1, num_boot = 1000, 
                     log_ratios = log_ratios, load_data = FALSE)[['df']]
  if (i == 1){
    df3 = df3t
  } else{
    df3 = rbind(df3,df3t)
  }
}
save(df3, file='R/runs/2D_MF_Gauss/metrics_reps.rda')


b = read_csv("R/runs/2D_MF_Gauss/2D_MFG_intercept_samples.csv.gz")
w = read_csv("R/runs/2D_MF_Gauss/2D_MFG_w_samples.csv.gz")
s = read_csv("R/runs/2D_MF_Gauss/2D_MFG_sigma_samples.csv.gz")
TT = 1000
df = data.frame(beta1 = posts_mcmc$w[1:TT,1], 
                beta2=posts_mcmc$w[1:TT,2], 
                sigma=posts_mcmc$sigma[1:TT],
                intercept = posts_mcmc$b[1:TT],
                type='MCMC')

df = rbind(df, data.frame(
  beta1 = w$w1[1:TT], 
  beta2 = w$w2[1:TT], 
  sigma= s$b[1:TT],
  intercept  = b$b[1:TT],
  type='Gauss-MF'))
library(ggplot2)
library(GGally)
p = ggpairs(df, 
            aes(color = type, alpha=0.5), 
            columns = c(1,2,3,4),
            columnLabels = c('beta1', "beta2", "sigma", 'intercept')
) + 
  scale_color_manual(values = c('blue', 'red')) + 
  scale_fill_manual(values = c('blue', 'red'))
p = p + theme_bw() + theme(
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank()
)
p
ggsave(p, filename = 'figures/2D.pairs.gauss_vs.mcmc.pdf', width = 7, height = 7)
#ggsave(p, filename = '~/Dropbox/Apps/Overleaf/bernvi/images/2D.pairs.gauss_vs.mcmc.pdf',width = 7, height=7)

  





