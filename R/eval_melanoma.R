#############################################################
# Melanoma experiment
library(keras)
library(tensorflow)
library(tfprobability)
library(ggplot2)
library(rstan)
library(loo)
library(rstan)

source("R/eval_utils_multi.R")

dir_name = 'R/runs/cpu_MELA_F1F2_Epo_100000_M_50_T_10/'
reps = 5
df = make_plots_and_stats(dir_name, reps)
df$k_bar #.04505376
df$k_rubin_lower #-0.2926889
df$k_rubin_upper #0.3827964

# load variational posterior samples from BF-VI M2 fit
load('R/runs/cpu_MELA_F1F2_Epo_100000_M_50_T_10/samples_1.rda')
w = samples$w

#From https://bookdown.org/content/2015/figures.html
apatheme=theme_bw(base_size = 22)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'),
        legend.title=element_blank(),
        legend.position = c(0.34,0.8),
        axis.text.y=element_text(size = 20),
        axis.text.x=element_text(size = 20))

#### Evaluation of the test-set
library(data.table)
testdata = fread("data/mela_tab_testData.csv")

x = testdata$standAge
y_obs = testdata$target

lls = rep(NA, length(x))
ps = rep(NA, length(x))
for (i in 1:length(x)) {
  h = w[,1] * x[i] + w[,2]
  p = 1./(1+exp(-h))
  p_mean = mean(p)
  lls[i] = y_obs[i] * log(p_mean) + (1 - y_obs[i]) * log(1 - p_mean)
  ps[i] = p_mean
}
hist(lls)
mean(lls) #-0.08501698
library(pROC)
ci.auc(y_obs, ps) #95% CI: 0.6086-0.7075 (DeLong)
auc(y_obs, ps) #0.66

#### MCMC
library(data.table)
if (FALSE){
  # lineare regression with 1 predictor and 4 data points
  trainData <- read.csv("data/mela_tab_data.csv")
  N = nrow(trainData)
  P = 1L #number of predictors
  D = P + 1L  #Dimension (number of weights) predictors + intercept + sigma
  melanom_data = list(N=N, M=P, X=matrix(trainData$standAge, ncol=1), y=trainData$target)
  m = stan_model(file='mcmc/mela_M2/mcmc_mela_m2.stan')
  fit <- sampling(m, data=melanom_data, iter=10000)
  posts_mcmc = extract(fit)
  fwrite(posts_mcmc, "mcmc/mela_M2/mcmc_M2.csv.gz")
}
mcmc_samples = fread("mcmc/mela_M2/mcmc_M2.csv.gz")
d = density(mcmc_samples$beta, adjust = 1.5) #2. larger bw for smoother plot
df = data.frame(method='M2: MCMC', slope = d$x, density = d$y)

##### BFVI
w1 = NULL
for (b in 1:11){
 w1 = append(w1, samples[(b-1)*4+4]$w[,1])
}
d = density(w1, adjust = 1) #2. larger bw for smoother plot
df = rbind(df, data.frame(method='M2: BF-VI', slope = d$x, density = d$y))

##### samples from semi-structured M3 model from Ivonne
tmvi_cnn = fread("Ivonne_MA/semi_posterior_slope_age.csv")
#plot(tmvi_cnn$slope_weight, tmvi_cnn$slope_density)
idx = sort(tmvi_cnn$slope_weight, index.return = TRUE)$ix
#plot(tmvi_cnn$slope_weight, tmvi_cnn$slope_density)
df = rbind(df, data.frame(method='M3: CNN+BF-VI', slope = c(tmvi_cnn$slope_weight[idx], 0.95), density = c(tmvi_cnn$slope_density[idx], 0)))

ggplot(df) +
  geom_line(aes(x=slope, y=density, col=method, linetype=method),size=1.5) +
  #geom_vline(xintercept=0.7057) + #ML solution
  xlim(-0.25,1) +
  scale_color_manual(values = c('blue', 'red', 'green')) +
  scale_linetype_manual(values = c('solid', 'dotted', 'solid')) +
  #scale_size_manual(values = 50*c(1.2,2,1.2)) +
  xlab(expression(beta[1])) +
  #xlab('slope [age]') +
  apatheme

ggsave('figures/age-slope.pdf',width = 7, height=4)
#ggsave('~/Dropbox/Apps/Overleaf/bernvi/images/age-slope.pdf',width = 7, height=4)
