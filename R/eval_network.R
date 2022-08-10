#############################################################
# Creating the 1-D NN Regression figure

library(keras)
library(tensorflow)
library(tfprobability)
library(ggplot2)
library(rstan)
library(loo)
library(rstan)

#source('eval_utils.R')
source("R/eval_utils_multi.R")
dir_name = 'R/runs/gpu_NETWORK_F1F2_Epo_100000_M_50_T_10/'
reps = 5
df = make_plots_and_stats(dir_name, reps)
df$k_bar #0.8665519
df$k_rubin_lower #0.6445761
df$k_rubin_upper #1.088528


#Reading of the MCMC samples from the python script
#https://github.com/tensorchiefs/bfvi_paper/blob/main/Python/networks.ipynb
post_mcmc = read.csv('Python/network/mcmc_y_quantiles.csv')
post_mcmc$X = NULL
names(post_mcmc) = c('X','Mean', 'qlow', 'qhigh', 'med')

#Loading of cached data
#y_preds_MFTM = read.csv(file.path(dir_name_mcmc,'y_preds_VIMLTS.csv'), header=FALSE)
xx = read.csv('Python/network/xx.csv', header=FALSE)
dim(xx)
xx = xx[,1]
y_med_g = read.csv('Python/network/y_predictions_50_VIGAUSS.csv', header=FALSE)
y_l_g = read.csv('Python/network/y_predictions_05_VIGAUSS.csv', header=FALSE)
y_u_g = read.csv('Python/network/y_predictions_95_VIGAUSS.csv', header=FALSE)

#load(file.path(dir_name, 'loss_hist_1.rda'))
load('R/runs/gpu_NETWORK_F1F2_Epo_100000_M_50_T_10/samples_1.rda')
P = 1L
x_r = c(-5.41239, -4.142973,-5.100401,-4.588446,-2.0570369,-2.0035906,-3.7404475,
        -5.344997,4.506781,5.9761415,4.073539,5.168227,4.1196156,2.4791312,2.0348845, 2.7495284)
y_r = c(0.973122  ,  0.96644104,  1.2311585 ,  0.5193988 , -1.2059958 ,
        -0.9434611 ,  0.8041748 ,  0.82996416, -1.3704962 , -0.3733918 ,
        -0.98566836, -1.1550032 , -1.0276004 ,  0.539029  ,  1.5336514 ,
        0.34641847)
x_r = x_r[1:9]
y_r = y_r[1:9]
N = length(x_r) # number data points
plot(x_r, y_r)
N = length(x_r)
P = 1L
D = 10L
sigma=0.2

mean_and_check_mc_error = function(a, atol=0.01, rtol=0.0){
  m = mean(a)
  s = sd(a)/sqrt(length(a))
  if (s > rtol*abs(m) + atol){
    print('There is something foule in the state of denmark')
  }
  return (m)
}

w = samples$w

# get the log-likelihood and the prior
get_post_predictive <- function(w, x) {
  w = tf$Variable(w, dtype='float32')
  x = tf$Variable(x, dtype='float32')
  #Determination of the likelihood
  w_11 =  tf$slice(w,c(0L,0L),c(T,1L))  
  b_11 =  tf$slice(w,c(0L,1L),c(T,1L))   
  w_12 =  tf$slice(w,c(0L,2L),c(T,1L))  
  b_12 =  tf$slice(w,c(0L,3L),c(T,1L))
  w_13 =  tf$slice(w,c(0L,4L),c(T,1L))  
  b_13 =  tf$slice(w,c(0L,5L),c(T,1L))
  #See https://stackoverflow.com/questions/33858021/outer-product-in-tensorflow
  x_rep = tf$reshape(x, shape=c(N,1L))
  
  w_11_rep = tf$reshape(w_11, shape=c(1L, T))
  h_1 = tf$sigmoid(x_rep * w_11_rep + tf$reshape(b_11, shape=c(1L,T)))
  w_12_rep = tf$reshape(w_12, shape=c(1L, T))
  h_2 = tf$sigmoid(x_rep * w_12_rep + tf$reshape(b_12, shape=c(1L,T)))
  w_13_rep = tf$reshape(w_13, shape=c(1L, T))
  h_3 = tf$sigmoid(x_rep * w_13_rep + tf$reshape(b_13, shape=c(1L,T)))
  
  w_21 =  tf$slice(w,c(0L,6L),c(T,1L))  
  w_22 =  tf$slice(w,c(0L,7L),c(T,1L))  
  w_23 =  tf$slice(w,c(0L,8L),c(T,1L))  
  b_2 =  tf$slice(w,c(0L,9L),c(T,1L))
  
  mu = h_1 * tf$squeeze(w_21)  + h_2 * tf$squeeze(w_22) + h_3 * tf$squeeze(w_23) + tf$squeeze(b_2)
  mu = tf$transpose(mu) #oups did it wrong way
  return (tfp$distributions$Normal(loc=mu, scale=sigma)$sample())
}





### Intercept
N = 100L
T = 50000L
x = seq(-10,10, length.out=N)
y = get_post_predictive(w, x)
#y = get_post_predictive_6(w, x)
y_med = apply(y, 2, quantile, probs=0.5)
y_l = apply(y, 2, quantile, probs=0.05)
y_u = apply(y, 2, quantile, probs=1-0.05)

df_samples = data.frame(x = x, y = y[1,]$numpy(), t = 1L)
for (t in 2:20) {
  df_samples = rbind(df_samples, data.frame(x = x, y = y[t,]$numpy(), t = t))
}
df_samples$t = as.factor(df_samples$t)
ggplot(df_samples) + geom_line(aes(x=x,y=y, group=t))

apatheme=theme_bw(base_size = 22)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'),
        legend.title=element_blank(),
        legend.position = c(0.75,0.85),
        axis.text.y=element_text(size = 20),
        axis.text.x=element_text(size = 20))

df_quantiles = data.frame(x = post_mcmc$X, y = post_mcmc$qlow, method = 'MCMC', d=1)
df_quantiles = rbind(df_quantiles, data.frame(x = post_mcmc$X, y = post_mcmc$qhigh, method = 'MCMC', d=2))
df_quantiles = rbind(df_quantiles, data.frame(x = post_mcmc$X, y = post_mcmc$med, method = 'MCMC', d=3))

df_quantiles = rbind(df_quantiles, data.frame(x = x, y = y_med, method = 'BF-VI',d=4))
df_quantiles = rbind(df_quantiles, data.frame(x = x, y = y_l, method = 'BF-VI',d=5))
df_quantiles = rbind(df_quantiles, data.frame(x = x, y = y_u, method = 'BF-VI',d=6))

df_quantiles = rbind(df_quantiles, data.frame(x = xx, y = y_u_g[,1], method = 'Gauss-MF',d=7))
df_quantiles = rbind(df_quantiles, data.frame(x = xx, y = y_l_g[,1], method = 'Gauss-MF',d=8))
df_quantiles = rbind(df_quantiles, data.frame(x = xx, y = y_med_g[,1], method = 'Gauss-MF',d=9))

ggplot(df_quantiles) + 
  geom_line(data = df_samples, aes(x=x,y=y,group=t), col='lightblue', alpha=0.5)+
  geom_point(data = data.frame(x_r, y_r), aes(x=x_r, y=y_r), col='black', size=2) + 
  geom_line(data = df_quantiles, aes(x=x,y=y, col=method, group = d)) +
  scale_color_manual(values = c('blue', 'orange', 'red')) + 
  apatheme 
ggsave('figures/network_ppd.pdf', width = 7, height=3.8)
# ggsave('~/Dropbox/Apps/Overleaf/bernvi/images/network_ppd.pdf', width = 7, height=3.8)

