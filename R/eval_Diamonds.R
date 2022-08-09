#############################################################
library(keras)
library(tensorflow)
library(tfprobability)
library(ggplot2)
library(rstan)
library(loo)
library(rstan)
source('R/eval_utils_multi.R')

#Loading of cached data
reps=5
dir_name = 'R/runs/Diamonds'
dir_name = 'R/runs/cpu_DIAMONDS_F1F2_Epo_100000_M_50_T_10/'
df2 = make_plots_and_stats(dir_name = dir_name, reps=reps)


if (FALSE){
  # we used provided MCMC samples
  library(posteriordb)
  library(posterior)
  my_pdb <- pdb_local(path='~/Documents/workspace/posteriordb')
  po <- posteriordb::posterior("diamonds-diamonds", my_pdb)
  draws = get_data(po)
  reference_posterior_draws_info(po)
  posteriordb::stan_code(po)
  rpd <- reference_posterior_draws(po)
  posterior::summarize_draws(rpd)
  #https://mc-stan.org/cmdstanr/reference/fit-method-draws.html
  chef_mcmc = as.data.frame(as_draws_df(rpd))[,1:26]
}
load(file="mcmc/Diamonds/chef_mcmc.rda")
load('data/diamonds.rda')

N = data$N # number data points
P = as.integer(data$K - 1L) #K contains the intercept
Xc = scale(data$X[,-1], center = TRUE, scale = FALSE)
Y = data$Y
ml = lm(Y ~ Xc, data.frame(Y = data$Y, Xc = Xc))
summary(ml)

mean_and_check_mc_error = function(a, atol=0.01, rtol=0.0){
  m = mean(a)
  s = sd(a)/sqrt(length(a))
  if (s > rtol*abs(m) + atol){
    print('There is something foule in the state of denmark')
  }
  return (m)
}

load('R/runs/cpu_DIAMONDS_F1F2_Epo_100000_M_50_T_10/samples_1.rda')
### Intercept
w = samples$w
hist(w[,1], freq = FALSE, 50, xlab='b', col='skyblue', main='VI-Samples for intercept b')
lines(density(posts_mcmc$b), col='red', lwd=2)
lines(density(chef_mcmc$Intercept), col='pink', lwd=2)
#hist(chef_mcmc$Intercept, xlim=c(5, 10))
#lines(density(w[,1]), col='skyblue', lwd=2)
abline(v=ml$coefficients[1], col='darkgreen', lwd=3)
legend('topleft', legend = c('VI', 'MCMC', 'ML'), lty=c(1,1),
       col=c('skyblue', 'red', 'darkgreen'), lwd=3, bty='n')

car::qqPlot(chef_mcmc$Intercept, main='MCMC Intercept')
car::qqPlot(chef_mcmc$sigma, main='MCMC Sigma')
for (i in 1:14){
  car::qqPlot(chef_mcmc[,i], main=paste0('MCMC b ', i))
}
cor(chef_mcmc)

sigma = tf$math$softplus(w[,P+2])$numpy()
sub_sample = 5000
pdf('figures/diamonds.pdf')
par(mfrow = c(2,1))
w_vi = as.data.frame(w[,1:25])
w_vi$sigma = sigma
names(w_vi) = c('Intercept', paste0('b',1:24), 'Sigma')
boxplot(w_vi[1:sub_sample,], ylim=c(-10,15), main='BF-VI', las=2, pch=21)
points(ml$coefficients[1:25], col='red', pch=4)
d = cbind(chef_mcmc$Intercept, chef_mcmc[1:24], chef_mcmc$sigma)
names(d) = c('Intercept', paste0('b',1:24), 'Sigma')
boxplot(d[1:sub_sample,], ylim=c(-10,15), main='MCMC', las=2)
points(ml$coefficients[1:25], col='red', pch=4)
par(mfrow = c(1,1))
dev.off()


