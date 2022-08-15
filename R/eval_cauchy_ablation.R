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

df_plot = read_csv('R/runs/Cauchy_1D/df_plot_F1F2.csv')
###########  produce cauchy_F1F2_M_comparison plot ###########
apatheme=theme_bw(base_size = 22)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'),
        #legend.title=element_blank(),
        legend.position = c(0.2,0.75),
        axis.text.y=element_text(size = 20),
        axis.text.x=element_text(size = 20))

d =  fread("Python/Cauchy/Gauss-VI_densities.csv.gz")[-1,]
kl_gauss = eval_kl_cauchy(d$V1, log(d$V2))
df_p = rbind(df_plot, data.frame(M=0.5, method='Gauss-VI', x=d$V1, density=d$V2, seed=1))
dd = df_p %>% filter(method %in% c('F1F2', 'Gauss-VI', 'MCMC')) %>%
  filter(M %in% c(0.5,2,6,10,30,50,0)) %>% #0.5 is Gauss, 0 is MCMC
  #group_by(method) %>%
  #mutate(num_rows=n()) %>%
  #sample_frac(0.4, weight=num_rows) %>%
  #ungroup %>%
  filter(as.numeric(seed) < 11) %>% arrange(method, factor(M), x)

library(dplyr)
p = ggplot(dd) +
  geom_line(aes(x=x, y=density, col=factor(M)),size=0.75) +
  xlab(expression(xi)) +
  scale_color_manual(
    name = 'Method',
    values = c(
      "0.5" = "green",
      "2" = "darkgreen",
      "6" = "steelblue",
      "10" = "pink",
      "30" = "grey",
      "50" = "darkblue",
      "0" = "red"),
    labels = c('Gauss-VI', 'M=2','M=6','M=10','M=30','M=50', 'MCMC')
  ) +
  scale_linetype(guide="none") +
  apatheme
p
ggsave('figures/cauchy_F1F2_M_comparison.pdf', p,width = 7, height=4.8)
#ggsave(p, width = 7, height=4.8,filename = '~/Dropbox/Apps/Overleaf/bernvi/images/cauchy_F1F2_M_comparison.pdf')

########## produce kl_vs_M_methods_chauchy plot #######
df_kl = read_csv('R/runs/Cauchy_1D/df_kl_F1F2.csv')
df_kl = rbind(df_kl, read_csv('R/runs/Cauchy_1D/df_kl_SigmoidF2.csv'))
df_kl = rbind(df_kl, read_csv('R/runs/Cauchy_1D/df_kl_TruncF2.csv'))

####### KL - Dependence on M #####
apatheme=theme_bw(base_size = 22)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'),
        legend.title=element_blank(),
        legend.position = c(0.8,0.8),
        axis.text.y=element_text(size = 20),
        axis.text.x=element_text(size = 20))



## Adding Gaus #####
df_kl = rbind(df_kl,
              data.frame(seed = 1, kl=kl_gauss, M=0.5, method='Gauss-VI', num_samples=1))

dat = df_kl %>%
  filter(M < 61009)
  #filter(method == 'F1F2' | method == 'Gauss-VI')

gg_chauchy = ggplot(dat) +
  geom_jitter(aes(x=M, y=kl, col=method),alpha=0.9) +
  #geom_boxplot(aes(x=M, y=kl, fill=method, group=cut_interval(x=M, length=0.5))) +
  #geom_line(data=df2, aes(x=x,y=y), linetype=2, size=1.2) +
  scale_x_log10() +
  ylab(expression(paste('KL(q(',theta,') || p(',theta,'|D))')))+
  #scale_y_log10(name=expression(paste('KL(q(',theta,') || p(',theta,'|D))'))) +
  #annotate("text", x = 20, y = 2/20, label = "~1/M", size=6) +
  scale_color_manual(
    name = 'Method',
    values = c(
      "Gauss-VI" = "green",
      "F1F2" = "darkblue",
      "SigmoidF2" = "turquoise",
      "TruncF2" = "brown"),
    labels = c('Gauss-VI', 'Gauss-F1F2','Gauss-SigmoidF2','TruncGauss-F2')
  )  +
  apatheme
gg_chauchy
ggsave('figures/kl_vs_M_methods_chauchy.pdf',gg_chauchy)
#ggsave(gg_chauchy, width = 7, height=4.8, filename = '~/Dropbox/Apps/Overleaf/bernvi/images/kl_vs_M_methods_chauchy.pdf')


