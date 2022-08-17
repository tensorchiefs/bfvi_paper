#############################################################
# Creating the 1-D Conjugate Prior Figure
library(ggplot2) 
library(dplyr)
library(magrittr)

df = read.csv('Python/Bernoulli/bernoulli_1D_F1F2.csv.gz')
df %>% select(Method, KL) %>% group_by(Method) %>% unique()

#### TODO #### Finish the plot
#From https://bookdown.org/content/2015/figures.html
apatheme=theme_bw(base_size = 22)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'),
        legend.title=element_blank(),
        legend.position = c(0.25,0.85),
        axis.text.y=element_text(size = 20),
        axis.text.x=element_text(size = 20))

df_p = df %>% 
  #filter(M %in% c(-1,0.5,1,2,10,30,100))  %>% #0.5 is Gauss, 0 is MCMC
  filter(M %in% c(-1,0.5,1,10,30,100))  %>% 
  arrange(desc(M),pi) 
ggplot(df_p) + geom_line(aes(x=pi, y=p, col=Method, linetype=Method),size=1.2) + 
  xlab(expression(pi)) +
  apatheme
  
#Saving
ggsave('figures/conjugate_prior.pdf',width = 7, height=7/sqrt(2))
# ggsave('~/Dropbox/Apps/Overleaf/bernvi/images/conjugate_prior.pdf',width = 7, height=3.8)

if (FALSE){
  #### Additional Figure
  library(readr)
  runs_M <- read_csv("runs/Bernoulli_1D/00_conj_prior_kl_eval_kl_all_random.csv")
  runs_M$X1 = NULL
  colnames(runs_M) = c('rep', 1,2,3,4,5,6,7,8,9,10,15,30,100,300,0.5)
  library(reshape2)
  df = melt(runs_M, id="rep")
  df$M = as.numeric(as.numeric(as.character(df$variable)))
  df$variable = NULL
  df$KL = df$value
  df$value = NULL
  df$type = 'BF-VI'
  df$type[df$M == 0.5] = 'Gauss'
  
  
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
  
  ggplot(df) + 
    geom_boxplot(aes(x=M, y=KL, fill=type, group=cut_interval(x=M, length=0.25)))+
    #geom_point(aes(x=M, y=KL, col=type)) + 
    scale_x_log10() + 
    scale_y_log10() +
    apatheme
} 

