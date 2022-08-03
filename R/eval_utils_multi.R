get_losses_metrics = function(dir_name_r, reps, num_boot = 1000, log_ratios=NULL, load_data=TRUE){
  loss_hists = all_samples = {}
  KL_est = KL_var = ELBO_est = ELBO_var = k_hat_var = upper_slope_loss = lower_slope_loss = evidence_est = evid_var = losses_sd = losses_mean = k_hats = rep(NA, reps)
  k_CI = matrix(rep(NA, 3*reps), ncol = 3)
  evid_CI = matrix(rep(NA, 3*reps), ncol = 3)
  ELBO_CI = matrix(rep(NA, 3*reps), ncol = 3)
  KL_CI = matrix(rep(NA, 3*reps), ncol = 3)
  for (rep in 1:reps) {
    if (load_data) {
      load(file.path(dir_name_r, paste0('loss_hist_',rep,'.rda')))
      loss_hists[[rep]] = loss_hist   
      load(file.path(dir_name_r, paste0('samples_',rep,'.rda')))
      all_samples[[rep]] = samples  
      
      
      xd = 80000:100000
      yd = loss_hists[[rep]][xd]
      lower_slope_loss[rep] = confint(lm(yd ~ xd))[2,1]
      upper_slope_loss[rep] = confint(lm(yd ~ xd))[2,2]
      
      N = length(loss_hist)
      sd(loss_hist[(N-1e3):N]) 
      losses_mean[rep] = mean(loss_hist[(N-1e3):N])
      losses_sd[rep] = sd(loss_hist[(N-1e3):N])
      
      #Estimation of k_hat
      log_joint = samples$LLs + samples$L_priors
      log_ratios = as.vector(log_joint - samples$log_qs)
      
      
    } 
    res_psis = psis(log_ratios)
    k_hats[rep] = res_psis$diagnostics$pareto_k
    
    #See e.g. https://bochang.me/blog/posts/iwae/ 
    evidence_est[rep] = pomp::logmeanexp(log_ratios)
   
    ELBO_est[rep] = mean(log_ratios)
    
     
    #Bootstrap CI
    boot_diff = boot_evid = boot_elbo = boot_khat = rep(NA, num_boot)
    print(rep)
    for (i in 1:num_boot){
      idx = sample(1:length(log_ratios), replace = TRUE)
      boot_khat[i] = psis(log_ratios[idx])$diagnostics$pareto_k
      boot_evid[i] = pomp::logmeanexp(log_ratios[idx])
      boot_elbo[i] =  mean(log_ratios[idx])
      boot_diff[i] =  boot_evid[i] - boot_elbo[i] 
    }
    k_CI[rep,] = quantile(boot_khat, c(0.05, 0.5, 0.95))
    k_hat_var[rep] = var(boot_khat)
    
    ELBO_CI[rep,] = quantile(boot_elbo, c(0.05, 0.5, 0.95))
    ELBO_var[rep] = var(boot_elbo)
    
    evid_CI[rep,] = quantile(boot_evid, c(0.05, 0.5, 0.95))
    evid_var[rep] = var(boot_evid)
    
    KL_CI[rep,] = quantile(boot_diff, c(0.05, 0.5, 0.95))
    KL_var[rep] = var(boot_diff)
    #plot(density(boot_elbo), xlim=c(-28.7,-28.5), main=rep)
    #lines(density(boot_evid), col='red')
    
    
    
    
  }#Reps
  df = data.frame(
    k_hats=k_hats,
    k_hat_q05  =  k_CI[,1],
    k_hat_q5  =  k_CI[,2],
    k_hat_q95  =  k_CI[,3],
    k_hat_var = k_hat_var,
    losses_mean = losses_mean,
    losses_sd = losses_sd,
    ELBO_est = ELBO_est,
    ELBO_q05 = ELBO_CI[,1],
    ELBO_q5 = ELBO_CI[,2],
    ELBO_q95 = ELBO_CI[,3],
    ELBO_var = ELBO_var,
    evidence_est = evidence_est,
    evid_q05 = evid_CI[,1],
    evid_q5 = evid_CI[,2],
    evid_q95 = evid_CI[,3],
    evid_var = evid_var,
    KL_est = KL_est,
    KL_q05 = KL_CI[,1],
    KL_q5 = KL_CI[,2],
    KL_q95 = KL_CI[,3],
    KL_var = KL_var,
    lower_slope_loss = lower_slope_loss,
    upper_slope_loss = upper_slope_loss
  )
  df$reps = 1:reps
  return(list(loss_hists=loss_hists, df=df))
}


make_plots_and_stats = function(dir_name, reps){
  print(paste0("Creating plots in ", dir_name))
  rets = get_losses_metrics(dir_name, reps=reps, num_boot = 1000)
  loss_hists = rets[['loss_hists']]
  df = rets[['df']]
  
  ## Plot of the loss function
  N_Loss = length(loss_hists[[1]])
  idx = seq(1, N_Loss, length.out = 1000) 
  l = quantile(loss_hists[[1]][idx], 0)
  u = quantile(loss_hists[[1]][idx], 0.85)
  #plot(loss_hists[[1]][idx], type='l', main=dir_name, col='skyblue', ylim=c(25,45))
  pdf(file.path(dir_name, 'loss.pdf'))
  plot(idx, loss_hists[[1]][idx], type='l', 
       main=dir_name, col='skyblue', ylim=c(l - 0.05*abs(l),u),
       xlab = 'Epochs',
       ylab='loss (ELBO)',
       cex.main=0.65
  )
  abline(h=-df$evidence_est[1],col='skyblue')
  for (rep in 2:reps){
    #lines(loss_hists[[rep]][idx], col=rep)  
    lines(idx,loss_hists[[rep]][idx], col=rep)  
    abline(h=-df$evidence_est[rep],col=rep)
  }
  dev.off()
  
  ## k-hat plot
  k_bar = mean(df$k_hats)
  var_boot = mean(df$k_hat_var)
  var_runs = var(df$k_hats)
  var_rubin = var_boot + (1 + 1/reps) * var_runs
  k_rubin_lower = k_bar - qnorm(0.95)*sqrt(var_rubin)
  k_rubin_upper = k_bar + qnorm(0.95)*sqrt(var_rubin)
  
  df2 = data.frame(x=0, ymin=k_rubin_lower, ymax=k_rubin_upper, y=k_bar)
  p = ggplot(df, aes(x=reps)) + 
    geom_errorbar(aes(y=k_hat_q5, ymin=k_hat_q05, ymax=k_hat_q95),col='red') + geom_point(aes(y=k_hat_q5)) + 
    ylab('k_hat') +
    geom_errorbar(data=df2, mapping=aes(x=x, ymin=ymin, ymax=ymax)) +
    geom_point(data=df2, aes(x=x, y=k_bar)) +
    labs(title = dir_name, 
         subtitle='Bootstrap (red) and Rubins (black) CI (90)')
  p
  ggsave(file.path(dir_name, 'k_hat.pdf'), p) 
  
  ggplot(df, aes(x=reps, y=evid_q5)) + 
    geom_errorbar(aes(ymin=evid_q05, ymax=evid_q95),col='red') + geom_point() + 
    ylab('Evidence') +
    labs(title = dir_name, 
         subtitle='Bootstrap CI (90) for 5 Different Runs 1K Bootstrap Samples')
  ggsave(file.path(dir_name, 'evidence.pdf')) 
  
  ggplot(df, aes(x=reps, y=ELBO_q5)) + 
    geom_errorbar(aes(ymin=ELBO_q05, ymax=ELBO_q95),col='red') + geom_point() + 
    ylab('ELBO') +
    labs(title = dir_name, 
         subtitle='Bootstrap CI (90) for 5 Different Runs 1K Bootstrap Samples')
  ggsave(file.path(dir_name, 'elbo.pdf')) 
  
  ggplot(df, aes(x=reps, y=KL_q5)) + 
    geom_errorbar(aes(ymin=KL_q05, ymax=KL_q95),col='red') + geom_point() + 
    ylab('KL(q||post) EVIDENCE - ELBO') +
    labs(title = dir_name, 
         subtitle='Bootstrap CI (90) for 5 Different Runs 1K Bootstrap Samples')
  ggsave(file.path(dir_name, 'kl.pdf'))
  
  return(
    data.frame(
      k_bar,
      k_rubin_lower,
      k_rubin_upper,
      dir_name
    )
  )
} 