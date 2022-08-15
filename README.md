# bfvi_paper

# How to Reproduce the results

## Cauchy
### Workflow
#### Raw Data
##### Gauss-VI: Python/Cauchy/Cauchy.ipnb → Python/Cauchy/Gauss-VI_densities.csv.gz
Gauss-VI_densities.csv contain samples w (first colom, in paper called xi) and log_p (second colom)

##### MCMC 
mcmc/cauchy/mcmc_cauchy.stan + R/eval_cauchy_ablation_create_kl_plots.R -> mcmc/cauchy/mcmc_samples_cauchy.rda 

##### BF-VI
run bfvi/bfvi/vimlt_fast.py: define methods and layers
run bfvi/bfvi/cauchy_ablation.py: definitions of functions
Parameter-Files 
Many ablation runs have been done with python (bfvi/bfvi/cauchy_eval_ablation.py) on GPU and CPU enviromemnts. The results are stored the from of the parameters defining the Bernstein flows in bfvi/bfvi/runs/cauchy_ablation/

#### Data for plotting
The script R/eval_cauchy_ablation_create_kl_plot_data.R creates the data files for the plotting in (R/runs/Cauchy1D)
#### Plotting
The figures (Fig 3a, 3b  at data 150822)  are created with  R/eval_cauchy_ablation.R
