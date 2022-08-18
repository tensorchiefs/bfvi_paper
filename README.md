# bfvi_paper

# How to Reproduce the results

## Bernoulli:
### Workflow
#### Raw Data (`Python/Bernoulli/bernoulli_1D_F1F2.csv.gz`)
This file includes samples form the analytical solution, the Gaussian-VI, and the BF-Vi approximation. It is created by [Python/Bernoulli/Bernoulli_1D_F1F2.ipnb](https://github.com/tensorchiefs/bfvi_paper/blob/main/Python/Bernoulli/Bernoulli_1D_F1F2.ipynb) together with `Python/Bernoulli/vimlts_fast.py`

#### Plotting
The figure is created by `R/eval_Bernoulli_1D.R` using the data from `Python/Bernoulli/bernoulli_1D_F1F2.csv.gz`

## Cauchy

### Workflow
#### Raw Data
##### Gauss-VI: `Python/Cauchy/Cauchy.ipnb` → `Python/Cauchy/Gauss-VI_densities.csv.gz`
`Gauss-VI_densities.csv` contain samples `w` (first colom, in paper called xi) and `log_p` (second colom)

##### MCMC 
[mcmc/cauchy/mcmc_cauchy.stan](https://github.com/tensorchiefs/bfvi_paper/blob/main/mcmc/cauchy/mcmc_cauchy.stan) + R/eval_cauchy_ablation_create_kl_plots.R -> mcmc/cauchy/mcmc_samples_cauchy.rda 

##### BF-VI
Parameter-Files 
Many ablation runs have been done with python (`bfvi/bfvi/cauchy_eval_ablation.py`) on GPU and CPU enviromemnts. The results are stored in from of the parameters defining the Bernstein flows in `bfvi/bfvi/runs/cauchy_ablation/`. For the creation of the paramter files manipulations (like commenting in/out the truncated normal) in:

* bfvi/bfvi/vimlt_fast.py: define methods and layers
* bfvi/bfvi/cauchy_ablation.py: definitions of functions

#### Data for plotting
The script R/eval_cauchy_ablation_create_kl_plot_data.R creates the data files for the plotting in (R/runs/Cauchy1D)
#### Plotting
The figures (Fig 3a, 3b at data 150822)  are created with  R/eval_cauchy_ablation.R

# Models with multiple parameters
The results for the remaining experiments (except Melanoma, which is described below) have created as follows. First using the script `multidimensional_script` with the follwing command line parameters T=10 M=50 method=’F1F2’ Epochs = 100000 reps=5 .

for example
```
#The argument order is 
#data, method, num_epochs, M, T
R CMD BATCH --vanilla "--args run 8SCHOOLS_CP F1F2 100000 50 10" multidimensional_script.R
```
produces samples from the posterior (`w` together with the variational posterior densities `log_qs`, the prior `L_prio` and likelihood `LLs` at those samples). Fore.g. `R/run/gpu_8SCHOOLS_F1F2_Epo_100000_M_50_T_10/samples_1.rda`) for 5 runs. Also storted in the directories are the loss histories. 

Using the corresponding eval-files (e.g `R/eval_8Schools_cp.R`) the metrics are calculated and the final plots are procuded. The metrics comprise `k-hat` with bootstraps CI. The Stan-files for the ground thruth can be found in the `mcmc` directory.


## Melanoma
Pre

### Workflow
#### Raw Data 


This includes samples for the age-slope posterior: MCMC samples (mcmc/mela_m2/mcmc_mela_m2.stan, R/eval_melanoma.r --> mcmc/mela_M2/mcmc_M2.csv.gz, M2 BF-VI samples (multidimensional_script.R store samples in R\runs\cpu_MELA_F1F2_Epo_100000_M_50_T_10, note that due to large data set size the fit was not possible on a standard GPU) , samples from semistructured M3 model (From Ivonne's master thesis, `Ivonne_MA/Semi-structured_NN.ipynb` (https://github.com/IvonneKo/ISIC_TMVI), samples stored in Ivonne_MA/semi_posterior_slope_age.csv) . 

#### Plotting
The figure is created by `R/eval_melanoma.r` which reads in posterior samples for the age-slope from MCMC, M2 BF-VI, M3 semi-structured.



