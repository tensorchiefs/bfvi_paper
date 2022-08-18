# bfvi_papern / How to Reproduce the results

# One dimensional examples.
For the one dimensional examples (Bernoulli/Cauchy) have been done using python. 

## Bernoulli:
### Workflow
#### Raw Data (`Python/Bernoulli/bernoulli_1D_F1F2.csv.gz`)
This file includes samples form the analytical solution, the Gaussian-VI, and the BF-Vi approximation. It is created by [Python/Bernoulli/Bernoulli_1D_F1F2.ipnb](https://github.com/tensorchiefs/bfvi_paper/blob/main/Python/Bernoulli/Bernoulli_1D_F1F2.ipynb) together with `Python/Bernoulli/vimlts_fast.py`

#### Plotting
The figure is created by `R/eval_Bernoulli_1D.R` using the data from `Python/Bernoulli/bernoulli_1D_F1F2.csv.gz`

## Cauchy
For Cauchy, we did a ablation study investigating the effects of different parts of the transformation. 
### Workflow
#### Raw Data
##### MCMC 
[mcmc/cauchy/mcmc_cauchy.stan](https://github.com/tensorchiefs/bfvi_paper/blob/main/mcmc/cauchy/mcmc_cauchy.stan) + R/eval_cauchy_ablation_create_kl_plots.R -> mcmc/cauchy/mcmc_samples_cauchy.rda 

##### Gauss-VI: `Python/Cauchy/Cauchy.ipnb` → `Python/Cauchy/Gauss-VI_densities.csv.gz`
`Gauss-VI_densities.csv` contain samples `w` (first colom, in paper called xi) and `log_p` (second colom)

##### BF-VI
Parameter-Files 
Many ablation runs have been done with python (`bfvi/bfvi/cauchy_eval_ablation.py`) on GPU and CPU enviromemnts. The results are stored in from of the parameters defining the Bernstein flows in respective directories of `bfvi/bfvi/runs/cauchy_ablation/`. For the creation of the paramter files manipulations (like commenting in/out the truncated normal) need to be done in addtion the following files:

* `bfvi/bfvi/vimlt_fast.py`: define methods and layers
* `bfvi/bfvi/cauchy_ablation.py`: definitions of functions

#### Data for plotting
The script `R/eval_cauchy_ablation_create_kl_plot_data.R` creates the data files for the plotting from the parameter files in the repectives files in (`R/runs/Cauchy1D`)

#### Plotting
The figures (Fig 3a, 3b)  are created with  `R/eval_cauchy_ablation.R`

# Models with multiple parameters
The BF-VI results for the remaining experiments (except Melanoma M2, which is described below) have created using `multidimensional_script.R`. In this script for the various models/data (like 8Schools_CP) the likelihood and the prior is defined. This script can be run with the follwing command line parameters T=10 M=50 method=’F1F2’ Epochs = 100000 reps=5 .

for example:
```
#The argument order is 
#data, method, num_epochs, M, T
R CMD BATCH --vanilla "--args run 8SCHOOLS_CP F1F2 100000 50 10" multidimensional_script.R
```
produces samples from the posterior (`w` together with the variational posterior densities `log_qs`, the prior `L_prio` and likelihood `LLs` at those samples). For e.g. (`R/run/gpu_8SCHOOLS_F1F2_Epo_100000_M_50_T_10/samples_1.rda`) for 5 runs. Also storted in the directories are the loss histories. 

#### Ground thruth
The Stan-files for the ground thruth can be found in the [mcmc](https://github.com/tensorchiefs/bfvi_paper/tree/main/mcmc) directory.

#### Ploting
Using the corresponding eval-files (e.g `R/eval_8Schools_cp.R`) the metrics are calculated and the final plots are procuded. The metrics comprise `k-hat` with bootstraps CI. 

# Semistructured Models (Melanoma)
Preprocessing M1, and M3 work on images. These images have been downscaled to 128x128 pixels. The code (getData and resizeImages) can be found at `bfvi_paper/Ivonne_MA/functions/` or downloaded from https://www.dropbox.com/s/n1jodnzb71l3j8w/trainRes.zip?dl=0 

### Workflow
#### Raw Data 

##### M2 (only tabular data)
* MCMC samples (`mcmc/mela_m2/mcmc_mela_m2.stan`, `R/eval_melanoma.r` --> `mcmc/mela_M2/mcmc_M2.csv.gz`
* BF-VI samples (`multidimensional_script.R` store samples in `R\runs\cpu_MELA_F1F2_Epo_100000_M_50_T_10`. Note that due to large data set size the fit was not possible on a standard GPU) 

###### M1 (only image data) and M3 (semistructed) 
Performance meassure for M1 can be produced with `Ivonne_MA/Semi-structured_NN.ipynb`
Performance meassure and samples for semistructured M3 model `Ivonne_MA/Semi-structured_NN.ipynb` (second half of the script)

#### Plotting
The figure is created by `R/eval_melanoma.r` which reads in posterior samples for the age-slope from MCMC, M2 BF-VI, M3 semi-structured.



