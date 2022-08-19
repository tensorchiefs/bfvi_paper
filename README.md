# bfvi_papern how to reproduce the results

# One dimensional examples.
The examples with only one parameter with a corresponding 1D posterior (Bernoulli/Cauchy) have been done using python. 

## Bernoulli:
### Workflow
#### Data generation 
The file `Python/Bernoulli/bernoulli_1D_F1F2.csv.gz` holds samples form the analytical solution, the Gaussian-VI, and the BF-Vi approximation. It is created by [Python/Bernoulli/Bernoulli_1D_F1F2.ipnb](https://github.com/tensorchiefs/bfvi_paper/blob/main/Python/Bernoulli/Bernoulli_1D_F1F2.ipynb) together with `Python/Bernoulli/vimlts_fast.py`

#### Plotting
The figure is created by `R/eval_Bernoulli_1D.R` using the data from `Python/Bernoulli/bernoulli_1D_F1F2.csv.gz`

## Cauchy
 
### Workflow
#### Data generation

##### MCMC 
To get a ground truth for the posterior we do MCMC
`mcmc/cauchy/mcmc_cauchy.stan`+ `R/eval_cauchy_ablation_create_kl_plots.R` -> `mcmc/cauchy/mcmc_samples_cauchy.rda` 

##### Gauss-VI: 
`Python/Cauchy/Cauchy.ipnb` → `Python/Cauchy/Gauss-VI_densities.csv.gz`
`Gauss-VI_densities.csv` contain samples `w` (first colom, in paper called xi) and `log_p` (second colom)

For Cauchy, we did a ablation study investigating the effects of different methods to transform between a predefined simple distribution and the variational posterior. Here we vary the simple transformation, and the transformations including the flexibility of the involved Bernstein polynomial.
##### BF-VI
Generated parameter files: 
Many ablation runs have been done with python (`bfvi/bfvi/cauchy_eval_ablation.py`) on GPU and CPU enviromemnts. The results are stored in from of the parameter files (columns are 'seed','epoch','az','bz', theta_1, ..., theat_M) holding information on the used seed for intialization and the number of epochs during training and the fitted parameters  of the transformation, i.e. the slope and intercept of the linear trafo and the parameters theta of the Bernstein polynomials, see e.g. `cauchy_eval_ablation_M_1_F1F2_params.csv` stored in  `bfvi/bfvi/runs/cauchy_ablation/`. For the creation of the paramter files manipulations (like commenting in/out the truncated normal) need to be done in addition the following files:

* `bfvi/bfvi/vimlt_fast.py`: define methods and layers
* `bfvi/bfvi/cauchy_ablation.py`: definitions of functions

#### Data for plotting
The script `R/eval_cauchy_ablation_create_kl_plot_data.R` creates the data files for the plotting from the parameter files in the repectives files in (`R/runs/Cauchy1D`)

#### Plotting
The figures (Fig 3a, 3b)  are created with  `R/eval_cauchy_ablation.R`

# Models with multiple parameters
The BF-VI results for the experiments with multi-dimensional posteriors have been created using `multidimensional_script.R` (the Melanoma experiments are an esception and are discussed below). In this script for the various models/data (like 8SCHOOLS_CP, DIAMONDS) the likelihood and the prior is defined. This script can be run with the follwing command line parameters:

* data: the model to be used like (`8SCHOOLS_CP` or `DIAMONDS`)
* method: the Bernstein Flow method that transforms from the simple predefined distribution to the variational posterior. We use `F1F2` in all experiments:  N(0,1) -> sigmoid -> Bernstein polynomial -> variational posterior distribution q. 
* M: the number of paramters in the Bernstein polynomial (we use 50 in all experiments) 
* T: The number of samples for the Monte-Carlo estimates during training (we use T=10 in all experiments).
* num_epochs: the number of epochs training (we use 100,000 in all experiments)

For example:
```
#The argument order is 
#data, method, num_epochs, M, T
R CMD BATCH --vanilla "--args run 8SCHOOLS_CP F1F2 100000 50 10" multidimensional_script.R
```
We run all experiment 5-times with random intializations. We sample from the fitted variational posterior q and store the samples `w` together with the variational log-posterior densities `log_qs`, the log-prior `L_prio` and log-likelihood `LLs` at those samples in a sample result file (e.g. `R/run/gpu_8SCHOOLS_F1F2_Epo_100000_M_50_T_10/samples_1.rda`) for all runs. In the same directory we also stroe the 5 loss histories (e.g. `R/run/gpu_8SCHOOLS_F1F2_Epo_100000_M_50_T_10/loss_hist_1.rda`). 

#### Ground thruth
The Stan-files for the ground thruth can be found in the [mcmc](https://github.com/tensorchiefs/bfvi_paper/tree/main/mcmc) directory.

#### Ploting
Using the corresponding eval-files (e.g `R/eval_8Schools_cp.R`) the metrics (e.g. `k-hat` with bootstraps CI) are calculated and plotted (e.g. `R/run/gpu_8SCHOOLS_F1F2_Epo_100000_M_50_T_10/k_hat.pdf`) and the final plots for the paper are procuded an stored in the `figures` directory.  

# Melanoma models M1, M2, and M3
Model M1 is based on images only, M2 on tabular data (age information) only, and M3 is semi-structured based on images and tabular data. The melanoma data have been downloaded from  https://challenge.isic-archive.com/data/#2020 and the images were downscaled to 128x128 pixels. The code (getData and resizeImages) can be found at `bfvi_paper/Ivonne_MA/functions/` or downloaded from https://www.dropbox.com/s/n1jodnzb71l3j8w/trainRes.zip?dl=0 

### Workflow
#### Raw Data 

##### M2 (only tabular data)
* MCMC posterior samples (`mcmc/mela_m2/mcmc_mela_m2.stan`, `R/eval_melanoma.r` --> `mcmc/mela_M2/mcmc_M2.csv.gz`
* BF-VI variational posterior samples (`multidimensional_script.R` store sample results in `R\runs\cpu_MELA_F1F2_Epo_100000_M_50_T_10`. Note that due to large data set size the fit was not possible on a standard GPU) 

###### M1 (only image data) and M3 (semistructed) 
The M1 and M3 models were fitted and M1's and M3's performance meassure were computed with `Ivonne_MA/Semi-structured_NN.ipynb`
With the same notebook, the samples from the fitted posterior for the age-slope were generated and stored in `Ivonne_MA/semi_posterior_slope_age.csv` (second half of the script)

#### Plotting
The melanoma figure in the paper is created by `R/eval_melanoma.r` which reads in posterior samples for the age-slope from MCMC, M2 BF-VI, M3 semi-structured and plots their estimated densities.



