# FIA_TR_fusion_tutorial

Code for a tutorial for fusing tree ring and forest inventory data in a Bayesian state-space model.

## Tutorial and Repository Structure
All analyses are run in the R/Rstudio environment, with STAN and rstan installed. For the tutorial during the workshop, we will be using POSIT Cloud to run through the tutorial, so you wont need to install anything on your own device, but if you want to work through this on your own, you will need a working installation of R and Rstudio (now Posit), as well as STAN and rstan. 

The repository has the folder `R`, which holds the R code for the tutorial. We will run through the following scripts to complete the tutorial, in the following order. 

1. `01_prep_TR_INV_DATA.R` -- The script compiles the tutorial data and formatting in a way to feed into the mode.
2. `02_Run_STAN_models.R` -- This script runs 4 state space models in stan (), and generates some diagnostic outputs.
3. `plot_pred_obs_data.R` -- We source this script for each model to generate and save tree-level predicted and observed plots.
4. `03_Model_Comparison.R` -- This script takes the outputs from steps 2 & 3 and does a simple model fit comparison
5. `04_Forecasting_SSM.R`-- This script is used to generate forecasts of tree growth from the posterior estimates of the "best-fit" model from 4.

All of the STAN model structures are stored in .stan files, in the `modelcode` folder. There are four models that we will call from `Run_STAN_models.R` and we will edit at least one of them. Briefly, we build up from a base state-space model with a tree-level random effect and a tree size fixed effect, then we include climate variables, and finally we include a competition effect. These models are simplified for the purpose of this tutorial. 

