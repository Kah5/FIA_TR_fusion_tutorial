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

All of the STAN model structures are stored in .stan files, in the `modelcode` folder. There are four models that we will call from `02_Run_STAN_models.R` and we will edit at least one of them. Briefly, we build up from a base state-space model with a tree-level random effect and a tree size fixed effect, then we include climate variables, and finally we include a competition effect. These models are simplified for the purpose of this tutorial. 

You can find few examples of random effect structures in the STAN code that we added during a workshop, which have not been run with the existing data, but are examples to work from when attempting to add in your own tree-ring and forest inventory data. These are in the 'modelcode' folder. 'model_3_plt_rand.stan' shows and example of adding a plot-level random effect with plot indexing, and 'model_3_plt_rand_MAP.stan' adds onto this with indexing for a climate variable where the data is provided at the plot-level rather than a data value for each tree. 

## Update after the Workshop

Based on discussions we had during the workshop in November 2023, we made some progress on improving speed and reducing redundancy in the state-space model. These updates are characterized in `model_3_nomissing_xscaled.stan`. 

Specifically we changed the following:
  1. removed the missing data model and rely on the generated quantities block for predicting 
  2. we are scaling each estimated tree size (x) by the mean and SD of the observed diameters (lines 59-60); if you are working with your own data, you will want to change this.
  3. the input tree ring and diameter data (y,x) are all read in as "long format," and indexed to the year and tree using indices in stan.

To see the datastructure and sample from the updated version of the model, look at `01_prep_TR_INV_DATA_longformat.R` and `02_Run_STAN_model_longformat.R`

