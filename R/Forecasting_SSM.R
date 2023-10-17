################################################################################
# Code to generate some forecasts from each model
################################################################################

# There are a few different ways to generate forecasts using this modeling approach:
# 1. Within the stan itself, either in the generated quantities block, or as missing data.
# 2. Using posterior estimates of the state space model parameters to make projections in R.

# We will walk through how to generate forecasts from the posterior estimates of the state space model
################################################################################

################################################################################
# Step 1: Compile covariate data to generate forecasts with 
################################################################################
# Here, we will assume that SDI is constant (though it should change with changing diameters)
# We will use downscaled future climate values from CMIP 5 to generate forecasts of
# tree growth under future climate emissions scenarios (and include climate projection uncertainty).
# We downloaded climate model runs for 21 different climate models, for 4 rcp emissions scenarios from:
# 