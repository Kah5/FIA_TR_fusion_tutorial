####################################################################################
# Compare the results of the STAN models
####################################################################################

# After each model, we ran the plot_pred_obs_data.R code, which generated some outputs in the outputs folder
# these files inlcude dataframes that have the actual values for diameter for each tree in each year, compared to
# the mean and 95% Credible intervals of diameter esitimates. There is also a file 
# containing predicted and observed values for increments. 

# Here, we will read these in and caluclate some metrics to evaulate each model based on 
# in-sample predictions. With a full model validation, we would also compare to some 
# held-out values (for either increments and/or diameters). For the purpose of this
# tutorial we will just do a comparison on the in-sample values, that is the values
# used to fit the model.

#######################################################################################
# Read in the predicted vs observed values of increment

# model 1
m1.dbh <- readRDS("tree.size.only.modelpred.obs.out.of.sample.dbh.rds")
m1.inc <- readRDS("tutorial/outputs/tree.size.only.modelpred.obs.out.of.sample.inc.rds")


# model 2


# model 3
m3.dbh <- readRDS("tutorial/outputs/TmaxPPTX.pred.obs.out.of.sample.dbh.rds")
m3.inc <- readRDS("tutorial/outputs/TmaxPPTXpred.obs.out.of.sample.inc.rds")

# model 4

#######################################################################################
# Combine all data frames & calculate model validation summary stats

#######################################################################################
# Generate some plots to compare model fits

