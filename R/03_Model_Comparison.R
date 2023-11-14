####################################################################################
# Compare the results of the STAN models
####################################################################################
library(tidyverse)
library(dplR)
library(reshape2)
library(rstan)

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
#######################################################################################

# list all model files with .dbh.rds or .inc.rds ending
dbh.files <- paste0("outputs/", list.files(path = "outputs/",".dbh.rds"))
inc.files <- paste0("outputs/", list.files(path = "outputs/",".inc.rds"))

# read all of them in:
dbh.file.list <- lapply(dbh.files, readRDS)
inc.file.list <- lapply(inc.files, readRDS)

dbh.file.df <- do.call(rbind, dbh.file.list)
inc.file.df <- do.call(rbind, inc.file.list)

# make sure that each member of the list has the appropriate modelname
dbh.model.names.list <- lapply(list.files(path = "outputs/",".dbh.rds"), function(x){substr(x,0,  nchar(x) - 31)})
inc.model.names.list <- lapply(list.files(path = "outputs/",".inc.rds"), function(x){substr(x,0,  nchar(x) - 30)})


dbh.file.df$model <- rep(unlist(dbh.model.names.list), sapply(dbh.file.list, nrow)) # add a model column for each DBH fit
inc.file.df$model <- rep(unlist(inc.model.names.list), sapply(inc.file.list, nrow)) # add a model column for each inc fit

#######################################################################################
# Calculate model validation summary stats
#######################################################################################
# MSPE-mean squared predictive error
# RMSE-root mean squared predictive error

# PPL- calculating posterior predictive loss for model comparison:

#PPL = sum((Zobs - predZ)^2) - sum(var(predZ))
# for diameters
in.sample.validation.dbh.metrics <- dbh.file.df %>% group_by(model) %>% summarise(DBH_MSPE = mean((z.data-mean.ci)^2, na.rm =TRUE),
                                                                   DBH_RMSPE = sqrt(mean((z.data-mean.ci)^2, na.rm =TRUE)),
                                                                   DBH_MAPE = mean(abs(z.data-mean.ci), na.rm =TRUE) ,
                                                                   #V1 = mean(inc.data-mean.ci, na.rm =TRUE)/(sum(predvar)^(1/2))/n(), # estimate of bias in predictors over time (close to 0 = unbiased)
                                                                   #V2 = (mean((inc.data-mean.ci)^2, na.rm =TRUE)/(sum(predvar)/n()^(1/2))),  # estimate of accuracy of MSPEs (close to 1 = accurate)
                                                                   #V3 = (mean((inc.data-mean.ci)^2,na.rm =TRUE)^(1/2)),
                                                                   DBH_PPL = sum((z.data - mean.ci)^2, na.rm = TRUE) + sum(predvar, na.rm = TRUE)) # posterior predictive loss# goodness of fit estimate (small = better fit)
in.sample.validation.dbh.metrics$validation <- "in-sample"


# for increments
in.sample.validation.inc.metrics <- inc.file.df %>% group_by(model) %>% summarise(INC_MSPE = mean((inc.data-mean.ci)^2, na.rm =TRUE),
                                                                                  INC_RMSPE = sqrt(mean((inc.data-mean.ci)^2, na.rm =TRUE)),
                                                                                  INC_MAPE = mean(abs(inc.data-mean.ci), na.rm =TRUE), 
                                                                                  INC_PPL = sum((inc.data - mean.ci)^2, na.rm = TRUE) + sum(predvar, na.rm = TRUE)) # posterior predictive loss# goodness of fit estimate (small = better fit)
in.sample.validation.inc.metrics$validation <- "in-sample"

in.sample.validation <- left_join(in.sample.validation.dbh.metrics, in.sample.validation.inc.metrics)

#######################################################################################
# Generate some plots to compare model fits
#######################################################################################
ggplot(in.sample.validation, aes(x = DBH_MSPE, y = DBH_PPL, color = model))+geom_point()+theme_bw()
ggsave(height = 4, width = 6,"outputs/DBH_model_comparison.png")
ggplot(in.sample.validation, aes(x = INC_MSPE, y = INC_PPL, color = model))+geom_point()+theme_bw()
ggsave(height = 4, width = 6, "outputs/INC_model_comparison.png")


# Lower values of MSPE indicate a better model fit. Which model best predictions increment, and which is best at diameter?
# Based on this, which model should we choose?

# Note: These are all in-sample validation metrics. Additional DBH measurements could be used to do an
# "out-of-sample" model comparison. Out of sample increment data would also be helpful to do an out-of sample validation. 
