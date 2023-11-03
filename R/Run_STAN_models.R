#-------------------------------------------------------------------------------
# RUN STAN MODELS
#-------------------------------------------------------------------------------
# We will run three versions of the state-space model, two of which are already coded up for you: 
# `model_noyearRE.stan` and `model_noyearRE_X.stan`

options(mc.cores = parallel::detectCores())
# --------------------------------
# Run model 1: TREE size only
# --------------------------------
# first run a model that only has a tree-level random effect and a tree size effect
model.name <- "tree.size.only.model"

# there are a few different ways to actually "run" the stan model, here we demonstrate compiliing the model:
model1 <- stan_model(file = "modelcode/model_treesize_only.stan")

# next we do the sampling steps of the model
model.1 <- sampling(model1, # path to the .STAN file, which contains the actual model code
                    data = model.data, # the data to use in the model
                    iter = 100, # the number of iterations to run the model for--Here we will only run for 100 iterations
                    # but for the model parameters to be trusted we would want to run the model for longer ~>3000 iterations
                    # for the purpose of timing, we will have you just read in the fully converged model outputs 
                    chains = 2, # the number of chains to run 
                    verbose=FALSE, 
                    control =  list(max_treedepth = 15), # control is where you can dig into specifics of the sampler
                    # more info here: 
                    #sample_file = model.name, #the sample_file argument will save samples in .csv files as you go
                    #adapt_delta = 0.99,
                    
                    # important: STAN won't include estimates of a parameter in the output if you 
                    # don't tell it to in the pars argument
                    pars = c("mu", "sigma_inc", "sigma_add", "sigma_dbh","betaX", "alpha_TREE","x", "inc"))


# we could also just compile and run the model using the stan() call:
# model.1 <- stan(file = "modelcode/model_treesize_only.stan", # path to the .STAN file, which contains the actual model code
#                 data = model.data, # the data to use in the model
#                 iter = 100, # the number of iterations to run the model for
#                 chains = 2, # the number of chains to run 
#                 verbose=FALSE, 
#                 control =  list(max_treedepth = 15), # control is where you can dig into specifics of the sampler
#                 # more info here: 
#                 #sample_file = model.name, #the sample_file argument will save samples in .csv files as you go
#                 #adapt_delta = 0.99,
#                 
#                 # important: STAN won't include estimates of a parameter in the output if you 
#                 # don't tell it to in the pars argument
#                 pars = c("mu", "sigma_inc", "sigma_add", "sigma_dbh","betaX", "alpha_TREE","x", "inc"))
# 
# #save the outputs so we don't need to rerun later if r crashes or we need to close out
saveRDS(model.1, "outputs/model.1.short.RDS")

# Check out the model diagnostics for this short sample:

# get model diagnostics and save these to look at
sampler_params <- get_sampler_params(model.1, inc_warmup = FALSE)

mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
sum_divergent_transitions_by_chain <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
sampler_diag <- data.frame(model.name = rep(model.name, 2),
                           chain = 1:2, 
                           accept = mean_accept_stat_by_chain, 
                           sum_divergent_transitions_by_chain = sum_divergent_transitions_by_chain, 
                           n.samples = nrow(sampler_params[[1]]))
# we want few (no) divergent transitions so this is good
sampler_diag
# there are no divergent transistions in either chain, and acceptance rate > 0.89
write.csv(sampler_diag,paste0("outputs/", model.name, "not_converged_sample_diag.csv"), row.names = FALSE)


# get the convergence statistics of the model:
fit_ssm_df <- as.data.frame(model.1) # takes awhile to convert to df
Rhats <- apply(fit_ssm_df, 2, Rhat)
hist(Rhats)
# most of the R hat values are below 1.01
ESS_bulks <- apply(fit_ssm_df, 2, ess_bulk)
hist(ESS_bulks)
ESS_tails <- apply(fit_ssm_df, 2, ess_tail)
hist(ESS_tails)

convergence.stats <- as.data.frame(rbind(Rhats, ESS_bulks, ESS_tails))
convergence.stats$Statistic <- c("Rhat", "ESS_bulk", "ESS_tail")

write.csv(convergence.stats, paste0("outputs/",model.name, "not_converged_convergence_stats.csv"))

# also generate some plots of the parameters

posterior <- as.array(model.1)
par.names = c("mu", "sigma_inc", "sigma_add", "sigma_dbh", "betaX","alpha_TREE[1]")

# traceplots show the values of each parameter over the number of samples. Values should vary, but center around a single value
#png(height = 4, width = 7, units = "in", res = 100, paste0("outputs/not_converged_traceplots_tutorial_", model.name, ".png"))
traceplot (model.1, pars = par.names, nrow = 8, ncol = 4, inc_warmup = FALSE) 
#dev.off()

# Now that we have seen what traceplots look like when they have not converged, lets read in the output from the same model, but with 3000 samples
model.1 <- readRDS("outputs/model.1.RDS")

png(height = 4, width = 7, units = "in", res = 100, paste0("outputs/traceplots_tutorial_", model.name, ".png"))
traceplot (model.1, pars = par.names, nrow = 8, ncol = 4, inc_warmup = FALSE) 
dev.off()
# check the pairs plots for posterior correlations
# posterior correlations in the model might be expected, for example between an intercept and a slope, 
# or they could cause problems with sampling/identifiablility
png(height = 7, width = 7, units = "in", res = 100, paste0("tutorial/outputs/pairs_plot_tutorial_", model.name, ".png"))
pairs(model.1, pars = c("mu", "sigma_inc", "sigma_add", "sigma_dbh", "betaX","alpha_TREE[1]", "alpha_TREE[2]"))
dev.off()
# sigma_inc and sigma_add have somme posterior correlations, as do the alpha_TREE values with the tree size effect

# create parameter plots to visualize the values of the effects:
covariates = c( "betaX")

cov.estimates <- fit_ssm_df %>% dplyr::select(covariates)
sigma <- fit_ssm_df[,"sigma_add"] # get additive process error

# Tree-level random effects 
alpha_trees <- dplyr::select(fit_ssm_df, "alpha_TREE[1]":paste0("alpha_TREE[", model.data$Nrow, "]"))


# plot year and tree random effects:
alpha_tree.m <- reshape2::melt(alpha_trees)
tree.quant <- alpha_tree.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                                ci.lo = quantile(value, 0.025, na.rm =TRUE),
                                                                ci.hi = quantile(value, 0.975, na.rm =TRUE))


ggplot()+geom_point(data = tree.quant, aes(x = variable, y = median))+
  geom_errorbar(data =tree.quant, aes(x = variable, ymin = ci.lo, ymax = ci.hi), linewidth = 0.1)+theme_bw()+
  theme(axis.text = element_text(angle= 45, hjust = 1), panel.grid = element_blank())+
  ylab("Estimated effect")

ggsave(here("tutorial/outputs", paste0("tree_random_effects_", model.name, ".png")))


# Make a violin plot of the fixed effects. 
covariates.m <- reshape2::melt(cov.estimates)
head(covariates.m)

# make violin plots of the samples to look at the medians & the uncertainty around each covariate:
ggplot(covariates.m, aes(x = variable, y = value))+geom_violin(draw_quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), fill = "grey")+theme_bw()
ggsave(here("tutorial/outputs", paste0("covariate_violin_plots_", model.name, ".png")))

# make plots of predicted vs observed for each tree/year
# comparing predicted vs. observed
# lets compare the predicted vs observed for the increments and the diameters:

# seelect the diameter estimates (x.pred) from the model.3 dataframe
x.pred <- dplyr::select(as.data.frame(model.1),"x[1,1]":paste0("x[", model.data$Nrow, ",", model.data$Ncol,"]"))
# select the increment estimates (inc) from the model.3 dataframe
inc.pred <- dplyr::select(as.data.frame(model.1),"inc[1,1]":paste0("inc[", model.data$Nrow, ",", model.data$Ncol,"]"))

model.out <- cbind(x.pred, inc.pred) # get this to make plots
mod.data <- model.data
source("tutorial/plot_pred_obs_data.R") # run script to make predicted vs obs plots
# check the outputs folder, you should have two new pdfs, two new png files, and two new RDS files
# this code plotted the predicted vs observed for this model 


# --------------------------------------------------------------------------------------------------------
# Run model 2: TREE size + two time-varying climate variables (Maximum temperature + total precipitation)
# --------------------------------------------------------------------------------------------------------
# Second, run a model that has tree size and a single climate variable in in it
# Normally, you would only add in one variable at a time so if something goes wrong, you know which addition to the model is causing an issue
# this model takes a little longer to run, so we will just run it with 2 chains to get some results
model.name <- "TmaxPPTX"
# add in the tree size effect with out the year random effects but with climate
model.2 <- stan(file = "tutorial/modelcode/model_noyearRE_X.stan",
                data = model.data,
                iter = 3000, chains = 2, verbose=FALSE, control =  list(max_treedepth = 15),
                #sample_file = model.name,
                #adapt_delta = 0.99,
                pars = c("mu", "sigma_inc", "sigma_add", "sigma_dbh","betaTmax", "betaPrecip", "betaX","alpha_TREE","x", "inc"))

saveRDS(model.2, "tutorial/outputs/model.2.RDS")
# get model diagnostics and save these to look at
sampler_params <- get_sampler_params(model.2, inc_warmup = FALSE)

mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
sum_divergent_transitions_by_chain <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
sampler_diag <- data.frame(model.name = rep(model.name, 2),
                           chain = 1:2, 
                           accept = mean_accept_stat_by_chain, 
                           sum_divergent_transitions_by_chain = sum_divergent_transitions_by_chain, 
                           n.samples = nrow(sampler_params[[1]]))
# we want few (no) divergent transitions so this is good
sampler_diag
# there are no divergent transistions in either chain, and acceptance rate > 0.89
write.csv(sampler_diag, here("tutorial/outputs", paste0(model.name, "_sample_diag.csv")), row.names = FALSE)


# get the convergence statistics of the model:
fit_ssm_df <- as.data.frame(model.2) # takes awhile to convert to df
Rhats <- apply(fit_ssm_df, 2, Rhat)
hist(Rhats)
# most of the R hat values are below 1.01
ESS_bulks <- apply(fit_ssm_df, 2, ess_bulk)
hist(ESS_bulks)
ESS_tails <- apply(fit_ssm_df, 2, ess_tail)
hist(ESS_tails)

convergence.stats <- as.data.frame(rbind(Rhats, ESS_bulks, ESS_tails))
convergence.stats$Statistic <- c("Rhat", "ESS_bulk", "ESS_tail")

write.csv(convergence.stats, here("tutorial/outputs", paste0(model.name, "_convergence_stats.csv")))

# also generate some plots of the parameters

posterior <- as.array(model.2)
par.names = c("mu", "sigma_inc", "sigma_add", "sigma_dbh", "betaX","betaTmax", "betaPrecip", "alpha_TREE[1]")

# traceplots show the values of each parameter over the number of samples. Values should vary, but center around a single value
png(height = 4, width = 7, units = "in", res = 100, paste0("tutorial/outputs/traceplots_tutorial_", model.name, ".png"))
traceplot (model.2, pars = par.names, nrow = 8, ncol = 4, inc_warmup = FALSE) 
dev.off()

# check the pairs plots for posterior correlations
# posterior correlations in the model might be expected, for example between an intercept and a slope, 
# or they could cause problems with sampling/identifiablility
png(height = 7, width = 7, units = "in", res = 100, paste0("tutorial/outputs/pairs_plot_tutorial_", model.name, ".png"))
pairs(model.2, pars = c("mu", "sigma_inc", "sigma_add", "sigma_dbh", "betaX","betaPrecip", "alpha_TREE[1]", "alpha_TREE[2]"))
dev.off()
# sigma_inc and sigma_add have somme posterior correlations, as do the alpha_TREE values with the tree size effect

# create parameter plots to visualize the values of the effects:
covariates = c( "betaX", "betaTmax", "betaPrecip")

cov.estimates <- fit_ssm_df %>% dplyr::select(covariates)
sigma <- fit_ssm_df[,"sigma_add"] # get additive process error

# Tree-level random effects 
alpha_trees <- dplyr::select(fit_ssm_df, "alpha_TREE[1]":paste0("alpha_TREE[", model.data$Nrow, "]"))


# plot year and tree random effects:
alpha_tree.m <- reshape2::melt(alpha_trees)
tree.quant <- alpha_tree.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                                ci.lo = quantile(value, 0.025, na.rm =TRUE),
                                                                ci.hi = quantile(value, 0.975, na.rm =TRUE))


ggplot()+geom_point(data = tree.quant, aes(x = variable, y = median))+
  geom_errorbar(data =tree.quant, aes(x = variable, ymin = ci.lo, ymax = ci.hi), linewidth = 0.1)+theme_bw()+
  theme(axis.text = element_text(angle= 45, hjust = 1), panel.grid = element_blank())+
  ylab("Estimated effect")

ggsave(here("tutorial/outputs", paste0("tree_random_effects_", model.name, ".png")))


# Make a violin plot of the fixed effects. 
covariates.m <- reshape2::melt(cov.estimates)
head(covariates.m)

# make violin plots of the samples to look at the medians & the uncertainty around each covariate:
ggplot(covariates.m, aes(x = variable, y = value))+geom_violin(draw_quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), fill = "grey")+theme_bw()
ggsave(here("tutorial/outputs", paste0("covariate_violin_plots_", model.name, ".png")))

# make plots of predicted vs observed for each tree/year
# comparing predicted vs. observed
# lets compare the predicted vs observed for the increments and the diameters:

# seelect the diameter estimates (x.pred) from the model.2 dataframe
x.pred <- dplyr::select(as.data.frame(model.2),"x[1,1]":paste0("x[", model.data$Nrow, ",", model.data$Ncol,"]"))
# select the increment estimates (inc) from the model.2 dataframe
inc.pred <- dplyr::select(as.data.frame(model.2),"inc[1,1]":paste0("inc[", model.data$Nrow, ",", model.data$Ncol,"]"))

model.out <- cbind(x.pred, inc.pred) # get this to make plots
mod.data <- model.data
source("tutorial/plot_pred_obs_data.R") # run script to make predicted vs obs plots
# check the outputs folder, you should have two new pdfs, two new png files, and two new RDS files
# this code plotted the predicted vs observed for this model 


#--------------------------------------------------------------------------------------
# Run a third model, but this time add in a plot level competition variable SDI
#--------------------------------------------------------------------------------------
# Make a copy of the stan code for the previous model, and save as "model_treesize_climate_sdi.stan"
# save it in the modelcode folder

# in the data block, tell stan to expect SDI by pasting the following, without the #:
# real SDI[Nrow]; //Vector of SDI values for each plot 

# in the parameters block, define the parameter betaSDI by pasting the following, without the #:
# real betaSDI;

# in the model block, set the prior by pasting the following, without the # just under the priors for betaPrecip:
# betaSDI ~ normal(0, 10);

# also in the model block, we need to add betaSDI to the equation for inc[i,t], so add
# + betaSDI*SDI[i] to the end of the equation, so that the whole line defining inc[i,t] reads:
# inc[i,t] ~ lognormal(alpha_TREE[i] + betaTmax*tmaxAprMayJunscaled[i,t]+ betaPrecip*wateryrscaled[i,t] + betaSDI*SDI[i], sigma_add);

# before running the model you can hit the "check" button to check that the stan code syntax is correct

# run the model here
model.name <- "TmaxPPTX_SDI"
# add in the tree size effect with out the year random effects but with climate
model.3 <- stan(file = "tutorial/modelcode/model_treesize_climate_sdi.stan",
                data = model.data,
                iter = 3000, chains = 2, verbose=FALSE, control =  list(max_treedepth = 15),
                #sample_file = model.name,
                #adapt_delta = 0.99,
                pars = c("mu", "sigma_inc", "sigma_add", "sigma_dbh","betaTmax", "betaPrecip", "betaX",
                         "betaSDI","alpha_TREE","x", "inc"))


# get model diagnostics and save these to look at
sampler_params <- get_sampler_params(model.3, inc_warmup = FALSE)

mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
sum_divergent_transitions_by_chain <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
sampler_diag <- data.frame(model.name = rep(model.name, 2),
                           chain = 1:2, 
                           accept = mean_accept_stat_by_chain, 
                           sum_divergent_transitions_by_chain = sum_divergent_transitions_by_chain, 
                           n.samples = nrow(sampler_params[[1]]))
# we want few (no) divergent transitions so this is good
sampler_diag
# there are no divergent transistions in either chain, and acceptance rate > 0.89
write.csv(sampler_diag, here("tutorial/outputs", paste0(model.name, "_sample_diag.csv")), row.names = FALSE)


# get the convergence statistics of the model:
fit_ssm_df <- as.data.frame(model.3) # takes awhile to convert to df
Rhats <- apply(fit_ssm_df, 2, Rhat)
hist(Rhats)
# most of the R hat values are below 1.01
ESS_bulks <- apply(fit_ssm_df, 2, ess_bulk)
hist(ESS_bulks)
ESS_tails <- apply(fit_ssm_df, 2, ess_tail)
hist(ESS_tails)

convergence.stats <- as.data.frame(rbind(Rhats, ESS_bulks, ESS_tails))
convergence.stats$Statistic <- c("Rhat", "ESS_bulk", "ESS_tail")

write.csv(convergence.stats, here("tutorial/outputs", paste0(model.name, "_convergence_stats.csv")))

# also generate some plots of the parameters

posterior <- as.array(model.3)
par.names = c("mu", "sigma_inc", "sigma_add", "sigma_dbh", "betaX","betaTmax", "betaPrecip", "alpha_TREE[1]")

# traceplots show the values of each parameter over the number of samples. Values should vary, but center around a single value
png(height = 4, width = 7, units = "in", res = 100, paste0("tutorial/outputs/traceplots_tutorial_", model.name, ".png"))
traceplot (model.3, pars = par.names, nrow = 8, ncol = 4, inc_warmup = FALSE) 
dev.off()

# check the pairs plots for posterior correlations
# posterior correlations in the model might be expected, for example between an intercept and a slope, 
# or they could cause problems with sampling/identifiablility
png(height = 7, width = 7, units = "in", res = 100, paste0("tutorial/outputs/pairs_plot_tutorial_", model.name, ".png"))
pairs(model.3, pars = c("mu", "sigma_inc", "sigma_add", "sigma_dbh", "betaX","betaPrecip", "alpha_TREE[1]", "alpha_TREE[2]"))
dev.off()
# sigma_inc and sigma_add have somme posterior correlations, as do the alpha_TREE values with the tree size effect

# create parameter plots to visualize the values of the effects:
covariates = c( "betaX", "betaTmax", "betaPrecip", "betaSDI")

cov.estimates <- fit_ssm_df %>% dplyr::select(covariates)
sigma <- fit_ssm_df[,c("sigma_add", "sigma_inc", "sigma_dbh")] # get additive process error

# Tree-level random effects 
alpha_trees <- dplyr::select(fit_ssm_df, "alpha_TREE[1]":paste0("alpha_TREE[", model.data$Nrow, "]"))


# plot year and tree random effects:
alpha_tree.m <- reshape2::melt(alpha_trees)
tree.quant <- alpha_tree.m %>% group_by(variable) %>% summarise(median = quantile(value, 0.5, na.rm =TRUE),
                                                                ci.lo = quantile(value, 0.025, na.rm =TRUE),
                                                                ci.hi = quantile(value, 0.975, na.rm =TRUE))


ggplot()+geom_point(data = tree.quant, aes(x = variable, y = median))+
  geom_errorbar(data =tree.quant, aes(x = variable, ymin = ci.lo, ymax = ci.hi), linewidth = 0.1)+theme_bw()+
  theme(axis.text = element_text(angle= 45, hjust = 1), panel.grid = element_blank())+
  ylab("Estimated effect")

ggsave(here("tutorial/outputs", paste0("tree_random_effects_", model.name, ".png")))


# Make a violin plot of the fixed effects. 
covariates.m <- reshape2::melt(cov.estimates)
head(covariates.m)

# make violin plots of the samples to look at the medians & the uncertainty around each covariate:
ggplot(covariates.m, aes(x = variable, y = value))+geom_violin(draw_quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), fill = "grey")+theme_bw()
ggsave(here("tutorial/outputs", paste0("covariate_violin_plots_", model.name, ".png")))

# make plots of predicted vs observed for each tree/year
# comparing predicted vs. observed
# lets compare the predicted vs observed for the increments and the diameters:

# seelect the diameter estimates (x.pred) from the model.3 dataframe
x.pred <- dplyr::select(as.data.frame(model.3),"x[1,1]":paste0("x[", model.data$Nrow, ",", model.data$Ncol,"]"))
# select the increment estimates (inc) from the model.3 dataframe
inc.pred <- dplyr::select(as.data.frame(model.3),"inc[1,1]":paste0("inc[", model.data$Nrow, ",", model.data$Ncol,"]"))

model.out <- cbind(x.pred, inc.pred) # get this to make plots
mod.data <- model.data
source("tutorial/plot_pred_obs_data.R") # run script to make predicted vs obs plots
# check the outputs folder, you should have two new pdfs, two new png files, and two new RDS files
# this code plotted the predicted vs observed for this model 

# Next, open the Model_Comparison.R Script to do a model evalution & comparison of the outputs of these models.
