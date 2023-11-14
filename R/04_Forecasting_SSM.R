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
# We downloaded climate model runs for 21 different climate models, for 4 rcp emissions scenarios from:rea
#-------------------------------------------------------------------
# set up to forecast through time from just the last X condition
# using future climate time series (2018-2099)
# note that this version also includes driver uncertainty
#-------------------------------------------------------------------


# lets read in the outputs from model.3
# now load in the longer output to do the model checking:
model.3 <- readRDS ("longOutputs/model.3.RDS")
model.name <- "TmaxPPTX_SDI"
fit_ssm_df <- as.data.frame(model.3) # takes awhile to convert to df

# get the covariates from the posterior estimates
covariates = c( "betaX", "betaTmax", "betaPrecip", "betaSDI")

cov.estimates <- fit_ssm_df %>% dplyr::select(covariates)
sigma <- fit_ssm_df[,c("sigma_add", "sigma_inc", "sigma_dbh")] # get additive process error

# Tree-level random effects 
alpha_trees <- dplyr::select(fit_ssm_df, "alpha_TREE[1]":paste0("alpha_TREE[", model.data$Nrow, "]"))


# plot year and tree random effects:
alpha_tree.m <- reshape2::melt(alpha_trees)

# covariates.m
covariates.m <- reshape2::melt(cov.estimates)
head(covariates.m)

# get the posterior estimates of tree diameters 
x.pred <- dplyr::select(as.data.frame(model.3),"x[1,1]":paste0("x[", model.data$Nrow, ",", model.data$Ncol,"]"))
# select the increment estimates (inc) from the model.3 dataframe
inc.pred <- dplyr::select(as.data.frame(model.3),"inc[1,1]":paste0("inc[", model.data$Nrow, ",", model.data$Ncol,"]"))

model.out <- cbind(x.pred, inc.pred) # get this to make plots



# driver uncertainty:
# for each scenario, in each year, for each plot, we sample from the distribution of ppt and tmax
# that has a mean of the ensemble mean, and variance of the ensemble variance

# read the rds with future climate scenarios:
clim.ts <- readRDS("data/future_climate_list.rds")

# get the historical climate to scale the future climate by:
clim.ts.df <- clim.ts$future.climate.ts
clim.ts.df$tmax.fall.spr[is.nan(clim.ts.df$tmax.fall.spr)] <- NA
clim.ts.df$tmax.fall.spr[is.nan(clim.ts.df$tmax.fall.spr)] <- NA
tmax.fallspr.df <- clim.ts$time_data$tmax.fallspr[!duplicated(clim.ts$time_data$tmax.fallspr),]
ppt.full.df <- clim.ts$time_data$wintP.wateryr[!duplicated(clim.ts$time_data$wintP.wateryr),]


# need to scale future climate data on the same scale as the past climate (but need to do it by plot?)
clim.ts.df$ppt.scale <-(clim.ts.df$year.ppt-mean(as.matrix(ppt.full.df)))/sd(as.matrix(ppt.full.df))
clim.ts.df$tmax.scaled <-(clim.ts.df$tmax.fall.spr-mean(as.matrix(tmax.fallspr.df)))/sd(as.matrix(tmax.fallspr.df))


# filter out for just the trees selected in this tutorial
head(clim.ts.df)

# read in the random tree id indices so we can match them up with the climate data:
tree.id <- saveRDS("data/tree.id.rand.sample.rds")

tree.id.numeric <- substring(as.character(tree.id), 2, nchar(as.character(tree.id))) # use the tree id from the first script to select ID in future.clim.subset to scale

# set up the 4 different emissions scenarios data, but we will just use RCP2.6
future.clim.subset.26 <- clim.ts.df %>% filter(rcp %in% "rcp26" & ID %in% unique(tree.id.numeric))
future.clim.subset.45 <- clim.ts.df %>% filter(rcp %in% "rcp45" & ID %in% unique(tree.id.numeric))
future.clim.subset.60 <- clim.ts.df %>% filter(rcp %in% "rcp60" & ID %in% unique(tree.id.numeric))
future.clim.subset.85 <- clim.ts.df %>% filter(rcp %in% "rcp85" & ID %in% unique(tree.id.numeric))
rm(clim.ts, clim.ts.df, model.1, model.2)

# generate forecasts:
# get a matrix of all samples from the betas
covariates.m$sample <- rep(1:3000, length(unique(covariates.m$variable)))
covariates.m$rowid <- 1:length(covariates.m$variable)

# get a datafraem with all of the samples from each of the 
betas.all <- covariates.m %>% select(-rowid) %>% group_by(sample) %>% spread(variable, value) %>% ungroup()%>% select( -sample)
head(betas.all)

# get a matrix of all the tree-level random effects
alpha.tree <- alpha_trees[,1]
head(alpha.tree)

# create a list of the covariates


# start with the samples of X, these are saved in x.pred from our last model run
x.pred[,1] # need to select the X value

# iterate_statespace.inc is a function that takes the posterior estimates through betas.all and alpha.tree and generates
# samples of diameter increment
iterate_statespace.inc <- function( x = x.pred[,1],  betas.all, alpha.tree,  SDdbh = 0, SDinc = 0, covariates) {
  
  # change this equation to take the betas.all matrix and a matrix of covariates so it can be applied flexibly
  # across all models
  if(length(betas.all) > 1){
    # if the model has more than just the tree size fixed effect
    if(nrow(covariates) > 1){ # if there are mulitple values for each future covariate
      
      # predict tree growth from tree values and posterior covariates
      tree.growth.list <- lapply(X = 1:nrow(covariates), FUN = function(p){alpha.tree +# sampled from tree level alpha randome effect
        betas.all[,1]*(x) + # Tree Size effect
        as.matrix(betas.all[,2:length(betas.all)]) %*% as.numeric(covariates[p,])})  # multiply by the rest of the covariate values
      tree.growth <- do.call(rbind, tree.growth.list)
      increment <- apply(tree.growth, 1, function(n){rlnorm(n = 1, meanlog = n, sdlog = SDinc)})
      #increment <- do.call(rbind, lapply(tree.growth, function(n){rlnorm(n = 1, meanlog = n, sdlog = SDinc)}))
      
      
      }else{ 
        #if(nrow(betas.all)>1)
        if(nrow(betas.all)>1){
        # if covariates are just a single/mean value, but betas have uncertainty
        
    # predict tree growth from tree values and posterior covariates
    tree.growth <- alpha.tree +# sampled from tree level alpha randome effect
      betas.all[,1]*(x) + # Tree Size effect
      as.matrix(betas.all[,2:length(betas.all)]) %*% as.numeric(covariates) # multiply by the rest of the covariate values
    increment <- apply(tree.growth, 1, function(n){rlnorm(n = 1, meanlog = n, sdlog = SDinc)})
    
     }else{ # if we are using beta means
          
          # predict tree growth from tree values and posterior covariates
          tree.growth <- alpha.tree +# sampled from tree level alpha randome effect
            betas.all[,1]*(x) + # Tree Size effect
            sum(betas.all[,2:length(betas.all)] *covariates)  # multiply by the rest of the covariate values
          increment <- do.call(rbind, lapply(tree.growth, function(n){rlnorm(n = 1, meanlog = n, sdlog = SDinc)}))
          
        }

    
    }}else{ # if the model just as one single value
  # predict tree growth from tree values and poseterior covariates
  tree.growth <- alpha.tree +# sampled from tree level alpha randome effect
    betas.all[,1]*(x) # Tree Size effect
  
  increment <- apply(tree.growth, 1, function(n){rlnorm(n = 1, meanlog = n, sdlog = SDinc)})
  
    }
  
  
  #tree.growth <-  ifelse(tree.growth < 0, 0, tree.growth) # we should actually sample from lognormal distribution here
 
  #xpred <- rnorm(length(tree.growth), (tree.growth + x), SDdbh) 
  
  increment
  
}

# take this function and run through it to make forecasts for the next century for the following 50 trees:



#---------------------------------------------------------------------------
##  Initial condition uncertainty: make forecasts from all MCMC iterations of
##    the final year, but use mean parameter values, mean climate, and no process error.
#---------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------
# Forecast function to walk tree diameters forward
#-------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------
##  Including climate model uncertainty: make forecasts from all MCMC iterations of
##    the final year, but use mean parameter values, GCM climate samples, and no process error.
#---------------------------------------------------------------------------
# set up the data frame to match trees and the alpha numbers
ID.match <- data.frame(ID = tree.id.numeric, 
           alpha.id = 1:length(unique(tree.id.numeric)))


forecast.tree <- function(idx){
    # now filter for the full ensemble
     m <- as.vector(data.frame(ID.match %>% filter(alpha.id == idx) %>% select(ID))$ID)
     future.ens <- future.clim.subset.26 %>% filter(ID == m)
   
    ens.proj.ordered <-  future.ens[order( future.ens$year),]
    
    ens.proj.ordered <-  future.ens[order( future.ens$year),]
    
    ggplot(ens.proj.ordered, aes(x = as.numeric(year), y = tmax.scaled, group=modelrun))+geom_line()
    unique(ens.proj.ordered$modelrun)
    # set up the output matrices
    time_steps <- length(2018:2098)
    nMCMC <- length(x.pred[,"x[1,1]"])
    forecast <- matrix(data = NA, nrow = nMCMC*21, ncol = time_steps)
    dbh <- matrix(data = NA, nrow = nMCMC*21, ncol = time_steps)
    
    # just for 1 tree # note that we start at x = 53
    for(t in 1:time_steps){
      if(t == 1){
        inc.pred <- iterate_statespace.inc(x = x.pred[,paste0("x[", idx,",", 53,"]")],  betas.all = betas.all, alpha.tree = alpha_trees[,idx], 
                                           SDdbh = median(sigma$sigma_dbh), 
                                           SDinc = median(sigma$sigma_inc),
                                           covariates = data.frame(TMAX = as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter( year %in% (2017 + t) &  rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(tmax.scaled)))), 
                                                                   PPT = as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(  year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale))%>% distinct() %>% select(ppt.scale))), 
                                                                   SDI = rep(model.data$SDI[idx], length(as.vector(data.matrix(ens.proj.ordered %>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct()%>% select(ppt.scale)))))))
        forecast[,t] <- inc.pred
        dbh[,t] <- forecast[,t]+x.pred[,paste0("x[", idx,",", 53,"]")]
        
      }else{
        
        inc.pred <- iterate_statespace.inc(x = dbh[,t-1],  betas.all = betas.all, alpha.tree = alpha_trees[,idx], 
                                           SDdbh = median(sigma$sigma_dbh), 
                                           SDinc = median(sigma$sigma_inc),
                                           covariates = data.frame(TMAX = as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(year %in% (2017 + t) &  rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(tmax.scaled)))), 
                                                                   PPT = as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(ppt.scale))), 
                                                                   SDI = rep(model.data$SDI[idx], length(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(ppt.scale)))))))
        
        forecast[,t] <- inc.pred
        dbh[,t] <- forecast[,t]+dbh[,t-1]
        
      }
      
      
    }
    
    
    varianceIC <- apply(forecast,2,var) # get variance from IC
    forecast.ic <- apply(forecast, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
    forecast.m <- reshape2::melt(forecast.ic) %>% group_by(Var2) %>% spread(Var1, value)
    
    inc.gg <- ggplot()+geom_line(data = forecast.m, aes(x = Var2, y = `50%`))+
      geom_ribbon(data = forecast.m, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
      ylab("forecasted diameter increment")+theme_bw()+xlab("Years after 2018")
    
    var.inc.ic <- apply(dbh, 2, var)
    dbh.ic <- apply(dbh, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
    dbh.forecast.m <- reshape2::melt(dbh.ic) %>% group_by(Var2) %>% spread(Var1, value)

    dbh.gg <- ggplot()+geom_line(data = dbh.forecast.m, aes(x = Var2, y = `50%`))+
      geom_ribbon(data = dbh.forecast.m, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
      ylab("forecasted diameter")+theme_bw()+xlab("Years after 2018")

    #png(height = 4, width = 8, units = "in", res = 150, paste0("outputs/increment_forecast_",m, ".png"))
    cowplot::plot_grid(dbh.gg,inc.gg, align = "hv", ncol = 2)
    ggsave(height = 4, width = 8, units = "in", paste0("outputs/increment_forecast_",m, ".png"))
    #dev.off()
    
}

# geneate the forecast for one tree

# not all trees in this dataset have matching future climate here:
# Get the unique tree.ids that do match so we can make forecasts of them
forecast.tree.id <- unique(ID.match %>% filter(ID %in% unique (future.clim.subset.26$ID)) %>% select(alpha.id))

forecast.tree(idx = forecast.tree.id[1,])
forecast.tree(idx = forecast.tree.id[2,])

# run the forecasts for all the trees--This will save images in the outputs folder
#lapply(forecast.tree.id[,1], forecast.tree)



#-------------------------------------------------------------------------------------------
# Forecast with uncertainty partitioning
#-------------------------------------------------------------------------------------------
# This function makes 4 forecasts, adding a new type of uncertainty in each new forecast--
# Initial condition uncertainty (uncertainty in tree size at 2018), 
# Parameter Uncertainty (associated with fixed effects & random effects)
forecast.tree.uncertainty.partitioning <- function(idx, future.climate = future.clim.subset.26){
  
  m <- as.vector(data.frame(ID.match %>% filter(alpha.id == idx) %>% select(ID))$ID)
  
   # now filter for the full ensemble
  future.ens <- future.climate %>% filter(ID == m)
  ens.proj.ordered <-  future.ens[order( future.ens$year),]
  
  ens.proj.means <- ens.proj.ordered %>% group_by(year) %>% summarise(meantmax = mean(tmax.scaled))
  
  #------------------------------------------------
  # Tree Size (Initial conditions uncertainty) only
  time_steps <- length(2018:2099)
  nMCMC <- length(x.pred[,"x[1,1]"])
  forecast <- matrix(data = NA, nrow = nMCMC, ncol = time_steps)
  dbh <- matrix(data = NA, nrow = nMCMC, ncol = time_steps)
  betas.point <- matrix(colMeans(betas.all), nrow = 1, ncol = length(colMeans(betas.all)))
  alphas.point <- matrix(colMeans(alpha_trees), nrow = 1, ncol = length(colMeans(alpha_trees)))
  
  # just for 1 tree # note that we start at x = 53
  for(t in 1:time_steps){
    if(t == 1){
      inc.pred <- iterate_statespace.inc(x = x.pred[,paste0("x[", idx,",", 53,"]")],  
                                         betas.all = betas.point, 
                                         alpha.tree = alphas.point[,idx], 
                                         SDdbh = 0, 
                                         SDinc = 0,
                                         #SDadd = 0,
                                         covariates = data.frame(TMAX = mean(as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(tmax.scaled)))), na.rm =TRUE), 
                                                                 PPT= mean(as.vector(data.matrix(ens.proj.ordered %>% ungroup() %>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(ppt.scale))), na.rm =TRUE), 
                                                                 SDI = model.data$SDI[idx]))#covariates[t,])
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+x.pred[,paste0("x[", idx,",", 53,"]")]
      
    }else{
      
      inc.pred <- iterate_statespace.inc(x = dbh[,t-1],  betas.all = betas.point, alpha.tree = alphas.point[,idx], 
                                         SDdbh = 0, 
                                         SDinc = 0,
                                         #SDadd = 0,
                                         covariates = data.frame(TMAX = mean(as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(tmax.scaled)))), na.rm =TRUE), 
                                                                 PPT= mean(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(ppt.scale))), na.rm =TRUE), 
                                                                 SDI = model.data$SDI[idx]))
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+dbh[,t-1]
      
    }
    
    
  }
  
  
  varianceIC <- apply(forecast,2,var) # get variance from IC
  forecast.ic <- apply(forecast, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
  forecast.m <- reshape2::melt(forecast.ic) %>% group_by(Var2) %>% spread(Var1, value)
  
  ggplot()+geom_line(data = forecast.m, aes(x = Var2, y = `50%`))+
    geom_ribbon(data = forecast.m, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
    ylab("forecasted diameter increment")+theme_bw()+xlab("Years after 2018")
  
  var.dbh.ic <- apply(dbh, 2, var)
  dbh.ic <- apply(dbh, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
  dbh.forecast.m <- reshape2::melt(dbh.ic) %>% group_by(Var2) %>% spread(Var1, value)
  
  ggplot()+geom_line(data = dbh.forecast.m, aes(x = Var2, y = `50%`))+
    geom_ribbon(data = dbh.forecast.m, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
    ylab("forecasted diameter")+theme_bw()+xlab("Years after 2018")
  
  #---------------------------------------
  # Parameter Uncertainty--Uncertainty around fixed effects 
  #---------------------------------------
  time_steps <- length(2018:2099)
  nMCMC <- length(x.pred[,"x[1,1]"])
  forecast <- matrix(data = NA, nrow = nMCMC, ncol = time_steps)
  dbh <- matrix(data = NA, nrow = nMCMC, ncol = time_steps)
  
  # just for 1 tree # note that we start at x = 53
  for(t in 1:time_steps){
    if(t == 1){
      inc.pred <- iterate_statespace.inc(x = x.pred[,paste0("x[", idx,",", 53,"]")],  
                                         betas.all = betas.all, 
                                         alpha.tree = alpha_trees[,idx], 
                                         SDdbh = 0, 
                                         SDinc = 0,
                                         #SDadd = 0,
                                         covariates = data.frame(TMAX = mean(as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(tmax.scaled)))), na.rm =TRUE), 
                                                                 PPT= mean(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(ppt.scale))), na.rm =TRUE), 
                                                                 SDI = model.data$SDI[idx]))#covariates[t,])
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+x.pred[,paste0("x[", idx,",", 53,"]")]
      
    }else{
      
      inc.pred <- iterate_statespace.inc(x = dbh[,t-1],  betas.all = betas.all, alpha.tree = alpha_trees[,idx], 
                                         SDdbh = 0, 
                                         SDinc = 0,
                                         #SDadd = 0,
                                         covariates = data.frame(TMAX = mean(as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(tmax.scaled)))), na.rm =TRUE), 
                                                                 PPT= mean(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(ppt.scale))), na.rm =TRUE), 
                                                                 SDI = model.data$SDI[idx]))
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+dbh[,t-1]
      
    }
    
    
  }
  
  
  variancePARAM <- apply(forecast,2,var) # get variance from IC
  forecast.PARAM <- apply(forecast, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
  forecastPARAM.m <- reshape2::melt(forecast.PARAM) %>% group_by(Var2) %>% spread(Var1, value)
  
  ggplot()+geom_line(data = forecastPARAM.m, aes(x = Var2, y = `50%`))+
    geom_ribbon(data = forecast.m, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
    ylab("forecasted diameter increment")+theme_bw()+xlab("Years after 2018")
  
  var.dbh.PARAM <- apply(dbh, 2, var)
  dbh.PARAM <- apply(dbh, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
  dbh.forecastPARAM.m <- reshape2::melt(dbh.PARAM) %>% group_by(Var2) %>% spread(Var1, value)
  
  ggplot()+geom_line(data = dbh.forecastPARAM.m, aes(x = Var2, y = `50%`))+
    geom_ribbon(data = dbh.forecast.m, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
    ylab("forecasted diameter")+theme_bw()+xlab("Years after 2018")
  
  #---------------------------------------
  # Process uncertainty (process errors)
  #---------------------------------------
  time_steps <- length(2018:2099)
  nMCMC <- length(x.pred[,"x[1,1]"])
  forecast <- matrix(data = NA, nrow = nMCMC, ncol = time_steps)
  dbh <- matrix(data = NA, nrow = nMCMC, ncol = time_steps)
  
  # just for 1 tree # note that we start at x = 53
  for(t in 1:time_steps){
    if(t == 1){
      inc.pred <- iterate_statespace.inc(x = x.pred[,paste0("x[", idx,",", 53,"]")],  betas.all = betas.all, alpha.tree = alpha_trees[,idx], 
                                         SDdbh = 0, 
                                         SDinc = mean(sigma[,"sigma_add"]),
                                         #SDadd = 0,
                                         covariates = data.frame(TMAX = mean(as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(tmax.scaled)))), na.rm =TRUE), 
                                                                 PPT= mean(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(ppt.scale))), na.rm =TRUE), 
                                                                 SDI = model.data$SDI[idx]))#covariates[t,])
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+x.pred[,paste0("x[", idx,",", 53,"]")]
      
    }else{
      
      inc.pred <- iterate_statespace.inc(x = dbh[,t-1],  betas.all = betas.all, alpha.tree = alpha_trees[,idx], 
                                         SDdbh = 0, 
                                         SDinc = mean(sigma[,"sigma_add"]),
                                         #SDadd = 0,
                                         covariates = data.frame(TMAX = mean(as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(tmax.scaled)))), na.rm =TRUE), 
                                                                 PPT= mean(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(ppt.scale))), na.rm =TRUE), 
                                                                 SDI = model.data$SDI[idx]))
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+dbh[,t-1]
      
    }
    
    
  }
  
  
  variancePROCESS <- apply(forecast,2,var) # get variance from IC
  forecast.PROCESS <- apply(forecast, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
  forecastPROCESS.m <- reshape2::melt(forecast.PROCESS) %>% group_by(Var2) %>% spread(Var1, value)
  
  ggplot()+geom_line(data = forecastPROCESS.m, aes(x = Var2, y = `50%`))+
    geom_ribbon(data = forecast.m, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
    ylab("forecasted diameter increment")+theme_bw()+xlab("Years after 2018")
  
  var.dbh.PROCESS <- apply(dbh, 2, var)
  dbh.PROCESS <- apply(dbh, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
  dbh.forecastPROCESS.m <- reshape2::melt(dbh.PROCESS) %>% group_by(Var2) %>% spread(Var1, value)
  
  ggplot()+geom_line(data = dbh.forecastPROCESS.m, aes(x = Var2, y = `50%`))+
    geom_ribbon(data = dbh.forecast.m, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
    ylab("forecasted diameter")+theme_bw()+xlab("Years after 2018")
  
  #---------------------------------------
  # Driver uncertainty (future climate driver uncertainty)
  #---------------------------------------
  # Since we have been adding in uncertainty this is the Full Forecast:
  future.ens <- future.clim.subset.26 %>% filter(ID == m)
  ens.proj.ordered <-  future.ens[order( future.ens$year),]
  
  ens.proj.ordered <-  future.ens[order( future.ens$year),]
  
  ggplot(ens.proj.ordered, aes(x = as.numeric(year), y = tmax.scaled, group=modelrun))+geom_line()
  unique(ens.proj.ordered$modelrun)
  # set up the output matrices
  time_steps <- length(2018:2098)
  nMCMC <- length(x.pred[,"x[1,1]"])
  forecast <- matrix(data = NA, nrow = nMCMC*21, ncol = time_steps)
  dbh <- matrix(data = NA, nrow = nMCMC*21, ncol = time_steps)
  
  # just for 1 tree # note that we start at x = 53
  for(t in 1:time_steps){
    if(t == 1){
      inc.pred <- iterate_statespace.inc(x = x.pred[,paste0("x[", idx,",", 53,"]")],  
                                         betas.all = betas.all, 
                                         alpha.tree = alpha_trees[,idx], 
                                         SDdbh = 0, 
                                         SDinc = 0,
                                         #SDadd = median(sigma$sigma_add),
                                         covariates = data.frame(TMAX = as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter( year %in% (2017 + t) &  rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(tmax.scaled)))), 
                                                                 PPT = as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(  year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale))%>% distinct() %>% select(ppt.scale))), 
                                                                 SDI = rep(model.data$SDI[idx], length(as.vector(data.matrix(ens.proj.ordered %>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct()%>% select(ppt.scale)))))))
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+x.pred[,paste0("x[", idx,",", 53,"]")]
      
    }else{
      
      inc.pred <- iterate_statespace.inc(x = dbh[,t-1],  
                                         betas.all = betas.all, 
                                         alpha.tree = alpha_trees[,idx], 
                                         SDdbh = 0, 
                                         SDinc = 0,
                                         #SDadd = median(sigma$sigma_add),
                                         covariates = data.frame(TMAX = as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(year %in% (2017 + t) &  rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(tmax.scaled)))), 
                                                                 PPT = as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(ppt.scale))), 
                                                                 SDI = rep(model.data$SDI[idx], length(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(ppt.scale)))))))
      
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+dbh[,t-1]
      
    }
    
    
  }
  
  Varianceclimate <- apply(forecast,2,var) # get variance from IC
  forecast.clim <- apply(forecast, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
  forecast.m.clim <- reshape2::melt(forecast.clim) %>% group_by(Var2) %>% spread(Var1, value)
  
  inc.gg <- ggplot()+geom_line(data = forecast.m.clim, aes(x = Var2, y = `50%`))+
    geom_ribbon(data = forecast.m.clim, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
    ylab("forecasted diameter increment")+theme_bw()+xlab("Years after 2018")
  inc.gg 
  
  var.dbh.clim <- apply(dbh, 2, var)
  dbh.clim <- apply(dbh, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
  dbh.forecast.m.clim<- reshape2::melt(dbh.clim) %>% group_by(Var2) %>% spread(Var1, value)
  
  dbh.gg <- ggplot()+geom_line(data = dbh.forecast.m.clim, aes(x = Var2, y = `50%`))+
    geom_ribbon(data = dbh.forecast.m.clim, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
    ylab("forecasted diameter")+theme_bw()+xlab("Years after 2018")
  dbh.gg 
  
  #---------------------------------------
  # Process and Driver Uncertainty!
  #---------------------------------------
  # Since we have been adding in uncertainty this is the Full Forecast:
  future.ens <- future.clim.subset.26 %>% filter(ID == m)
  ens.proj.ordered <-  future.ens[order( future.ens$year),]
  
  ens.proj.ordered <-  future.ens[order( future.ens$year),]
  
  ggplot(ens.proj.ordered, aes(x = as.numeric(year), y = tmax.scaled, group=modelrun))+geom_line()
  unique(ens.proj.ordered$modelrun)
  # set up the output matrices
  time_steps <- length(2018:2098)
  nMCMC <- length(x.pred[,"x[1,1]"])
  forecast <- matrix(data = NA, nrow = nMCMC*21, ncol = time_steps)
  dbh <- matrix(data = NA, nrow = nMCMC*21, ncol = time_steps)
  
  # just for 1 tree # note that we start at x = 53
  for(t in 1:time_steps){
    if(t == 1){
      inc.pred <- iterate_statespace.inc(x = x.pred[,paste0("x[", idx,",", 53,"]")],  
                                         betas.all = betas.all, 
                                         alpha.tree = alpha_trees[,idx], 
                                         SDdbh = 0, 
                                         SDinc = median(sigma$sigma_add),
                                         #SDadd = median(sigma$sigma_add),
                                         covariates = data.frame(TMAX = as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter( year %in% (2017 + t) &  rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(tmax.scaled)))), 
                                                                 PPT = as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(  year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale))%>% distinct() %>% select(ppt.scale))), 
                                                                 SDI = rep(model.data$SDI[idx], length(as.vector(data.matrix(ens.proj.ordered %>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct()%>% select(ppt.scale)))))))
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+x.pred[,paste0("x[", idx,",", 53,"]")]
      
    }else{
      
      inc.pred <- iterate_statespace.inc(x = dbh[,t-1],  
                                         betas.all = betas.all, 
                                         alpha.tree = alpha_trees[,idx], 
                                         SDdbh = 0, 
                                         SDinc = median(sigma$sigma_add),
                                         #SDadd = median(sigma$sigma_add),
                                         covariates = data.frame(TMAX = as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(year %in% (2017 + t) &  rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(tmax.scaled)))), 
                                                                 PPT = as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(ppt.scale))), 
                                                                 SDI = rep(model.data$SDI[idx], length(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(ppt.scale)))))))
      
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+dbh[,t-1]
      
    }
    
    
  }
  
  
  VariancePROCESS <- apply(forecast,2,var) # get variance from IC
  forecast.PROCESS <- apply(forecast, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
  forecast.m.PROCESS <- reshape2::melt(forecast.PROCESS) %>% group_by(Var2) %>% spread(Var1, value)
  
  inc.gg <- ggplot()+geom_line(data = forecast.m.PROCESS, aes(x = Var2, y = `50%`))+
    geom_ribbon(data = forecast.m.clim, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
    ylab("forecasted diameter increment")+theme_bw()+xlab("Years after 2018")
  inc.gg 
  
  var.dbh.PROCESS <- apply(dbh, 2, var)
  dbh.PROCESS <- apply(dbh, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
  dbh.forecastPROCESS.m <- reshape2::melt(dbh.clim) %>% group_by(Var2) %>% spread(Var1, value)
  
  dbh.gg <- ggplot()+geom_line(data = dbh.forecastPROCESS.m, aes(x = Var2, y = `50%`))+
    geom_ribbon(data = dbh.forecastPROCESS.m, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
    ylab("forecasted diameter")+theme_bw()+xlab("Years after 2018")
  
  dbh.gg
  
  
  
  #-------------------------------------------------------------------------------------------
  # Combine all the forecasts and make some partitioning plots
  #-------------------------------------------------------------------------------------------
  # make all the forecast and uncertainty partitioning 
  

  
    V.pred.sim.dbh     <- rbind(var.dbh.clim[1:81], var.dbh.PROCESS[1:81], var.dbh.PARAM[1:81], var.dbh.ic[1:81])
    
    # combine forecasts:
    pred.sims.dbh     <- data.frame(IPP.0 = dbh.clim[1, 1:81],
                                IPD.0 = dbh.PROCESS[1, 1:81],
                                IPA.0 = dbh.PARAM[1, 1:81],
                                IP.0 = dbh.ic[1, 1:81],
                                
                                IPP.100 = dbh.clim[4, 1:81],
                                IPD.100 = dbh.PROCESS[4, 1:81],
                                IPA.100 = dbh.PARAM[4, 1:81],
                                IP.100 = dbh.ic[4, 1:81],
                              
                                year = 2018:2098)
    
   
    

    V.pred.sim.inc     <- rbind(Varianceclimate[1:81], 
                                variancePROCESS[1:81], 
                                variancePARAM[1:81], 
                                varianceIC[1:81])
    
    
    pred.sims.inc     <- data.frame(IPP.0 = forecast.clim[1,1:81],
                                IPD.0 = forecast.PARAM[1,1:81],
                                IPA.0 = forecast.PROCESS[1,1:81],
                                IP.0 = forecast.ic[1,1:81],
                                
                                IPP.100 = forecast.clim[4,1:81],
                                IPD.100 = forecast.PARAM[4,1:81],
                                IPA.100 = forecast.PROCESS[4,1:81],
                                IP.100 = forecast.ic[4,1:81],
                                
                                year = 2018:2098)
    #axis.name <- "Increment"
 
    #-------------------------------------------------------------------
    # Make uncertainty partitioning plots 
    #-------------------------------------------------------------------
    
    ## DBH plots
  V.pred.sim.rel.inc <- apply(V.pred.sim.inc, 2, function(x) {x/max(x)})
  
  V.pred.sim.rel.dbh <- apply(V.pred.sim.dbh, 2, function(x) {x/max(x)})
  
 
  pred.sims.m.dbh <- reshape2::melt(pred.sims.dbh, id.vars = "year")
  pred.sims.class <- pred.sims.m.dbh %>% separate(col = variable, sep = "[.]", into = c("unc","lo")) %>%
    spread(key = lo, value = value)
  colnames(pred.sims.class) <- c("year", "uncertainty", "Low", "High")
  my_cols <- c( "#d95f02",
                         "#1b9e77",
                       
                        "black",
                        "#7570b3", 
                        "grey")
                        
  pred.sims.class$uncertainty <- factor(pred.sims.class$uncertainty, levels = c("IPP","IPD","IPA", "IP"))
  
  
  predY_plot_dbh <- ggplot(data = pred.sims.class, aes(x=year, fill = uncertainty))+
    geom_ribbon(aes(ymin=Low, ymax=High))+
    ylab("Forecasted Diameter (cm)")+
    xlab("Year")+theme_bw()+
    scale_fill_manual(values = my_cols, name = NULL)+ theme(legend.position = "none", panel.grid = element_blank())#+ylim(-10, 1)
  
  predY_plot_dbh
  
  ##--------------------------------------------------------------
  #  Plot the proportion of uncertainty
  #--------------------------------------------------------------
  var_rel_preds <- as.data.frame(t(V.pred.sim.rel.dbh*100))
  var_rel_preds$x <- 1:nrow(var_rel_preds)
  
  my_cols <- c(  "#d95f02",
                          "#1b9e77",
                      
                        "black",
                        "#7570b3", 
                        "grey")
  tmpvar <- var_rel_preds
  tmpvar$year <- 2018:2098
  colnames(tmpvar) <- c( "Process Error","Driver Uncertainty","Parameter Uncertainty", "Initial Conditions", "x", "year")
  variance.df <- tmpvar %>% gather(simtype, variance, -x, -year)
                        
  variance.df$simtype <- factor(variance.df$simtype, levels = c("Driver Uncertainty","Process Error","Parameter Uncertainty", "Initial Conditions"))
                        
  prop.var.dbh <- ggplot(variance.df, aes(x=year, fill = simtype))+
                     geom_ribbon(aes(ymin=0, ymax=variance), color = "grey")+
                     ylab(paste("% of total variance diameter"))+    xlab("Year")+
                     scale_fill_manual(values = my_cols, name = NULL)+#, 
                          #labels = c("Process error", "Driver uncertainty",  "Parameter uncertainty","alpha uncertainty", "Initial conditions"))+
                     scale_y_continuous(labels=paste0(seq(0,100,25),"%"), expand = c(0, 0))+
                     theme_bw()+theme(panel.grid = element_blank())
  prop.var.dbh  
  
  
  #----------------------------------------
  # Make the same plots but for increments
  #----------------------------------------
  ## INCREMENT plots
  V.pred.sim.rel.inc <- apply(V.pred.sim.inc, 2, function(x) {x/max(x)})
  
  pred.sims.m.inc <- reshape2::melt(pred.sims.inc, id.vars = "year")
  pred.sims.class <- pred.sims.m.inc %>% separate(col = variable, sep = "[.]", into = c("unc","lo")) %>%
    spread(key = lo, value = value)
  colnames(pred.sims.class) <- c("year", "uncertainty", "Low", "High")
  my_cols <- c( "#d95f02",
                         "#1b9e77",
                       
                        "black",
                        "#7570b3", 
                        "grey")
                        
  pred.sims.class$uncertainty <- factor(pred.sims.class$uncertainty, levels = c("IPA","IPP","IPD", "IP"))
  
  
  predY_plot <- ggplot(data = pred.sims.class, aes(x=year, fill = uncertainty))+
    geom_ribbon(aes(ymin=Low, ymax=High))+
    ylab("Forecasted Increment (cm)")+
    xlab("Year")+theme_bw()+
    scale_fill_manual(values = my_cols, name = NULL)+ theme(legend.position = "none", panel.grid = element_blank())#+ylim(-10, 1)
  
  predY_plot
  
  ##--------------------------------------------------------------
  #  Plot the proportion of uncertainty
  #--------------------------------------------------------------
  var_rel_preds <- as.data.frame(t(V.pred.sim.rel.inc*100))
  var_rel_preds$x <- 1:nrow(var_rel_preds)
  
  my_cols <- c( "#d95f02",
                         "#1b9e77",
                       
                        "black",
                        "#7570b3", 
                        "grey")
                        tmpvar <- var_rel_preds
                        tmpvar$year <- 2018:2098
                        colnames(tmpvar) <- c("Process Error", "Driver Uncertainty","Parameter Uncertainty", "Initial Conditions", "x", "year")
                        variance.df <- tmpvar %>% gather(simtype, variance, -x, -year)
                        
                        variance.df$simtype <- factor(variance.df$simtype, levels = c("Driver Uncertainty","Process Error","Parameter Uncertainty","Initial Conditions"))
                        
                        prop.var <- ggplot(variance.df, aes(x=year, fill = simtype))+
                          geom_ribbon(aes(ymin=0, ymax=variance), color = "white")+
                          ylab(paste("% of total variance increment"))+    xlab("Year")+
                          scale_fill_manual(values = my_cols, name = NULL)+#, 
                          #labels = c("Process error", "Driver uncertainty",  "Parameter uncertainty","alpha uncertainty", "Initial conditions"))+
                          scale_y_continuous(labels=paste0(seq(0,100,25),"%"), expand = c(0, 0))+
                          theme_bw()+theme(panel.grid = element_blank())
                        prop.var  
                        

cowplot::plot_grid(predY_plot, prop.var, 
          predY_plot_dbh, prop.var.dbh, 
          ncol = 2, align = "hv", rel_widths = c(0.75, 1))
ggsave(height = 6, width = 8, units = "in", filename=paste0("outputs/variance_partitioning_tree_",idx,".jpg") ) 
}
#forecast.tree.uncertainty.partitioning(idx=2, future.climate = future.clim.subset.26)
# idx = 

forecast.tree.uncertainty.partitioning(idx=forecast.tree.id[1,], future.climate = future.clim.subset.26)


