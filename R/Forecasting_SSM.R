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

# driver uncertainty:
# for each scenario, in each year, for each plot, we sample from the distribution of ppt and tmax
# that has a mean of the ensemble mean, and variance of the ensemble variance

# read the rds with climate scenarios:
clim.data <- readRDS("/Users/kellyheilman/Downloads/PRISM_non_scaled.rds")

clim.ts <- readRDS("/Users/kellyheilman/Downloads/pipo.cores.ds.mean.correct.climate_2018_2099.RDS")
colnames(clim.ts)[6:7] <-c("year.ppt", "tmax.fall.spr") 


clim.ts.df <- clim.ts #$future.climate.ts
clim.ts.df$tmax.fall.spr[is.nan(clim.ts.df$tmax.fall.spr)] <- NA
#tmax.fallspr.df <- tmax.fallspr


# need to scale future climate data on the same scale as the past climate
clim.ts.df$ppt.scale <-(clim.ts.df$year.ppt-mean(as.matrix(clim.data$wintP.wateryr)))/sd(as.matrix(clim.data$wintP.wateryr))
clim.ts.df$tmax.scaled <-(clim.ts.df$tmax.fall.spr-mean(as.matrix(clim.data$tmax.fallspr)))/sd(as.matrix(clim.data$tmax.fallspr))


# need to fix this so that it aligns with trees, not plots
unique.ll <- unique(clim.ts.df[,c("lat", "lon")])
unique.ll$ID<- 1:length(unique.ll$lat)

clim.ts.df <- left_join(clim.ts.df, unique.ll)

# filter out for just the trees selected in this tutorial
tree.id %in% unique(future.clim$ID)
head(clim.ts.df)

tree.id <- substring(as.character(tree.id), 2, nchar(as.character(tree.id))) # use the tree id from the first script to select ID in future.clim.subset to scale

# set up the 4 different emissions scenarios data
future.clim.subset.26 <- clim.ts.df %>% filter(rcp %in% "rcp26" & ID %in% unique(tree.id))
future.clim.subset.45 <- clim.ts.df %>% filter(rcp %in% "rcp45" & ID %in% unique(tree.id))
future.clim.subset.60 <- clim.ts.df %>% filter(rcp %in% "rcp60" & ID %in% unique(tree.id))
future.clim.subset.85 <- clim.ts.df %>% filter(rcp %in% "rcp85" & ID %in% unique(tree.id))




# generate forecasts:

# get a matrix of all samples from the betas
covariates.m$sample <- rep(1:3000, length(unique(covariates.m$variable)))
betas.all <- covariates.m %>% group_by(sample) %>% spread(variable, value) %>% ungroup()%>% select( -sample)



# get a matrix of all the tree-level random effects
alpha.tree <- alpha_trees[,1]

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
      tree.growth.list <- lapply(X = 1:nrow(covariates), FUN = function(x){alpha.tree +# sampled from tree level alpha randome effect
          betas.all[,1]*(x) + # Tree Size effect
          as.matrix(betas.all[,2:length(betas.all)]) %*% as.numeric(covariates[x,])})  # multiply by the rest of the covariate values
      tree.growth <- do.call(rbind, tree.growth.list)
      
    }else{ # if its just a single/mean value
      
      # predict tree growth from tree values and posterior covariates
      tree.growth <- alpha.tree +# sampled from tree level alpha randome effect
        betas.all[,1]*(x) + # Tree Size effect
        as.matrix(betas.all[,2:length(betas.all)]) %*% as.numeric(covariates)  # multiply by the rest of the covariate values
      
      
    }}else{ # if the model just as one single value
      # predict tree growth from tree values and poseterior covariates
      tree.growth <- alpha.tree +# sampled from tree level alpha randome effect
        betas.all[,1]*(x) # Tree Size effect
    }
  #tree.growth <-  ifelse(tree.growth < 0, 0, tree.growth) # we should actually sample from lognormal distribution here
  increment <- apply(tree.growth, 1, function(x){rlnorm(n = 1, meanlog = x, sdlog = SDinc)})
  
  #xpred <- rnorm(length(tree.growth), (tree.growth + x), SDdbh) 
  
  increment
  
}

# take this function and run through it to make forecasts for the next century for the following 50 trees:



#---------------------------------------------------------------------------
##  Initial condition uncertainty: make forecasts from all MCMC iterations of
##    the final year, but use mean parameter values, mean climate, and no process error.
#---------------------------------------------------------------------------

# now filter for the full ensemble
m <- 1
future.ens <- future.clim.subset.26 %>% filter(ID == m)
ens.proj.ordered <-  future.ens[order( future.ens$year),]
#ggplot(proj.ordered, aes(as.numeric(year), mean.ppt))+geom_point()+stat_smooth()
# just use the ensemble means (w.out driver uncertainty)

#covariates <- data.frame(SDI,SICOND, ppt, tmax)

time_steps <- length(2018:2099)
nMCMC <- length(x.pred[,"x[1,1]"])
forecast <- matrix(data = NA, nrow = nMCMC, ncol = time_steps)
dbh <- matrix(data = NA, nrow = nMCMC, ncol = time_steps)

# just for 1 tree # note that we start at x = 53
for(t in 1:time_steps){
  if(t == 1){
    inc.pred <- iterate_statespace.inc(x = x.pred[,paste0("x[", m,",", 53,"]")],  betas.all = betas.all, alpha.tree = alpha_trees[,m], 
                                       SDdbh = 0, 
                                       covariates = data.frame(TMAX = mean(as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(tmax.scaled)))), na.rm =TRUE), 
                                                               PPT= mean(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(ppt.scale))), na.rm =TRUE), 
                                                               SDI = model.data$SDI[m]))#covariates[t,])
    forecast[,t] <- inc.pred
    dbh[,t] <- forecast[,t]+x.pred[,paste0("x[", m,",", 53,"]")]
    
  }else{
    
    inc.pred <- iterate_statespace.inc(x = dbh[,t-1],  betas.all = betas.all, alpha.tree = alpha_trees[,m], 
                                       SDdbh = 0, covariates = data.frame(TMAX = mean(as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(tmax.scaled)))), na.rm =TRUE), 
                                                                          PPT= mean(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(ID %in% m & year %in% (2017 + t) & rcp %in% "rcp26") %>% select(ppt.scale))), na.rm =TRUE), 
                                                                          SDI = model.data$SDI[m]))
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

var.inc.ic <- apply(dbh, 2, var)
dbh.ic <- apply(dbh, 2, function(x){quantile(x, c(0.025, 0.25, 0.5, 0.975), na.rm = TRUE)})
dbh.forecast.m <- reshape2::melt(dbh.ic) %>% group_by(Var2) %>% spread(Var1, value)

ggplot()+geom_line(data = dbh.forecast.m, aes(x = Var2, y = `50%`))+
  geom_ribbon(data = dbh.forecast.m, aes(x = Var2, ymin= `2.5%`, ymax = `97.5%`), fill = "brown", alpha = 0.5)+
  ylab("forecasted diameter")+theme_bw()+xlab("Years after 2018")
#-------------------------------------------------------------------------------------------
# Do the same thing but add in climate uncertainty
#-------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------
##  Including climate model uncertainty: make forecasts from all MCMC iterations of
##    the final year, but use mean parameter values, GCM climate samples, and no process error.
#---------------------------------------------------------------------------

forecast.tree <- function(m){
  # now filter for the full ensemble
  #unique(future.clim.subset.26$ID) # only have future climate for a fraction of trees??
  m.clim <- unique(future.clim.subset.26$ID)[m]
  future.ens <- future.clim.subset.26 %>% filter(ID == m.clim)
  ens.proj.ordered <-  future.ens[order( future.ens$year),]
  unique(ens.proj.ordered$modelrun)
  # set up the output matrices
  time_steps <- length(2018:2098)
  nMCMC <- length(x.pred[,"x[1,1]"])
  forecast <- matrix(data = NA, nrow = nMCMC*21, ncol = time_steps)
  dbh <- matrix(data = NA, nrow = nMCMC*21, ncol = time_steps)
  
  # just for 1 tree # note that we start at x = 53
  for(t in 1:time_steps){
    if(t == 1){
      inc.pred <- iterate_statespace.inc(x = x.pred[,paste0("x[", m,",", 53,"]")],  betas.all = betas.all, alpha.tree = alpha_trees[,m], 
                                         SDdbh = 0, 
                                         covariates = data.frame(TMAX = as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter( year %in% (2017 + t) &  rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(tmax.scaled)))), 
                                                                 PPT = as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter(  year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale))%>% distinct() %>% select(ppt.scale))), 
                                                                 SDI = rep(model.data$SDI[m], length(as.vector(data.matrix(ens.proj.ordered %>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct()%>% select(ppt.scale)))))))
      forecast[,t] <- inc.pred
      dbh[,t] <- forecast[,t]+x.pred[,paste0("x[", m,",", 53,"]")]
      
    }else{
      
      inc.pred <- iterate_statespace.inc(x = dbh[,t-1],  betas.all = betas.all, alpha.tree = alpha_trees[,m], 
                                         SDdbh = 0, 
                                         covariates = data.frame(TMAX = as.numeric(as.matrix(data.frame(ens.proj.ordered %>% ungroup()%>% filter(year %in% (2017 + t) &  rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(tmax.scaled)))), 
                                                                 PPT = as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(ppt.scale))), 
                                                                 SDI = rep(model.data$SDI[m], length(as.vector(data.matrix(ens.proj.ordered%>% ungroup() %>% filter( year %in% (2017 + t) & rcp %in% "rcp26" & !is.na(ppt.scale)) %>% distinct() %>% select(ppt.scale)))))))
      
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
  
  png(height = 4, width = 8, units = "in", res = 150, paste0("tutorial/outputs/increment_forecast_",m, ".png"))
  plot_grid(dbh.gg,inc.gg, align = "hv", ncol = 2)
  dev.off()
  
}

forecast.tree(m = 8)
lapply(unique(future.clim.subset.26$ID), forecast.tree)


## Improvements still needed here
# 1. Checking future climate & matching to the trees (non of the trees with the new model have a growth decline)
# 2. Save the forecast plots --not saving within the function
# 3. using sigma_inc for the uncertainty
# 4. Uncertainty partitioning function
# 5. Write this up so students can just call the functions and run with it for a given tree
