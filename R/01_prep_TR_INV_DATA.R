
# load the packages used by this script
library(tidyverse)
library(dplR)
library(reshape2)
library(rstan)



#------------------------Overview------------------------------------------
# This code reads in example tree ring increment, diameter data, and covariate data for 
# running the tutorial version of the state-space model. This code is just an example of
# how to read in and format rwl/tree ring data. 

# Example tree ring data come from the Forest Inventory and Analysis Program. 
# They are from Ponderosa Pine trees in Arizona and have been previously published at: https://doi.org/10.25739/sd83-rk24
# These data are in a slightly different format from those published at the above DOI.
# We have included associated covariates (already scaled at the site or region scale) and 
# there is not any associated location information (lat-long coordinates) shared and the PLOT and treeid
# values do not match publicly available FIA database data.


#--------------------------------------------------------------------------------
# Read in the tree ring data
#--------------------------------------------------------------------------------
# if the data is in .rwl (tucson format). the function read.rwl from dplR 
rw.format <- read.rwl("data/INC_tutorial_rw.rwl")
head(rw.format) # rwl files are read in with years as rows/rownames and trees as columns

summary(rw.format) # dplR summary gives you more infomation on the ring width data for each tree
spag.plot(rw.format) # plot the extent of all the tree ring time series

# lets convert this to the format that we need:
rw.df <- data.frame(rw.format)

# add a year column since the rownames are years and the column names are tree numbers

rw.df$year <- row.names(rw.df)
head(rw.df)# this still has the same format, but now year is a column and its a dataframe

rw.df.m <- reshape2:: melt(rw.df, id.vars = "year") # using reshape2 package to "melt" the dataframe into 3 columns: year, tree, and increment measurement
colnames(rw.df.m) <- c("year", "tree", "inc") # rename the columns
head(rw.df.m) # view the first few records

# lets plot all the tree rings a few ways
hist(rw.df.m$inc, xlab = "radial increment (cm)", main = "Radial Increments")

ggplot(rw.df.m, aes(x = as.numeric(year), y = inc, group = tree))+geom_line()+xlab("Year")+
  ylab("Radial Increment (cm)")

# note that our data is already truncated at 1966 to limit the # of years we are estimating & 
# to minimize issues related to the fading record and missing competition on the forest plot
# but you could truncate your data using something like the following (make sure you do the same for DIAMTER)
rw.df.m <- rw.df.m %>% filter(year > 1965)

# our data radial increment in centimeters, so we need to convert it to diameter increment & 
# we need to reformat the data to have trees as rows and time as columns
rw.spread <- rw.df.m %>% mutate(diainc = inc*2) %>% # does the conversion to diameter increment
  mutate(year = as.numeric(year)) %>% 
  complete(year = 1966:2018, # use complete to fill in missing years from the database with NA values
           # note: the increment data only go until 1996, so we need to add in some NA values after
           
           fill = list( diainc = NA))%>%
  dplyr::select(year, tree, diainc) %>% # selects just the diameter increment, year, and tree
  group_by(tree) %>% # we will group by tree (our desired row variable)
  spread(year, diainc) # use spread to "spread" the diameter increments into columns by year

head(rw.spread)


#--------------------------------------------------------------------------------
# Read in the Diameter data
#--------------------------------------------------------------------------------
# our diameter data is in a .csv file with information about the tree id and the year the diameter was sampled
# you can use this as a guide to format your own diameter data
# Note on units: Our diameter data is already in centimeters!
DBH.data <- read.csv("data/DBH_tutorial_data.csv")
summary(DBH.data) # our data has 3 columns for tree, year sampled, and for Diameter
hist(DBH.data$DBH_cm, xlab = "Diameter (cm)")

# note that some trees only have 1 dbh measurement, while others have two
ggplot(DBH.data, aes(x = year, y = DBH_cm, group = tree))+geom_line()+geom_point()

# lets reformat these data so we have trees as rows and years as columns using spread
# note that this doesnt have years where we don't have DBH measurements
# also note that if you have additional columns in your data you will want to use `select(tree, year, DBH_cm)` to select only the year, tree, and diameter columns
DBH.data %>% group_by(tree) %>% spread(year, DBH_cm) 

# add in the missing years, so now we have a full dataframe with trees as rows, and years 1965-2018 as columns
DBH.spread <- DBH.data %>% mutate(year = as.numeric(year)) %>% 
  complete(year = 1966:2018, # use complete to fill in missing years from the database with NA values
           fill = list( DBH_cm = NA))  %>% 
  group_by(tree) %>% spread(year, DBH_cm) %>%
  filter(!is.na(tree)) # doing this added an NA value into the tree column for some reason, so remove it

summary(DBH.spread)

#--------------------------------------------------------------------------------
# Read in associated plot-level data
#--------------------------------------------------------------------------------
cov.data <- read.csv("data/tutorial_covariates.csv")
head(cov.data)
# note that these are all standardized variables already
# covariate data can be linked to the tree ring data and diametr data by treeid 
summary(cov.data)

# in the GCB model we included SDI (Stand Density Index) and SICOND (site index) as covariates, but 
# have started including Mean Annual Precipitation and Mean annual Temperature as covariates.
# we could also include other plot or tree-level variables of interest that you have associated with your data.
# You may not have many additional covariates to include either, which means you might fit a model 
# with a plot-level random effect to help explain some spatial variation in growth.

#--------------------------------------------------------------------------------
# Read in associated time-varying climate data
#--------------------------------------------------------------------------------
climate.long <- read.csv("data/tutorial_climate.csv")
head(climate.long)
# This example climate data was extracted from PRISM climate data for each FIA plot
# I have it in a dataframe that has a column for tree, time, and two of the climate variables we used
# see https://prism.oregonstate.edu for downloads

# note that these variables are already scaled

summary(climate.long)

# lets create a list of two matrices, one for precipitation and one for the temperature variable
time_data <- list()
water_yr_precip <- climate.long %>% dplyr::select(tree, year, Water_year_precip) %>%
  group_by(tree) %>% 
  spread(year, Water_year_precip) 

fall_spr_Tmax <- climate.long %>% dplyr::select(tree, year, Fall_spring_tmax) %>%
  group_by(tree) %>% 
  spread(year, Fall_spring_tmax) 
head(fall_spr_Tmax)

time_data$water_yr_precip <- water_yr_precip[, 2:length(water_yr_precip)]
time_data$fall_spr_Tmax <- fall_spr_Tmax[, 2:length(fall_spr_Tmax)]


#--------------------------------------------------------------------------------
# Assemble all the data into a model.data object that can be used by STAN
#--------------------------------------------------------------------------------
# note that this isnt the only way to format data for STAN, it depends on how you 
# specify the model. We specified the model in a way that has an index for tree and year
# so that it is easier to read. 
# STAN does not directly like to take NA values in the data, so we need to create a missing
# data model for both the increment and diameter data. 

# make it pretty small--50 trees
set.seed(22)
row.id.samp <- sample(1:518, 50, replace = FALSE)
tree.id <- rw.spread[row.id.samp, ]$tree
dat <- unique(rw.spread[row.id.samp,2:length(rw.spread)])

# save tree.id for use in the forecasting
saveRDS(tree.id, "data/tree.id.rand.sample.rds")

# STAN does not by default take missing data, so we need indexing for a missing data model
# Extract the missing values into a VECTOR
dat_complete <- dat[!is.na(dat)]

# Extract the missing and present values as MATRICES
ind_pres <- which(!is.na(dat), arr.ind = TRUE)
ind_miss <- which(is.na(dat), arr.ind = TRUE)


# get missing z data:
datz <- DBH.spread[row.id.samp,2:length(DBH.spread)]



# Extract the missing values into a VECTOR
dat_completez <- datz[!is.na(datz)]
head(dat_completez) # vector of non-missing diameters

# Extract the missing and present values as MATRICES
ind_presz <- which(!is.na(datz), arr.ind = TRUE)
ind_missz <- which(is.na(datz), arr.ind = TRUE)

head(ind_missz) # matrix of indices for the mssing values

# put everything into a list for STAN input
model.data <- list(Nrow = nrow(dat), # Nrow is the # of trees
                   Ncol = ncol(dat), # Ncol is the # of years
                   
                   # data and indices for the missing data and complete increment data
                   Ncomp = length(dat_complete),
                   Nmiss = sum(is.na(dat)),
                   dat_complete = dat_complete,
                   ind_pres = ind_pres,
                   ind_miss = ind_miss , 
                   
                   # data and indices for the missing data and complete diameter data
                   Nrow_z = nrow(datz),
                   Ncol_z = ncol(datz),
                   Ncomp_z = length(dat_completez),
                   Nmiss_z = sum(is.na(datz)),
                   dat_completez = dat_completez,
                   ind_presz = ind_presz,
                   ind_missz = ind_missz, 
                   
                   # read in the covariate data
                   tmaxAprMayJunscaled = data.frame(time_data$fall_spr_Tmax[row.id.samp,]), 
                   wateryrscaled = data.frame(time_data$water_yr_precip[row.id.samp,]),  
                   SDI = cov.data[row.id.samp,]$SDI)


model.data

# View all of the data inputs
names(model.data)

# save model.data
saveRDS(model.data, "outputs/model.data.RDS")
saveRDS(dat, "outputs/dat.rds") # save the increment data
saveRDS(datz, "outputs/datz.rds") # save the diameter data
# Next 
