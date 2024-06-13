
# load the packages used by this script
library(tidyverse)
library(dplR)
library(reshape2)
library(rstan)

#---
# note, some of the data conversions are a little strange (e.g. converting to wide and then back to long) because
# I was adapting tutorial code to fit the model structure developed for our larger project after the fact 
# and I was trying to quickly get you from the tutorial dataset to an example dataset structure.

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


# put this in long format:
y.long <- reshape2::melt(rw.spread) %>% mutate(tree.id = str_remove(tree, "X")) %>% 
  rename(`year` = "variable", 
         `inc` = "value") %>% select(tree.id, year, inc) %>% filter(!is.na(inc))
colnames(y.long) <- c("treeid", "year", "inc")
# does the conversion to diameter increment 
#y.long <- y.long %>% mutate(diainc = (inc)*2) %>% filter(!is.na(diainc))
y.long$year <- as.numeric(as.character(y.long$year))
y.long$diainc
summary(y.long$inc)

# separate y into training and testing data:
N_train <- nrow(y.long)*0.80
N_test <- nrow(y.long)*0.20
train_ind <- sample(c(1:nrow(y.long)), size = N_train, replace = FALSE)

train_y <- y.long[train_ind,]
test_y <- y.long[-train_ind,]

summary(train_y)
summary(test_y)



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
  complete(year = 1966:2010, # use complete to fill in missing years from the database with NA values
           fill = list( DBH_cm = NA))  %>% 
  group_by(tree) %>% spread(year, DBH_cm) %>%
  filter(!is.na(tree)) # doing this added an NA value into the tree column for some reason, so remove it

summary(DBH.spread)


# set up z diameter data frame in long format and leave out some diameter remeasuremetns for validation:

#row.names(data$z) <- 1:1046 # 1 to # of trees
z.long <- reshape2::melt(as.matrix(DBH.spread[,2:ncol(DBH.spread)]))

colnames(z.long) <- c("tree.id", "year", "DIA")
summary(z.long)

time.df <- data.frame(year = 1966:2010, 
                      time = 1:length(1966:2010))
z.long <- z.long %>% filter(!is.na(z.long$DIA))
class(z.long$year)
z.repeated <- z.long %>% filter(!is.na(DIA)) %>% group_by(tree.id) %>% summarise(n = n()) %>% filter(n > 1)


z.long <- left_join(z.long, time.df) %>% filter(time %in% 1:36)
train_y <- left_join(train_y, time.df)
test_y <- left_join(test_y, time.df)


# if SDI is NA for the first year, assign the next year's SDI
data$SDIscaled[is.na(data$SDIscaled[,1]), ] <- data$SDIscaled[is.na(data$SDIscaled[,1]),2]





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

mod.data <- list(Nrow = length(unique(y.long$treeid)),
                 Ncol = 36,
                
                 
                 # data and indices for increment
                 y = train_y$inc, 
                 treeinc = as.numeric(train_y$treeid), 
                 yearinc = as.numeric(train_y$time),
                 Ninc = length(train_y$inc),
                 
                 
                 # data and indices for the complete diameter data
                 z = z.long$DIA, # diameter data 
                 treedia = z.long$tree.id, # index for the tree
                 yeardia = z.long$time, # index for the year
                 Ndia = length(z.long$DIA), 
                 
                 
                 tmaxAprMayJunscaled = time_data$fall_spr_Tmax[,1:36], 
                 wateryrscaled = time_data$water_yr_precip[,1:36], 
                 
                 SDI = cov.data$SDI)

saveRDS(mod.data, "data/mod.data.long.format.rds")


