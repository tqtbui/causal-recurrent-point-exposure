#################################
#Simulation for estimator
#################################

#Load source code files
source("0. general functions.R")
source("1. estimator_functions.R")

#-----------------------------
#Preparing the data
#-----------------------------

#survival dataset
dat <- survival::bladder1
t_fits <- 1:max(dat$stop) #observation period

# construct the appropriate dataset for estimation from raw data
dat_ben <- dat_transf(dat = dat, 
                      t_fits = t_fits, 
                      covnames = c("number", "size"), #initial number of tumors and size of largest initial tumor
                      event.col = "rtumor", 
                      treat.col = "treatment", 
                      treat.name = c("pyridoxine", "thiotepa")) #compare placebo vs pyriodoxine and thiotepa

#-----------------------------
#Analysis
#-----------------------------

# use our method to do analysis
set.seed(123456)
fit <- BBestimators(dat_ben, t_fits = 1:60, tau = 64, kfolds = 5, 
                    verbose = TRUE, 
                    pi.library = c("SL.glm"), 
                    event.library = c("survSL.km", "survSL.coxph"), 
                    cens.library = c("survSL.km", "survSL.coxph"), 
                    c.library = c("SL.glm"),
                    d.library = c("SL.glm"), 
                    covnames = c("number", "size"))

# visualizing the results
library(dplyr)
library(ggplot2)
res <- plotBB(fit)
