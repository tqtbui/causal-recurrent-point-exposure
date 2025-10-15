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
nsample <- length(unique(dat$id)) #number of individualss in the data

# construct the appropriate dataset for estimation from raw data
dat_ben <- matrix(0, ncol = 5+length(t_fits), nrow = nsample)
colnames(dat_ben) <- c("delta", "x", "a", "l1", "l2", 
                       paste0("NE_", t_fits))

#loop through each individual to collect data
for (i in 1:nsample) {
  
  #covariates
  tmp <- dat[dat$id == i,]
  dat_ben[i,"l1"] <- tmp$number[1] #initial number of tumors
  dat_ben[i,"l2"] <- tmp$size[1] #size of largest initial tumor
  
  #exposure
  dat_ben[i,"a"] <- ifelse(tmp$treatment[1] == "placebo", 0, 1) #compare placebo vs pyriodoxine and thiotepa
  
  #survival
  dat_ben[i,"x"] <- tmp$stop[nrow(tmp)] #end of observation period
  if (tmp$status[nrow(tmp)] == 0) {
    dat_ben[i,"delta"] <- 0
  } else if ((tmp$status[nrow(tmp)] == 1) & (tmp$stop[nrow(tmp)] < max(t_fits))) {
    dat_ben[i,"delta"] <- 0
  } else {
    dat_ben[i,"delta"] <- 1
  }
  
  #number of events
  for (j in 1:nrow(tmp)) {
    if (tmp$status[j] == 1) {
      dat_ben[i,paste0("NE_", tmp$stop[j]:max(t_fits))] <- dat_ben[i,paste0("NE_", tmp$stop[j]:max(t_fits))] + ifelse(tmp$rtumor[j] == ".", 0, as.numeric(tmp$rtumor[j]))
    }
  }
}

# remove individual that die at time 0
dat_ben <- dat_ben[!((dat_ben[,"x"] == 0) & (dat_ben[,"delta"] == 1)),]

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
                    covnames = c("l1", "l2"))

# visualizing the results
library(dplyr)
library(ggplot2)
res <- plotBB(fit)
