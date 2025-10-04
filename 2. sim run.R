#################################
#Simulation for estimator
#################################

#Load source code files
source("0. general functions.R")
source("1. estimator_functions.R")

#---------------------
# Global variables
#---------------------

# check some system parameters
# cat(Sys.getenv('OMP_NUM_THREADS'))

# Pass array value, sample size, and model number
s <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) # simulation replicate (1-1000)
scenario <- as.numeric(Sys.getenv('SCENARIO')) 

# parameters
d.zCoef <- c(-2, -0.5, 0.1, 0.1, -0.5, -0.3, 0.1) 
e.zCoef <- c(1, -0.5, 0.1, 0.1, -0.5, -0.1, -0.5)

p.zCoef1 <- c(-0.5, 0.1, 0.1, 0.5)
p.zCoef2 <- c(-0.5, 0, 0, 0.5)

c.zCoef1 <- c(-2, 0.5, 0, 0, 0)
c.zCoef2 <- c(-3, 0, 0.1, 0.1, 0.5)

# get the data from the scenario
if (scenario == 1) {
  p.zCoef <- p.zCoef1
  c.zCoef <- c.zCoef1
} else if (scenario == 2) {
  p.zCoef <- p.zCoef1
  c.zCoef <- c.zCoef2
} else if (scenario == 3) {
  p.zCoef <- p.zCoef2
  c.zCoef <- c.zCoef2
}

#---------------------
# Generate data
#---------------------

nsamples <- 2500
t_fits <- seq(0.5, 6, by = 0.5)
tau <- 12

set.seed(s)
dat.now <- array(0, dim = c(nsamples, 8 + length(t_fits)),
                 dimnames = list(1:nsamples,
                                 c("c", "d", "x", "delta",
                                   "a", "l1", "l2", "l3", paste0("NE_", t_fits))))

for (i in 1:nsamples) {
  tmp <- gen_instance(end = tau, d.zCoef = d.zCoef, e.zCoef = e.zCoef, 
                      c.zCoef = c.zCoef, p.zCoef = p.zCoef)
  dat.now[i,] <- c(tmp$c, tmp$d, tmp$x, tmp$delta, 
                   tmp$a, tmp$l1, tmp$l2, tmp$l3,  
                   n_func(t = t_fits, data = tmp$e))
}

cat("\n Finished generating data \n")

# cov names for IPW estimator
if (scenario == 3) {
  ipwcov <- c("l1","l2")
} else {
  ipwcov <- c("l1","l2","l3")
}

#---------------------
# Use different methods to fit the data
#---------------------

# number of sample sizes
nsamples <- c(500, 1000, 1500, 2000, 2500)

# specify pointwise timepoint to fit estimator: each time point where N(t) is estimated
t_fits <- 1:6 

# output holder
out <- data.frame(size = double(),
                  method = character(), 
                  cutoff = character(), 
                  time = double(),
                  estimand = character(), 
                  est = double(), 
                  sd = double())

# start simulation
start <- Sys.time()
cat("\n Simulation starts \n")
for (i in 1:length(nsamples)) {
  
  # parametric only
  cat("\n Calculating IPW parametric estimators \n")
  tmp <- try(IPWestimators(dat.now[1:nsamples[i],], t_fits = t_fits, 
                           covnames = ipwcov), 
             silent = FALSE)
  if (!("try-error" %in% class(tmp))) {
    out <- rbind.data.frame(out, tmp)
  }
  
  # SuperLearner
  cat("\n Calculating AIPW superlearner estimators \n")
  tmp <- try(BBestimators(dat.now[1:nsamples[i],], t_fits = t_fits, 
                          kfolds = 5, 
                          pi.library = c("SL.glm", "SL.lgb", "SL.randomForest"),
                          event.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                          cens.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                          c.library = c("SL.glm", "SL.lgb", "SL.randomForest"), 
                          d.library = c("SL.glm", "SL.lgb", "SL.randomForest"),
                          verbose = TRUE), 
             silent = FALSE)
  if (!("try-error" %in% class(tmp))) {
    out <- rbind.data.frame(out, tmp)
  }
  
  # Schaubel 2010
  cat("\n Calculating Schaubel estimators \n")
  tmp <- try(DSestimators(dat.now[1:nsamples[i],], t_fits = t_fits, 
                          pi.library = c("SL.glm")), 
             silent = FALSE)
  if (!("try-error" %in% class(tmp))) {
    out <- rbind.data.frame(out, tmp)
  }
  
  # Westling 2023
  cat("\n Calculating Westling estimators \n")
  tmp <- try(TWestimators(dat.now[1:nsamples[i],], t_fits = t_fits, 
                          kfolds = 5, verbose = TRUE, 
                          pi.library = c("SL.glm", "SL.lgb", "SL.randomForest"), 
                          event.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                          cens.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam")), 
             silent = FALSE)
  if (!("try-error" %in% class(tmp))) {
    out <- rbind.data.frame(out, tmp)
  }
  
}

# Finish simulation
end <- Sys.time()
cat("Simulation finished after ", hms_span(start, end), "\n")

# save the results
save(out, dat.now,
     file = paste0("scenario ", scenario, "/sim_", s, ".RData"))




