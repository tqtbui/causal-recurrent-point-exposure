#################################
#Simulation for estimator
#################################

#Load source code files
source("0. general functions.R")
source("1. estimator_functions.R")

#---------------------
# Global variables
#---------------------

# parameters
d.zCoef <- c(-2, -0.5, 0.1, 0.1, -0.5, -0.3, 0.1) 
e.zCoef <- c(1, -0.5, 0.1, 0.1, -0.5, -0.1, -0.5)
p.zCoef <- c(-0.5, 0.1, 0.1, 0.5)
c.zCoef <- c(-2, 0.5, 0, 0, 0)

nsamples <- 500 # sample size
t_fits <- seq(0.5, 6, by = 0.5) # time points
tau <- 12 # end of study
nsim <- 1 # number of simulations 

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
for (i in 1:nsim) {
  
  #generate data
  dat.now <- array(0, dim = c(nsamples, 8 + length(t_fits)),
                   dimnames = list(1:nsamples,
                                   c("c", "d", "x", "delta",
                                     "a", "l1", "l2", "l3", paste0("NE_", t_fits))))
  
  for (i in 1:nsamples) {
    one <- gen_instance(end = tau, d.zCoef = d.zCoef, e.zCoef = e.zCoef, 
                        c.zCoef = c.zCoef, p.zCoef = p.zCoef)
    dat.now[i,] <- c(one$c, one$d, one$x, one$delta, 
                     one$a, one$l1, one$l2, one$l3,  
                     n_func(t = t_fits, data = one$e))
  }
  
  # IPW estimator
  cat("\n Calculating IPW parametric estimators \n")
  tmp <- try(IPWestimators(dat.now, t_fits = 1:6), 
             silent = FALSE)
  if (!("try-error" %in% class(tmp))) {
    tmp$sim <- i
    out <- rbind.data.frame(out, tmp)
  }
  
  # AIPW with SuperLearner
  cat("\n Calculating AIPW superlearner estimators \n")
  tmp <- try(BBestimators(dat.now, t_fits = 1:6, 
                          kfolds = 5, verbose = TRUE,
                          pi.library = c("SL.glm"),
                          event.library = c("survSL.coxph", "survSL.km"), 
                          cens.library = c("survSL.coxph", "survSL.km"), 
                          c.library = c("SL.glm"), 
                          d.library = c("SL.glm")), 
             silent = FALSE)
  if (!("try-error" %in% class(tmp))) {
    tmp$sim <- i
    out <- rbind.data.frame(out, tmp)
  }
  
  # Schaubel 2010
  cat("\n Calculating Schaubel estimators \n")
  tmp <- try(DSestimators(dat.now, t_fits = 1:6, 
                          pi.library = c("SL.glm")), 
             silent = FALSE)
  if (!("try-error" %in% class(tmp))) {
    tmp$sim <- i
    out <- rbind.data.frame(out, tmp)
  }
  
  # Westling 2023
  cat("\n Calculating Westling estimators \n")
  tmp <- try(TWestimators(dat.now, t_fits = 1:6, 
                          kfolds = 5, verbose = TRUE, 
                          pi.library = c("SL.glm", "SL.lgb"), 
                          event.library = c("survSL.coxph", "survSL.km"), 
                          cens.library = c("survSL.coxph", "survSL.km")), 
             silent = FALSE)
  if (!("try-error" %in% class(tmp))) {
    tmp$sim <- i
    out <- rbind.data.frame(out, tmp)
  }
  
}
# Finish simulation
end <- Sys.time()
cat("Simulation finished after ", hms_span(start, end), "\n")

# save the results
save(out, file = "simple_simulation.RData")




