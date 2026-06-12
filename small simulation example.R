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
nsim <- 10 # number of simulations 

# output holder
out <- data.frame(size = double(),
                  method = character(), 
                  cutoff = character(), 
                  time = double(),
                  estimand = character(), 
                  est = double(), 
                  sd = double())

#---------------------
# Run simulations
#---------------------

# start simulation
start <- Sys.time()
cat("\n Simulation starts \n")
for (i in 1:nsim) {
  
  #generate data
  dat.now <- array(0, dim = c(nsamples, 8 + length(t_fits)),
                   dimnames = list(1:nsamples,
                                   c("c", "d", "x", "delta",
                                     "a", "l1", "l2", "l3", paste0("NE_", t_fits))))
  
  for (k in 1:nsamples) {
    one <- gen_instance(end = tau, d.zCoef = d.zCoef, e.zCoef = e.zCoef, 
                        c.zCoef = c.zCoef, p.zCoef = p.zCoef)
    dat.now[k,] <- c(one$c, one$d, one$x, one$delta, 
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
  tmp <- try(DSestimators(dat.now, t_fits = 1:6), 
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
  
  cat("Simulation:", i, "/", nsim, "\n")
  
}
# Finish simulation
end <- Sys.time()
cat("Simulation finished after ", hms_span(start, end), "\n")

# save the results
save(out, file = "simple_simulation.RData")

#---------------------
# Summarize simulation results
#---------------------

# take a look at the results of the simple simulation
load("simple_simulation.RData")
View(out)

#libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

#add counterfactual data
load("simulation_data.RData") #obtained from running the first part of file '3. visualization.R'
res.df <- out[-c(1:nrow(out)),,drop=FALSE]
for (i in 1:6) {
  tmp <- out %>%
    filter(time == i) %>%
    mutate(true = case_when(
      estimand == "mu_1" ~ mu_1[2*i], 
      estimand == "mu_0" ~ mu_0[2*i], 
      estimand == "eta_1" ~ eta_1[2*i], 
      estimand == "eta_0" ~ eta_0[2*i]
    ))
  res.df <- rbind.data.frame(res.df, tmp)
}

# select the columns
res.df <- res.df %>% 
  mutate(estimator = method) %>%
  select(size, time, estimand, true, sim, estimator, est, sd)

res.df$sd[res.df$est > 1 & res.df$estimand %in% c("eta_1", "eta_0") & res.df$estimator == "One-step AIPW"] <- NA
res.df$est[res.df$est > 1 & res.df$estimand %in% c("eta_1", "eta_0") & res.df$estimator == "One-step AIPW"] <- NA

# confidence intervals
eps <- 1e-06
res.df.1 <- res.df %>%
  filter(estimand %in% c("mu_1", "mu_0")) %>%
  mutate(estp = pmax(est, eps),
         g_mu = log(estp),
         g_mu_sd = sqrt(((1/estp) * sd)^2),
         lb = exp(g_mu - 1.96*g_mu_sd),
         ub = exp(g_mu + 1.96*g_mu_sd))

res.df.2 <- res.df %>%
  filter(estimand %in% c("eta_1", "eta_0")) %>%
  mutate(estp = pmax(est, eps),
         estp = pmin(estp, 1-eps),
         g_mu = log(-log(estp)),
         g_mu_sd = sqrt(((1/(estp*log(estp))) * sd)^2),
         lb = exp(-exp(g_mu + 1.96*g_mu_sd)),
         ub = exp(-exp(g_mu - 1.96*g_mu_sd)))

res.df <- rbind.data.frame(res.df.1, res.df.2) %>%
  select(size, time, estimand, sim,
         estimator, true, est, sd, lb, ub)

# bias and coverage
res.df <- res.df %>% 
  mutate(scaled_bias = sqrt(size)*(est-true), 
         relative_bias = (est-true)/true*100, 
         sq_dist = (est-true)^2,
         coverage = as.numeric((true >= lb) & (true <= ub)))

# sumarize across simulation numbers
res.df.summarize <- res.df %>%
  group_by(size, time, estimand, estimator) %>%
  summarise(scaled_bias__mean = mean(scaled_bias, na.rm = TRUE),
            relative_bias__mean = mean(relative_bias, na.rm = TRUE),
            theoretical_sd__mean = mean(sd, na.rm = TRUE), 
            empirical_sd__mean = sd(est, na.rm = TRUE), 
            coverage__mean = mean(coverage, na.rm = TRUE), 
            sd_ratio__mean = mean(sd/empirical_sd__mean, na.rm = TRUE),
            rmsq__mean = sqrt(mean(sq_dist, na.rm = TRUE)),
            #
            scaled_bias__sd = sd(scaled_bias, na.rm = TRUE), 
            relative_bias__sd = sd(relative_bias, na.rm = TRUE),
            theoretical_sd__sd = sd(sd, na.rm = TRUE), 
            sd_ratio__sd = sd(sd/empirical_sd__mean, na.rm = TRUE), 
            rmsq__sd = sqrt(sd(sq_dist, na.rm = TRUE))) %>%
  mutate(coverage__sd = sqrt(coverage__mean*(1-coverage__mean)/size), 
         empirical_sd__sd = theoretical_sd__sd)
names(res.df.summarize)

# average over time points
res.df.summarize.time <- res.df.summarize %>%
  group_by(size, estimand, estimator) %>%
  summarise(scaled_bias__mean = mean(scaled_bias__mean), 
            relative_bias__mean = mean(relative_bias__mean), 
            theoretical_sd__mean = mean(theoretical_sd__mean), 
            empirical_sd__mean = mean(empirical_sd__mean), 
            coverage__mean = mean(coverage__mean), 
            sd_ratio__mean = mean(sd_ratio__mean),
            rmsq__mean = mean(rmsq__mean),
            #
            scaled_bias__sd = mean(scaled_bias__sd), 
            relative_bias__sd = mean(relative_bias__sd), 
            theoretical_sd__sd = mean(theoretical_sd__sd), 
            empirical_sd__sd = mean(empirical_sd__sd), 
            coverage__sd = mean(coverage__sd), 
            sd_ratio__sd = mean(sd_ratio__sd), 
            rmsq__sd = mean(rmsq__sd))

res.df.summarize.long <- res.df.summarize.time %>%
  select(size, estimand, estimator, 
         scaled_bias__mean, relative_bias__mean, empirical_sd__mean, theoretical_sd__mean, coverage__mean, sd_ratio__mean, rmsq__mean,
         scaled_bias__sd, relative_bias__sd, empirical_sd__sd, theoretical_sd__sd, coverage__sd, sd_ratio__sd, rmsq__sd) %>%
  pivot_longer(cols = scaled_bias__mean:rmsq__sd, 
               names_to = c("summary", "mean_or_sd"), 
               names_sep = c("__")) %>%
  pivot_wider(names_from = mean_or_sd, 
              values_from = value)

#View the results
View(res.df.summarize.long)
