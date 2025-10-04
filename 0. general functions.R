#################################
#Loading libraries
#################################

library(reda) # generate recurrent event data
library(survSuperLearner) # Westling - ensemble learning to estimate conditional survival functions from right-censored time-to-event data
library(SuperLearner) # super learner
library(randomForestSRC) # fast unified random forests for survival, regression, and classification
library(survival) # survival analysis
library(CFsurvival) # Westling2023 paper
library(lightgbm) #light GBM
library(quadprog) #quadratic programming

# library(devtools)
# install.packages("reda")
# install.packages("SuperLearner")
# install.packages("randomForestSRC")
# install.packages("survival")
# install.packages("xgboost")
# install.packages("randomForest")
# install.packages("lightgbm")
# install.packages("quadprog")
# devtools::install_github("tedwestling/CFsurvival")
# devtools::install_github("tedwestling/survSuperLearner")

#################################
#Helper functions
#################################

#' Counting process function: extract event times that are less or equal than t
n_func <- function(t, data) {
  # data: data.frame output from simEventData function
  # t: scalar value for time
  times <- data[data[, "event"] == 1, "time"]
  
  out <- sapply(1:length(t), function(i) sum(times <= t[i]))
  return(out)
}

#-------------------

#' Time printing
hms_span <- function(start, end) {
  dsec <- as.numeric(difftime(end, start, unit = "secs"))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600*hours - 60*minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}

#-------------------

# calculate interactions between A and X
interact <- function(A, X) {
  sapply(1:ncol(X), function(i) A*X[,i])
}

#################################
#Functions to generate data
#################################

#' Weibull rate function: rho(t) = f(t)/S(t-) for baseline hazard
rho_weibull <- function(t, lambda, k) {
  # lambda: scale
  # k: shape
  out <- k * t ^ (k - 1) * lambda ^ (-k)
  return(out)
}

#-------------------

#' Sample generation function
gen_instance <- function(a = NULL, endCensoring = 10^6, end = 12, 
                         counterfactual = FALSE, 
                         p.zCoef = c(-0.5, 0.1, 0.2, 0.3),
                         c.zCoef = c(-2, 0.5, 0, 0, 0), 
                         d.zCoef = c(-3, -0.5, 0.3, -0.5, 0.1, -0.3, -0.1), 
                         e.zCoef = c(0, -0.5, 0.1, 0.3, -0.1, -0.3, 0.3)) {
  
  # Sample covariates l
  l1 <- rbinom(n = 1, size = 1, prob = 0.5)
  l2 <- runif(n = 1, min = -1, max = 1) 
  # l3 <- -0.5 + 3*rbeta(n = 1, shape1 = 2, shape2 = 2)
  l3 <- 0.5 + 3*rbeta(n = 1, shape1 = 2, shape2 = 2)
  
  # Sample treatment allocation a if not provided
  if (is.null(a)) {
    eta <- sum(c(1, l1, l2, l3)*p.zCoef)
    pi <- exp(eta) / (1 + exp(eta))
    a <- rbinom(n = 1, size = 1, prob = pi)
  }
  
  # Death simulation
  d <- simEventData(
    endTime = end, #admin censoring
    rho = rho_weibull, #baseline rate
    arguments = list(rho = list(lambda = 1, k = 1.1)), #coefs for baseline rate
    z = c(1, a, l1, l2, l3, a*l3, l1*l2), #covariates
    zCoef = d.zCoef, #coefs for death
    recurrent = FALSE
  )
  
  # Censoring simulation
  c <- simEventData(
    endTime = endCensoring, #censoring virtually not "censored"
    rho = rho_weibull, #baseline rate
    arguments = list(rho = list(lambda = 1, k = 1)), #coefs for baseline rate
    z = c(1, a, l1, l2, l3), #covariates
    zCoef = c.zCoef, #coefs for censoring 
    recurrent = FALSE
  ) 
  
  # Observed survival data: X and delta
  if (!counterfactual) { 
    
    # with censoring
    x <- min(c$time, d$time) #x
    delta <- ifelse(d$time <= c$time, 1, 0) #D' = D \wedge tau 
    
  } else { 
    
    # counterfactuals
    x <- d$time
    delta <- 1
    
  } 
  
  #' Recurrent event simulation
  e <- simEventData(
    endTime = x, #ends at death or censoring
    rho = rho_weibull, #baseline rate
    arguments = list(rho = list(lambda = 1, k = 1.1)), #coefs for baseline rate
    z = c(1, a, l1, l2, l3, a*l2, l1*l3), #covariates
    zCoef = e.zCoef #coefs for recurrent events
  )
  
  #return
  out <- list(
    "c" = c$time,
    "d" = d$time,
    "x" = x,
    "delta" = delta,
    "a" = a,
    "l1" = l1,
    "l2" = l2,
    "l3" = l3, 
    "e" = e[, -1] #remove the first column
  )
  
  return(out)
}

#################################
#True probability functions
#################################

#' True propensity score
true_pi <- function(l1, l2, l3, p.zCoef = c(-0.5, 0.1, 0.2, 0.3)) {
  
  eta <- sum(c(1, l1, l2, l3)*p.zCoef)
  out <- exp(eta) / (1 + exp(eta))
  
  return(out)
}

#-------------------

#' True survival function for death
true_h <- function(t, a, l1, l2, l3, 
                   d.zCoef = c(-3, -0.5, 0.3, -0.5, 0.1, -0.3, -0.1), 
                   lambda=1, k=1.1) {
  
  S0 <- exp(-(t/lambda)^k)
  S <- S0^(exp(sum(c(1, a, l1, l2, l3, a*l3, l1*l2)*d.zCoef)))
  
  return(S)
}

#-------------------

#' True survival function for censoring
true_k <- function(t, a, l1, l2, l3, 
                   c.zCoef =  c(-2, 0.5, 0, 0, 0), 
                   lambda=1, k=1) {
  
  S0 <- exp(-(t/lambda)^k)
  S <- S0^(exp(sum(c(1, a, l1, l2, l3)*c.zCoef)))
  
  return(S)
}

#################################
#SuperLearner wrappers
#################################

SL.lgb <- function(Y, X, newX, family, obsWeights, id,
                   depth = 3, rate = 0.01, rounds = 200) {
  stopifnot(require("lightgbm"))
  
  params <- list(task = "train",
                 max_depth = depth,
                 learning_rate = rate)
  
  if (family$family == 'gaussian') {
    params$objective <- "regression"
  } else if (family$family == 'binomial') {
    params$objective <- "binary"
  }
  
  m <- lgb.train(params,
                 lgb.Dataset(data = as.matrix(X), label = Y, weight = obsWeights),
                 nrounds = rounds, verbose = -1)
  
  pred <- predict(m, as.matrix(newX))
  
  fit = list(object = m)
  class(fit) <- 'SL.lgb'
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.lgb <- function(object, newdata, ...) {
  stopifnot(require("lightgbm"))
  pred <- predict(object$object, as.matrix(newdata))
  return(pred)
}

# Define HAL function for SuperLearner
SL.hal <- function(Y, X, newX, family, obsWeights, id,
                   num_knots = 10, max_degree = 1) {
  stopifnot(require("hal9001"))
  
  if (family$family == 'gaussian') {
    params$objective <- "regression"
  } else if (family$family == 'binomial') {
    params$objective <- "binary"
  }
  
  m <- fit_hal(as.matrix(X), Y,
               num_knots = num_knots,
               max_degree = max_degree)
  
  pred <- predict(m, as.matrix(newX))
  
  fit = list(object = m)
  class(fit) <- 'SL.hal'
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.SL.hal <- function(object, newdata, ...) {
  stopifnot(require("hal9001"))
  pred <- predict(object$object, as.matrix(newdata))
  return(pred)
}



