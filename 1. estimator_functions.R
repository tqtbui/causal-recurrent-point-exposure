#################################
#Functions to produce estimator in Baer et al 2025
#################################

#-------------------
# Fit propensity score estimator 
#-------------------

pi_fit <- function(A, L, newL, 
                   pi.library = c("SL.glm", "SL.gam", "SL.lgb"), 
                   save.fit = FALSE, pi.trunc = 0, ...) {
  
  # output object
  ret <- list()
  
  # fit using Superlearner
  fit <- SuperLearner(
    Y = A, 
    X = L, 
    newX = newL, 
    family='binomial',
    SL.library = pi.library, ...
  )
  
  # output prediction and/or the fit if required
  ret$pi.hat <- c(fit$SL.predict)
  if(save.fit) ret$pi.fit <- fit
  
  # truncation
  ret$pi.hat <- pmax(ret$pi.hat, pi.trunc)
  ret$pi.hat <- pmin(ret$pi.hat, 1-pi.trunc)
  
  # return
  return(ret)
} 

#-------------------
# Fit survival estimator for terminal event and censoring k(t; a, l), h(t; a, l)
#-------------------

surv_fit <- function(X, Delta, A, L, newL, eval.times, 
                     event.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                     cens.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                     save.fit = FALSE, k.trunc = 0, h.trunc = 0, eps = 1e-06, ...) {
  
  # output object
  ret <- list(eval.times = eval.times)
  
  # prepare covariates
  AL <- cbind(A, L)
  newAL <- rbind(cbind(A = 0, newL), 
                 cbind(A = 1, newL))
  
  # fitting survival function to survival data
  fit <- survSuperLearner(
    time = X,
    event = Delta,
    X = AL,
    newX = newAL,
    new.times = eval.times,
    event.SL.library = event.library,
    cens.SL.library = cens.library, ...
  )
  
  # output fit if required
  if(save.fit) ret$surv.fit <- fit
  
  # return predictions for A = 0
  ret$h.hat.0 <- fit$event.SL.predict[1:nrow(newL),]
  ret$h.hat.0 <- pmax(ret$h.hat.0, h.trunc) # truncation
  if(any(ret$h.hat.0 == 0)) ret$h.hat.0[ret$h.hat.0 == 0] <- min(ret$h.hat.0[ret$h.hat.0 > 0]) # prevent 0
  
  ret$k.hat.0 <- fit$cens.SL.predict[1:nrow(newL),]
  ret$k.hat.0 <- pmax(ret$k.hat.0, k.trunc) # truncation
  if(any(ret$k.hat.0 == 0)) ret$k.hat.0[ret$k.hat.0 == 0] <- min(ret$k.hat.0[ret$k.hat.0 > 0]) # prevent 0
  
  # return predictions for A = 1
  ret$h.hat.1 <- fit$event.SL.predict[-(1:nrow(newL)),]
  ret$h.hat.1 <- pmax(ret$h.hat.1, h.trunc) # truncation
  if(any(ret$h.hat.1 == 0)) ret$h.hat.1[ret$h.hat.1 == 0] <- min(ret$h.hat.1[ret$h.hat.1 > 0]) # prevent 0
  
  ret$k.hat.1 <- fit$cens.SL.predict[-(1:nrow(newL)),]
  ret$k.hat.1 <- pmax(ret$k.hat.1, k.trunc) # truncation
  if(any(ret$k.hat.1 == 0)) ret$k.hat.1[ret$k.hat.1 == 0] <- min(ret$k.hat.1[ret$k.hat.1 > 0]) # prevent 0
  
  # return fitted values for K(X)
  pred <- predict.survSuperLearner(fit,
                                   newdata = AL,
                                   new.times = eval.times)
  k.fitted <- pred$cens.SL.predict
  k.fitted <- pmax(k.fitted, k.trunc) # truncation
  if(any(k.fitted == 0)) k.fitted[k.fitted == 0] <- min(k.fitted[k.fitted > 0]) # prevent 0
  ret$k.x.m <- sapply(1:nrow(AL), function(i) stepfun(eval.times, c(1,k.fitted[i,]), right = FALSE)(X[i]-eps))
  
  # return
  return(ret)
}

#-------------------
# Fit b function E[Delta/K * I(X>=u) * NE(t) | X>=u, A, L] 
#-------------------

b_fit <- function(X, Delta, A, L, newL, k.x.m, NE, 
                  eval.times, fit.times = quantile(eval.times, probs = seq(0, 0.95, by = 0.05)),
                  c.library = c("SL.glm", "SL.gam", "SL.lgb"), 
                  d.library = c("SL.glm", "SL.gam", "SL.lgb"),
                  eps = 1e-06, ...) {
  
  # output object
  nt <- ncol(NE)
  empt <- array(0, dim = c(nrow(newL), length(eval.times), nt))
  empt2 <- matrix(0, nrow = nrow(newL), ncol = length(fit.times))
  ret <- list(eval.times = eval.times, 
              b.hat.1 = empt, b.hat.0 = empt)
  
  # datasets
  idx.d <- Delta == 1
  dat.now <- cbind("A" = A, L)
  newAL <- rbind(cbind("A" = 1, newL), 
                 cbind("A" = 0, newL))
  
  # for each t fit a new b function
  for (t in 1:nt) {
    
    y <- Delta * NE[,t] / k.x.m
    b.tmp.1 <- empt2
    b.tmp.0 <- empt2
    remove.u <- c()
    
    c.pred.1 <- rep(0, nrow(newL))
    c.pred.0 <- rep(0, nrow(newL))
    d.pred.1 <- rep(0, nrow(newL))
    d.pred.0 <- rep(0, nrow(newL))
    
    for (u in 1:length(fit.times)) {
      
      idx.c <- X > fit.times[u]
      
      crit.c <- (length(unique(y[idx.c & idx.d])) == 1) | (sum(idx.c & idx.d) < 10)
      crit.d <- (length(unique(Delta[idx.c])) == 1) | (sum(idx.c) < 10)
      
      if (!crit.c) {
        
        # fit the c model
        fit.c <- SuperLearner(
          Y = y[idx.c & idx.d],          # I(X>u) = 1 and Delta = 1
          X = dat.now[idx.c & idx.d,],
          family = gaussian(),
          newX = newAL,
          method = "method.NNLS",
          SL.library = c.library, 
          ...
        )
        c.pred.1 <- fit.c$SL.predict[1:nrow(newL)] 
        c.pred.0 <- fit.c$SL.predict[-c(1:nrow(newL))]
      }
      
      if (!crit.d) {
        # fit the d model
        fit.d <- SuperLearner(
          Y = Delta[idx.c],              # I(X>u) = 1
          X = dat.now[idx.c,],
          family = binomial(),
          newX = newAL,
          method = "method.CC_LS",
          SL.library = d.library, 
          ...
        )
        d.pred.1 <- fit.d$SL.predict[1:nrow(newL)]
        d.pred.0 <- fit.d$SL.predict[-c(1:nrow(newL))]
      }
      
      b.tmp.1[,u] <- c.pred.1 * d.pred.1
      b.tmp.0[,u] <- c.pred.0 * d.pred.0
    }
    
    ret$b.hat.1[,,t] <- t(sapply(1:nrow(newL), function(i) {
      approx(x = fit.times, y = b.tmp.1[i,], xout = eval.times, rule = 2)$y
    }))
    ret$b.hat.0[,,t] <- t(sapply(1:nrow(newL), function(i) {
      approx(x = fit.times, y = b.tmp.0[i,], xout = eval.times, rule = 2)$y
    }))
    
  }
  
  # return
  return(ret)
}

#-------------------
# Estimate pi, k, h and b
#-------------------

estimate <- function(dat, t_fits = 2, 
                     kfolds = 5, eps = 1e-06, tau = 12, incr = 0.01,
                     pi.library = c("SL.glm", "SL.gam", "SL.lgb"), 
                     event.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                     cens.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                     c.library = c("SL.glm", "SL.gam", "SL.lgb"), 
                     d.library = c("SL.glm", "SL.gam", "SL.lgb"),
                     covnames = c("l1", "l2", "l3"),
                     verbose = TRUE) {
  
  # extract information from the data
  delta <- dat[,"delta"] # death indicator
  x <- dat[,"x"]         # death or censoring times
  a <- dat[,"a"]         # treatment allocation
  l <- dat[,covnames, drop = FALSE] # covariates
  nsamples <- nrow(dat)
  
  # number of events
  n_e <- dat[,paste0("NE_",t_fits), drop = FALSE] 
  
  # estimation grid
  u_grid_prime <- seq(0, tau, by = incr)
  u_grid <- sort(unique(c(u_grid_prime, x)))
  
  # result holder
  empt <- matrix(0, nrow = nsamples, ncol = length(u_grid))
  empt2 <- array(0, dim = c(nsamples, length(u_grid), length(t_fits)))
  ret <- list(delta = delta, x = x, a = a, l = l, 
              t_fits = t_fits, n_e = n_e,
              eval.times = u_grid,
              pi.hat = rep(0, nsamples), 
              k.hat.1 = empt, k.hat.0 = empt, 
              h.hat.1 = empt, h.hat.0 = empt, 
              b.hat.1 = empt2, b.hat.0 = empt2)
  
  # creating folds so that delta 1 and 0 appear in all folds
  ae.00 <- which((delta == 0) & (a == 0))
  ae.01 <- which((delta == 0) & (a == 1))
  ae.10 <- which((delta == 1) & (a == 0))
  ae.11 <- which((delta == 1) & (a == 1))
  
  folds.00 <- sample(rep(1:kfolds, length = length(ae.00)))
  folds.01 <- sample(rep(1:kfolds, length = length(ae.01)))
  folds.10 <- sample(rep(1:kfolds, length = length(ae.10)))
  folds.11 <- sample(rep(1:kfolds, length = length(ae.11)))
  
  idx <- rep(NA, nsamples)
  idx[ae.00] <- folds.00
  idx[ae.01] <- folds.01
  idx[ae.10] <- folds.10
  idx[ae.11] <- folds.11
  
  if (verbose) {
    cat("Start fitting ... \n")
  }
  
  # cross fitting
  for (fold in 1:kfolds) {
    
    # training and test indexes
    idx_tst <- idx == fold
    idx_trn <- idx != fold
    
    # training data
    delta_trn <- delta[idx_trn] # death indicator
    x_trn <- x[idx_trn]         # survival time
    a_trn <- a[idx_trn]         # treatment indicator
    l_trn <- l[idx_trn,,drop = FALSE]        # covariate
    n_e_trn <- n_e[idx_trn,,drop=FALSE]    # number of events (training set)
    
    delta_tst <- delta[idx_tst] # death indicator
    x_tst <- x[idx_tst]         # survival time
    a_tst <- a[idx_tst]         # treatment indicator
    l_tst <- l[idx_tst,,drop = FALSE]        # covariate
    n_e_tst <- n_e[idx_tst,,drop=FALSE]    # number of events (test set)
      
    # Fitting the propensity score model
    fit.pi <- pi_fit(A = a_trn, 
                     L = as.data.frame(l_trn), 
                     newL = as.data.frame(l_tst), 
                     pi.library = pi.library)
    
    ret$pi.hat[idx_tst] <- fit.pi$pi.hat 
    
    # Fit survival estimator for terminal event and censoring k(t; a, l), h(t; a, l)
    fit.surv <- surv_fit(X = x_trn, Delta = delta_trn, A = a_trn, 
                         L = as.data.frame(l_trn), 
                         newL = as.data.frame(l_tst), 
                         eval.times = u_grid, eps = eps, 
                         event.library = event.library, 
                         cens.library = cens.library)
    
    ret$k.hat.1[idx_tst,] <- fit.surv$k.hat.1
    ret$k.hat.0[idx_tst,] <- fit.surv$k.hat.0
    ret$h.hat.1[idx_tst,] <- fit.surv$h.hat.1
    ret$h.hat.0[idx_tst,] <- fit.surv$h.hat.0
    
    k.x.m <- fit.surv$k.x.m
  
    # Fit b function
    fit.b <- b_fit(X = x_trn, Delta = delta_trn, A = a_trn, 
                   L = as.data.frame(l_trn), 
                   newL = as.data.frame(l_tst), 
                   k.x.m = k.x.m, NE = n_e_trn, 
                   eval.times = u_grid,
                   fit.times = quantile(u_grid_prime, probs = seq(0, 0.95, by = 0.05)),
                   c.library = c.library, 
                   d.library = d.library,
                   eps = eps) 
    
    ret$b.hat.1[idx_tst,,] <- fit.b$b.hat.1
    ret$b.hat.0[idx_tst,,] <- fit.b$b.hat.0
    
    if (verbose) {
      cat("Fold ", fold, " finished ... \n")
    }
    
  } # end of cross fitting
  
  # return
  return(ret)
}

#-------------------
# Summarize one simulation into estimators
#-------------------

# function to obtain nuisance vectors for estimators
inputvecs <- function(obj, 
                      cutoff = 0, eps = 1e-06) {
  
  # get some basic info
  delta <- obj$delta
  x <- obj$x
  a <- obj$a
  l <- obj$l
  t_fits <- obj$t_fits
  n_e <- obj$n_e
  eval.times <- obj$eval.times
  
  # dimensions
  nsamples <- length(delta)
  nt <- length(t_fits)
  nu <- length(eval.times)
  i_x_t <- sapply(1:nt, function(j){
    x > t_fits[j]
  })
  
  # truncate the probabilities
  pi.hat <- pmax(obj$pi.hat, cutoff)
  pi.hat <- pmin(pi.hat, 1-cutoff)
  k.hat.1 <- pmax(obj$k.hat.1, cutoff) 
  k.hat.0 <- pmax(obj$k.hat.0, cutoff) 
  h.hat.1 <- pmax(obj$h.hat.1, cutoff) 
  h.hat.0 <- pmax(obj$h.hat.0, cutoff) 
  
  # calculate the necessary quantities
    
  k.x.1 <- sapply(1:nsamples, function(i) stepfun(eval.times, c(1,k.hat.1[i,]), right = FALSE)(x[i]))
  k.x.0 <- sapply(1:nsamples, function(i) stepfun(eval.times, c(1,k.hat.0[i,]), right = FALSE)(x[i]))
  k.x <- ifelse(a == 1, k.x.1, k.x.0)
  
  k.x.m.1 <- sapply(1:nsamples, function(i) stepfun(eval.times, c(1,k.hat.1[i,]), right = FALSE)(x[i]-eps))
  k.x.m.0 <- sapply(1:nsamples, function(i) stepfun(eval.times, c(1,k.hat.0[i,]), right = FALSE)(x[i]-eps))
  k.x.minus <- ifelse(a == 1, k.x.m.1, k.x.m.0)
  
  h.x.1 <- sapply(1:nsamples, function(i) stepfun(eval.times, c(1,h.hat.1[i,]), right = FALSE)(x[i]))
  h.x.0 <- sapply(1:nsamples, function(i) stepfun(eval.times, c(1,h.hat.0[i,]), right = FALSE)(x[i]))
  h.x <- ifelse(a == 1, h.x.1, h.x.0)
  
  # h(t)
  h.t.1 <- sapply(1:nsamples, function(i) stepfun(eval.times, c(1,h.hat.1[i,]), right = FALSE)(t_fits))
  h.t.0 <- sapply(1:nsamples, function(i) stepfun(eval.times, c(1,h.hat.0[i,]), right = FALSE)(t_fits))
  
  # h(x v t)
  h.x.v.t.1 <- sapply(1:nsamples, function(i) stepfun(eval.times, c(1,h.hat.1[i,]), right = FALSE)(sapply(1:nt, function(j) {max(t_fits[j], x[i])})))
  h.x.v.t.0 <- sapply(1:nsamples, function(i) stepfun(eval.times, c(1,h.hat.0[i,]), right = FALSE)(sapply(1:nt, function(j) {max(t_fits[j], x[i])})))
  
  if (nt == 1) {
    h.x.v.t <- sapply(1:nsamples, function(i) {
      ifelse(a[i] == 1, h.x.v.t.1[i], h.x.v.t.0[i]) 
    })
  } else {
    h.x.v.t <- sapply(1:nsamples, function(i) {
      if (a[i] == 1) {
        tmp <- h.x.v.t.1[,i]
      } else {
        tmp <- h.x.v.t.0[,i]
      }
      return(tmp)
    })
  }
  
  # reformat
  pi.hat <- ifelse(a == 1, pi.hat, 1-pi.hat) #pi(A)
  if (nt == 1) {
    h.t.1 <- as.matrix(h.t.1)
    h.t.0 <- as.matrix(h.t.0)
    h.x.v.t <- as.matrix(h.x.v.t)
  } else {
    h.t.1 <- t(h.t.1)
    h.t.0 <- t(h.t.0)
    h.x.v.t <- t(h.x.v.t)
  }
  
  # quantities related to F
  b.hat.1 <- obj$b.hat.1
  b.hat.0 <- obj$b.hat.0
  
  f.0.1 <- sapply(1:nt, function(j) {b.hat.1[,1,j]}) #because HK = 1 when u = 0
  f.0.0 <- sapply(1:nt, function(j) {b.hat.0[,1,j]})
  
  b.x.1 <- sapply(1:nt, function(j) {
    sapply(1:nsamples, function(i) stepfun(eval.times, c(b.hat.1[i,1,j],b.hat.1[i,,j]), right = FALSE)(x[i]))
  })
  b.x.0 <- sapply(1:nt, function(j) {
    sapply(1:nsamples, function(i) stepfun(eval.times, c(b.hat.0[i,1,j],b.hat.0[i,,j]), right = FALSE)(x[i]))
  })
  b.x <- sapply(1:nt, function(j) {
    ifelse(a==1, b.x.1[,j], b.x.0[,j])
  })
  
  # integrals for f
  int.f.1 <- sapply(1:nt, function(j) {
    sapply(1:nsamples, function(i) {
      vals <- diff(-log(k.hat.1[i,])) * ((b.hat.1[i,-nu,j] + b.hat.1[i,-1,j])/2) #trapezoidal rule
      if(any(eval.times[-1] > x[i])) vals[eval.times[-1] > x[i]] <- 0 #only integrate until x
      sum(vals)
    })
  })
  int.f.0 <- sapply(1:nt, function(j) {
    sapply(1:nsamples, function(i) {
      vals <- diff(-log(k.hat.0[i,])) * ((b.hat.0[i,-nu,j] + b.hat.0[i,-1,j])/2) #trapezoidal rule
      if(any(eval.times[-1] > x[i])) vals[eval.times[-1] > x[i]] <- 0 #only integrate until x
      sum(vals)
    })
  })
  int.f <-  sapply(1:nt, function(j) {
    ifelse(a == 1, int.f.1[,j], int.f.0[,j])
  })
  
  # integral for h
  int.h.1 <- sapply(1:nt, function(j) {
    uvt <- sapply(1:nu, function(k) {max(eval.times[k], t_fits[j])})
    sapply(1:nsamples, function(i) {
      huvt <- stepfun(eval.times, c(1,h.hat.1[i,]), right = FALSE)(uvt)
      vals <- diff(1/k.hat.1[i,]) * ((huvt[-nu] / h.hat.1[i,-nu] + huvt[-1] / h.hat.1[i,-1])/2)
      if(any(eval.times[-1] > x[i])) vals[eval.times[-1] > x[i]] <- 0
      sum(vals)
    })
  })
  int.h.0 <- sapply(1:nt, function(j) {
    uvt <- sapply(1:nu, function(k) {max(eval.times[k], t_fits[j])})
    sapply(1:nsamples, function(i) {
      huvt <- stepfun(eval.times, c(1,h.hat.0[i,]), right = FALSE)(uvt)
      vals <- diff(1/k.hat.0[i,]) * ((huvt[-nu] / h.hat.0[i,-nu] + huvt[-1] / h.hat.0[i,-1])/2)
      if(any(eval.times[-1] > x[i])) vals[eval.times[-1] > x[i]] <- 0
      sum(vals)
    })
  })
  int.h <-  sapply(1:nt, function(j) {
    ifelse(a == 1, int.h.1[,j], int.h.0[,j])
  })
  
  return(list(a = a, delta = delta, n_e = n_e, i_x_t = i_x_t, 
              pi.hat = pi.hat, k.x.minus = k.x.minus, k.x = k.x, 
              h.x = h.x, h.x.v.t = h.x.v.t, h.t.1 = h.t.1, h.t.0 = h.t.0, 
              f.0.1 = f.0.1, f.0.0 = f.0.0, b.x = b.x, 
              int.f = int.f, int.h = int.h, 
              cutoff = cutoff, 
              t_fits = t_fits, nsamples = nsamples))
}

# function to calculate the estimators and put them into a dataframe
outputs <- function(obj) {
  
  # output object
  out <- data.frame(size = double(),
                    method = character(), 
                    cutoff = character(), 
                    time = double(),
                    estimand = character(), 
                    est = double(), 
                    sd = double())
  
  # IPW weights
  ipwts_1 <- (as.numeric(obj$a == 1) * obj$delta) / (obj$pi.hat * obj$k.x.minus) 
  ipwts_0 <- (as.numeric(obj$a == 0) * obj$delta) / (obj$pi.hat * obj$k.x.minus)
  
  cnt <- 0
  for (j in 1:length(obj$t_fits)) {
    
    # IPW part
    ipwest <- cbind(ipwts_1 * obj$n_e[,j], 
                    ipwts_0 * obj$n_e[,j], 
                    ipwts_1 * obj$i_x_t[,j], 
                    ipwts_0 * obj$i_x_t[,j])
    
    aug_1 <- cbind(((as.numeric(obj$a == 1) - obj$pi.hat) * obj$f.0.1[,j]) / obj$pi.hat, 
                   ((as.numeric(obj$a == 0) - obj$pi.hat) * obj$f.0.0[,j]) / obj$pi.hat,
                   ((as.numeric(obj$a == 1) - obj$pi.hat) * obj$h.t.1[,j]) / obj$pi.hat,
                   ((as.numeric(obj$a == 0) - obj$pi.hat) * obj$h.t.0[,j]) / obj$pi.hat)
    
    aug_2 <- cbind((as.numeric(obj$a == 1) * (1-obj$delta) * obj$b.x[,j]) / obj$pi.hat,
                   (as.numeric(obj$a == 0) * (1-obj$delta) * obj$b.x[,j]) / obj$pi.hat,
                   (as.numeric(obj$a == 1) * (1-obj$delta) * obj$h.x.v.t[,j]) / (obj$pi.hat * obj$h.x * obj$k.x),
                   (as.numeric(obj$a == 0) * (1-obj$delta) * obj$h.x.v.t[,j]) / (obj$pi.hat * obj$h.x * obj$k.x))
    
    aug_3 <- cbind((as.numeric(obj$a == 1) * obj$int.f[,j]) / obj$pi.hat, 
                   (as.numeric(obj$a == 0) * obj$int.f[,j]) / obj$pi.hat, 
                   (as.numeric(obj$a == 1) * obj$int.h[,j]) / obj$pi.hat, 
                   (as.numeric(obj$a == 0) * obj$int.h[,j]) / obj$pi.hat)
    
    infl <- ipwest - aug_1 + aug_2 - aug_3
    
    #horvitz-thompson - AIPW
    est <- colSums(infl) / obj$nsamples
    vars <- cbind(infl[,1] - est[1], 
                  infl[,2] - est[2], 
                  infl[,3] - est[3], 
                  infl[,4] - est[4])
    EU <- diag(c(-obj$nsamples, -obj$nsamples, 
                 -obj$nsamples, -obj$nsamples))
    vars <- solve(EU) %*% t(vars) %*% vars %*% t(solve(EU)) 
    out[cnt+1:4,] <- list(size = obj$nsamples,
                          method = "One-step AIPW", 
                          cutoff = obj$cutoff, 
                          time = obj$t_fits[j],
                          estimand = c("mu_1", "mu_0", "eta_1", "eta_0"), 
                          est = est, 
                          sd = sqrt(diag(vars)))
    
    cnt <- cnt+4
  }
  
  return(out)
}

#-------------------
# Baer estimation function
#-------------------

BBestimators <- function(dat, t_fits = 2, kfolds = 5, tau = 12, 
                         verbose = TRUE, cutoff = 0, eps = 1e-06, 
                         pi.library = c("SL.glm", "SL.gam", "SL.lgb"), 
                         event.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                         cens.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                         c.library = c("SL.glm", "SL.gam", "SL.lgb"), 
                         d.library = c("SL.glm", "SL.gam", "SL.lgb"), 
                         covnames = c("l1", "l2", "l3")) {
  
  nuisance <- estimate(dat = dat, t_fits = t_fits, kfolds = kfolds, tau = tau,
                       verbose = verbose, eps = eps, 
                       pi.library = pi.library, 
                       event.library = event.library, 
                       cens.library = cens.library, 
                       c.library = c.library, 
                       d.library = d.library, 
                       covnames = covnames)
  inputs <- inputvecs(nuisance, cutoff = cutoff, eps = eps)
  out <- outputs(inputs)
  
  return(out)
}

#plot the BB estimator
plotBB <- function(obj, isotonize = TRUE, 
                   transf = TRUE, eps = 1e-06) {
  
  #adding time=0
  obj <- rbind(obj[1:4,], obj)
  obj
  obj[1:4,"time"] <- 0
  obj[1:4,"est"] <- c(0,0,1,1)
  obj[1:4,"sd"] <- 0
  obj[1:4,"estimand"] <- c("mu_1", "mu_0", "eta_1", "eta_0")
  
  #isotonic regression
  if (isotonize) {
    #only need isotonize the ub and lb if we do uniform confidence bands, but here we have pointwise CIs
    obj$est[obj$estimand == "mu_1"] <- isoreg(obj$time[obj$estimand == "mu_1"], obj$est[obj$estimand == "mu_1"])$yf
    obj$est[obj$estimand == "mu_0"] <- isoreg(obj$time[obj$estimand == "mu_0"], obj$est[obj$estimand == "mu_0"])$yf
    
    obj$est[obj$estimand == "eta_1"] <- 1-isoreg(obj$time[obj$estimand == "eta_1"], 1-obj$est[obj$estimand == "eta_1"])$yf
    obj$est[obj$estimand == "eta_0"] <- 1-isoreg(obj$time[obj$estimand == "eta_0"], 1-obj$est[obj$estimand == "eta_0"])$yf
  }
  
  #transform confidence intervals
  if (transf) {
    obj.1 <- obj %>%
      filter(estimand %in% c("mu_1", "mu_0")) %>%
      mutate(estp = pmax(est, eps),
             g_mu = log(estp),
             g_mu_sd = sqrt(((1/estp) * sd)^2),
             lb = exp(g_mu - 1.96*g_mu_sd),
             ub = exp(g_mu + 1.96*g_mu_sd),
             gest = exp(g_mu))
    obj.2 <- obj %>%
      filter(estimand %in% c("eta_1", "eta_0")) %>%
      mutate(estp = pmax(est, eps),
             estp = pmin(estp, 1-eps),
             g_mu = log(-log(estp)),
             g_mu_sd = sqrt(((1/(estp*log(estp))) * sd)^2),
             lb = exp(-exp(g_mu + 1.96*g_mu_sd)),
             ub = exp(-exp(g_mu - 1.96*g_mu_sd)),
             gest = exp(-exp(g_mu)))
    obj <- rbind.data.frame(obj.1, obj.2)
  } else {
    obj <- obj %>%
      mutate(lb = est - 1.96*sd, 
             ub = est + 1.96*sd)
  }
  
  print(obj %>%
    ggplot(aes(x = time, y = est)) +
    geom_line() +
    geom_line(aes(x = time, y = lb), linetype = "dashed") +
    geom_line(aes(x = time, y = ub), linetype = "dashed") + 
    facet_wrap(vars(estimand), scales = "free_y"))
  
  return(obj)
}

#################################
#Functions to calculate comparative estimators
#################################

#-------------------
#Westling for eta
#-------------------

TWestimators <- function(dat, t_fits = 2, kfolds = 5, 
                         verbose = TRUE, cutoff = 0, tau = 12, 
                         pi.library = c("SL.glm", "SL.gam", "SL.lgb"), 
                         event.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                         cens.library = c("survSL.rfsrc", "survSL.coxph", "survSL.gam"), 
                         covnames = c("l1", "l2", "l3")) {
  
  # extract information from the data
  delta <- dat[,"delta"] # death indicator
  x <- dat[,"x"]         # death or censoring times
  a <- dat[,"a"]         # treatment allocation
  l <- dat[,covnames] # covariates
  nsamples <- nrow(dat)
  
  # number of events
  n_e <- dat[,paste0("NE_",t_fits)] 
  if (length(t_fits) == 1) {
    n_e <- as.matrix(n_e)
  } 
  
  # result holder
  u_grid <- seq(0, tau, by = 0.01)
  
  fit <- CFsurvival(time = x, 
                    event = delta, 
                    treat = a, 
                    confounders = l, 
                    fit.times = t_fits, 
                    fit.treat = c(0,1), 
                    verbose = verbose,
                    conf.band = FALSE,
                    contrasts = NULL,
                    nuisance.options = list(prop.SL.library = pi.library,
                                            event.SL.library = event.library, 
                                            cens.SL.library = cens.library, 
                                            cens.trunc = cutoff, 
                                            prop.trunc = cutoff, 
                                            V = kfolds, 
                                            eval.times = u_grid,
                                            survSL.cvControl = list(V = kfolds), 
                                            save.nuis.fits = TRUE))
  
  # output object
  out <- data.frame(size = nsamples,
                    method = "Westling et al. 2023", 
                    cutoff = cutoff, 
                    time = fit$surv.df$time[fit$surv.df$time %in% t_fits],
                    estimand = paste0("eta_", fit$surv.df$trt[fit$surv.df$time %in% t_fits]), 
                    est = fit$surv.df$surv[fit$surv.df$time %in% t_fits], 
                    sd = fit$surv.df$se[fit$surv.df$time %in% t_fits])
  
  return(out)
}

#-------------------
#Schaubel & Zhang for mu
#-------------------

DSestimators <- function(dat, t_fits = 2, 
                         cutoff = 0, tau = 12, 
                         pi.library = c("SL.glm"), 
                         covnames = c("l1", "l2", "l3")) {
  
  # extract information from the data
  delta <- dat[,"delta"] # death indicator
  x <- dat[,"x"]         # death or censoring times
  a <- dat[,"a"]         # treatment allocation
  l <- dat[,covnames] # covariates
  nsamples <- nrow(dat)
  
  # number of events
  n_e <- dat[,paste0("NE_",t_fits)] 
  if (length(t_fits) == 1) {
    n_e <- as.matrix(n_e)
  } 
  
  # u grid
  u_grid <- seq(0, tau, by = 0.01)
  
  #interpolate N_E
  ts <- colnames(dat)[contains("NE_", vars = colnames(dat))]
  ts <- as.numeric(substr(ts, 4, nchar(ts)))
  N_E <- dat[,paste0("NE_",ts)]
  N_E <- cbind(0, N_E)
  ts <- union(c(0), ts)
  dN_E <- matrix(0, nrow = nsamples, ncol = length(u_grid))
  for (t in 1:(length(ts)-1)) {
    for (i in 1:nsamples) {
      dN_E[i,(u_grid >= ts[t]) & (u_grid < ts[t+1])] <- (N_E[i,t+1] - N_E[i,t])/(sum((u_grid >= ts[t]) & (u_grid < ts[t+1])))
    }
  }
  
  #Calculate N_C
  dN_C <- c()
  y <- c()
  for (u in u_grid) {
    iu <- x >= u
    y <- cbind(y, iu) #rows are indv
    dN_C <- cbind(dN_C, ifelse(delta == 1, 0, as.numeric((x >= u) & (x <= u + 0.01)))) #column is u
  }
  
  # output object
  out <- data.frame(size = double(),
                    method = character(), 
                    cutoff = character(), 
                    time = double(),
                    estimand = character(), 
                    est = double(), 
                    sd = double())
  cnt <- 0
  
  # get the weight w
  fit.pi <- pi_fit(A = a,
                   L = as.data.frame(l),
                   newL = as.data.frame(l),
                   pi.library = pi.library, 
                   pi.trunc = cutoff)
  p <- fit.pi$pi.hat
  w.1 <- as.numeric(a == 1)/p #w in the paper
  w.0 <- as.numeric(a == 0)/(1-p)
  w <- cbind(w.0, w.1)

  #estimation for each treatment
  for (trt in 0:1) {
    
    #compute G weights
    pi.hat <- c()
    for (u in u_grid) {
      iu <- x >= u
      pi.hat <- c(pi.hat, mean(as.numeric(a == trt) * iu)) #indexed by u
    }
    Lambda_C <- sapply(1:nsamples, function(i) { #loop over individuals
      cumsum(as.numeric(a[i] == trt)*dN_C[i,]/pi.hat)
    }) # columns are individuals
    G <- exp(-rowMeans(Lambda_C, na.rm = TRUE)) #consists of u
    
    #estimate
    mu_hat <- sapply(1:nsamples, function(i){
      cumsum(w[i,trt+1] * dN_E[i,]/G)
    }) #columns are indv
    mu_hat <- rowMeans(mu_hat, na.rm = TRUE) #elements are u
    
    #dM_C
    dM_C <- matrix(0, nrow = nsamples, ncol = length(u_grid))
    for (i in 1:nsamples) {
      for (u in 1:length(u_grid)) {
        dM_C[i,u] <- as.numeric(a[i] == trt)*(dN_C[i,u] - y[i,u] * mean(as.numeric(a == trt)*dN_C[,u]/pi.hat[u], na.rm=TRUE))
      }
    }
    
    #I
    I <- 0
    for (i in 1:nsamples) {
      I <- I + as.matrix(c(1,l[i,])) %*% t(as.matrix(c(1,l[i,]))) * p[i] * (1-p[i])
    }
    I <- I/nsamples
    
    #U
    U <- sapply(1:nsamples, function(i) {
      (a[i] - p[i])*c(1,l[i,])
    }) #columns are individuals
    U <- t(U) #rows are individuals
    
    #g
    if (trt == 0) {
      g <- (-1)^trt * as.numeric(a == trt) * p/(1-p)
    } else {
      g <- (-1)^trt * as.numeric(a == trt) * (1-p)/p
    }

    for (t in 1:length(t_fits)) {
      
      #mu_hat_G
      mu_hat_t <- mu_hat[max(which(u_grid <= t_fits[t]))]

      #variance
      d <- sapply(1:nsamples, function(i) {
        sum(dN_E[i,u_grid <= t_fits[t]] / G[u_grid <= t_fits[t]]) * g[i] * c(1,l[i,])
      }) #columns are individuals, rows are covariates
      d <- rowMeans(d, na.rm = TRUE) #vector with length = # of covariates
      
      sigma <- c()
      for (i in 1:nsamples) {
        p1 <- t(as.matrix(d)) %*% solve(I) %*% as.matrix(U[i,])
        p2 <- sum((mu_hat_t - mu_hat[u_grid<=t_fits[t]])*as.numeric(a[i] == trt)*dM_C[i,u_grid<=t_fits[t]]/pi.hat[u_grid<=t_fits[t]])
        p3 <- sum(w[i,trt+1]*dN_E[i,u_grid<=t_fits[t]]/G[u_grid<=t_fits[t]])
        p3 <- p3 - mu_hat_t
        sigma <- c(sigma, p1+p2+p3)
      }
      sigma <- mean(sigma^2, na.rm = TRUE)
      
      out[cnt+1,] <- list(size = nsamples,
                          method = "Schaubel & Zhang 2010", 
                          cutoff = cutoff, 
                          time = t_fits[t],
                          estimand = paste0("mu_", trt), 
                          est = mu_hat_t, 
                          sd = sqrt(sigma/nsamples))
      cnt <- cnt+1
    }
  }

  return(out)
}

#-------------------
#Functions to calculate IPW estimators
#-------------------

IPWestimators <- function(dat, t_fits = 2,
                          cutoff = 0, tau = 12, incr = 0.01,
                          pi.library = c("SL.glm"),
                          event.library = c("survSL.coxph", "survSL.coxph"),
                          cens.library = c("survSL.coxph", "survSL.coxph"),
                          covnames = c("l1", "l2", "l3"), 
                          covnames.surv = c("l1", "l2", "l3"),
                          p.zCoef = c(-0.5, 0.1, 0.2, 0.3),
                          c.zCoef = c(-2, 0.3, 0, 0, 0)) {
  
  # extract information from the data
  delta <- dat[,"delta"] # death indicator
  x <- dat[,"x"]         # death or censoring times
  a <- dat[,"a"]         # treatment allocation
  l <- dat[,covnames] # covariates
  nsamples <- nrow(dat)
  
  # number of events
  n_e <- dat[,paste0("NE_",t_fits)]
  if (length(t_fits) == 1) {
    n_e <- as.matrix(n_e)
  }
  
  # u grid
  u_grid <- seq(0, tau, by = incr)
  nus <- length(u_grid)
  
  #Calculate N_C
  dN_C <- c()
  for (u in u_grid) {
    dN_C <- cbind(dN_C, ifelse(delta == 1, 0, as.numeric((x >= u) & (x < u + incr)))) #column is u
  }
  dN_C_bar <- colSums(dN_C)
  
  # output object
  out <- data.frame(size = double(),
                    method = character(),
                    cutoff = character(),
                    time = double(),
                    estimand = character(),
                    est = double(),
                    sd = double())
  cnt <- 0
  
  #--------- PROPENSITY WEIGHTS
  fit.pi <- pi_fit(A = a,
                   L = as.data.frame(l),
                   newL = as.data.frame(l),
                   pi.library = pi.library,
                   pi.trunc = cutoff, save.fit = TRUE)
  
  p <- fit.pi$pi.hat
  w <- cbind(as.numeric(a == 0)/(1-p), as.numeric(a == 1)/p)
  
  #I
  I <- 0
  for (i in 1:nsamples) {
    I <- I + as.matrix(c(1,l[i,])) %*% t(as.matrix(c(1,l[i,]))) * p[i] * (1-p[i])
  }
  I <- I/nsamples
  
  #U
  U <- sapply(1:nsamples, function(i) {
    (a[i] - p[i])*c(1,l[i,])
  }) #columns are individuals
  
  #--------- SURVIVAL WEIGHTS
  lsurv <- dat[,covnames.surv]
  fit.K <- surv_fit(X = x,
                    Delta = delta,
                    A = a,
                    L = as.data.frame(lsurv),
                    newL = as.data.frame(lsurv),
                    event.library = event.library,
                    cens.library = cens.library,
                    eval.times = sort(x), save.fit = TRUE,
                    k.trunc = cutoff, h.trunc = cutoff, eps = 0)
  
  #fitted
  k.x.m <- fit.K$k.x.m
  gamma.hat <- fit.K$surv.fit$cens.fitLibrary$survSL.coxph_All$object$coefficients
  
  #covariate for this model
  Lstar <- cbind(a,lsurv)
  expterm <- as.vector(exp(Lstar %*% as.matrix(gamma.hat)))
  
  #indicator matrix
  ymat <- sapply(1:nsamples, function(i) {
    as.numeric(x[i] >= u_grid)
  }) #columns are individuals, rows are u
  
  #Y_i(u)exp(gamma*L_i)
  const <- sapply(1:nus, function(u) {
    ymat[u,]*expterm #columns are u, rows are indv
  })
  
  #sum[Y_i(u)exp(gamma*L_i)]
  s0 <- colSums(const) / nsamples #u
  
  #sum[Y_i(u)exp(gamma*L_i)*L_i]/sum[Y_i(u)exp(gamma*L_i)]
  s1 <- mapply(function(u) {
    t(Lstar) %*% as.matrix(const[,u]) / nsamples
  }, u = 1:nus, SIMPLIFY = FALSE)
  
  #sum[Y_i(u)exp(gamma*L_i)*L_i]/sum[Y_i(u)exp(gamma*L_i)]
  s2 <- mapply(function(u) {
    t(Lstar) %*% diag(const[,u]) %*% Lstar / nsamples
  }, u = 1:nus, SIMPLIFY = FALSE)
  
  #J
  J <- Reduce(`+`, Map(function(x,y,z,v) {
    (x/z - as.matrix(y/z)%*%t(as.matrix(y/z)))*v
  }, s2, s1, s0, dN_C_bar))/nsamples
  
  #W
  newdat <- cbind.data.frame(u_grid, 0, 0)
  for (i in 1:length(covnames)) {
    newdat <- cbind.data.frame(newdat, 0)
  }
  colnames(newdat) <- c("time", "event", "A", covnames)
  Lambda0 <- predict(fit.K$surv.fit$cens.fitLibrary$survSL.coxph_All$object, 
                     newdata = newdat, type = "expected")
  lambda0 <- c(0, diff(Lambda0))
  dM_C <- sapply(1:nsamples, function(i) {
    dN_C[i,] - ymat[,i]*lambda0*expterm[i]
  }) #columns are individuals, rows are u
  
  #V
  h_element <- Map(function(x,y,z) {
    x*z/y
  }, s1, s0, lambda0)
  lambda0.x <- sapply(1:nsamples, function(i) {
    Reduce(`+`, lambda0[which(u_grid <= x[i])])
  }) #columns are individuals
  V <- sapply(1:nsamples, function(i) {
    h.x <- Reduce(`+`, h_element[which(u_grid <= x[i])])
    h.x - lambda0.x[i] * Lstar[i,]
  }) #columns are individuals
  
  #s1/s0
  s1s0 <- sapply(1:nus, function(u) {
    s1[[u]]/s0[u]
  }) #columns are u, rows covariates
  
  for (trt in 0:1) {
    
    #IPW weights
    ipwts <- w[,trt+1]*delta/k.x.m
    
    #g
    if (trt == 0) {
      g <- (-1)^trt * as.numeric(a == trt) * p/(1-p)
    } else {
      g <- (-1)^trt * as.numeric(a == trt) * (1-p)/p
    }
    
    for (t in 1:length(t_fits)) {
      
      # n_e[,t] <- 1
      
      #mu_hat
      mu_vec <- ipwts*n_e[,t]
      mu_hat <- mean(mu_vec, na.rm = TRUE)
      
      #eta_hat
      eta_vec <- ipwts*as.numeric(x > t_fits[t])
      eta_hat <- mean(eta_vec, na.rm = TRUE)
      
      #d function
      d_mu <- sapply(1:nsamples, function(i) {
        g[i] * c(1,l[i,]) * delta[i] * n_e[i,t] / k.x.m[i]
      }) #columns are individuals, rows are covariates
      d_mu <- rowMeans(d_mu, na.rm = TRUE) #vector with length = # of covariates
      
      d_eta <- sapply(1:nsamples, function(i) {
        g[i] * c(1,l[i,]) * delta[i] * as.numeric(x[i] > t_fits[t]) / k.x.m[i]
      }) #columns are individuals, rows are covariates
      d_eta <- rowMeans(d_eta, na.rm = TRUE) #vector with length = # of covariates
      
      #f1 function
      f1_mu <- sapply(1:nsamples, function(i) {
        w[i,trt+1] * delta[i] * n_e[i,t] * ymat[,i] * expterm[i] / (k.x.m[i] * s0)
      }) #columns are individuals, rows are u
      f1_mu <- rowMeans(f1_mu, na.rm = TRUE)
      
      f1_eta <- sapply(1:nsamples, function(i) {
        w[i,trt+1] * delta[i] * as.numeric(x[i] > t_fits[t]) * ymat[,i] * expterm[i] / (k.x.m[i] * s0)
      }) #columns are individuals, rows are u
      f1_eta <- rowMeans(f1_eta, na.rm = TRUE)
      
      #f2 function
      f2_mu <- sapply(1:nsamples, function(i) {
        w[i,trt+1] * delta[i] * n_e[i,t] * expterm[i] * V[,i] / k.x.m[i]
      }) #columns are individuals, rows are covariates
      f2_mu <- rowMeans(f2_mu, na.rm = TRUE)
      
      f2_eta <- sapply(1:nsamples, function(i) {
        w[i,trt+1] * delta[i] * as.numeric(x[i] > t_fits[t]) * expterm[i] * V[,i] / k.x.m[i]
      }) #columns are individuals, rows are covariates
      f2_eta <- rowMeans(f2_eta, na.rm = TRUE) #covariates
      
      #calculate variance
      sigma_mu <- c()
      sigma_eta <- c()
      for (i in 1:nsamples) {
        
        p3_mu <- mu_vec[i] - mu_hat
        p3_eta <- eta_vec[i] - eta_hat
        
        p1_mu <- as.vector(t(as.matrix(d_mu)) %*% solve(I) %*% as.matrix(U[,i]))
        p1_eta <- as.vector(t(as.matrix(d_eta)) %*% solve(I) %*% as.matrix(U[,i]))
        
        p2_mu <- sum(f1_mu*dM_C[,i]) + t(as.matrix(f2_mu)) %*% solve(J) %*% sweep(s1s0, MARGIN = 1, Lstar[i,], "-") %*% as.matrix(dM_C[,i])
        p2_eta <- sum(f1_eta*dM_C[,i]) + t(as.matrix(f2_eta)) %*% solve(J) %*% sweep(s1s0, MARGIN = 1, Lstar[i,], "-") %*% as.matrix(dM_C[,i])
        
        sigma_mu <- c(sigma_mu, p1_mu + p3_mu + p2_mu)
        sigma_eta <- c(sigma_eta, p1_eta + p3_eta + p2_eta)
      }
      sigma_mu <- mean(sigma_mu^2, na.rm = TRUE)
      sigma_eta <- mean(sigma_eta^2, na.rm = TRUE)
      
      out[cnt+1,] <- list(size = nsamples,
                          method = "IPW parametric",
                          cutoff = cutoff,
                          time = t_fits[t],
                          estimand = paste0("mu_", trt),
                          est = mu_hat,
                          sd = sqrt(sigma_mu/nsamples))
      
      out[cnt+2,] <- list(size = nsamples,
                          method = "IPW parametric",
                          cutoff = cutoff,
                          time = t_fits[t],
                          estimand = paste0("eta_", trt),
                          est = eta_hat,
                          sd = sqrt(sigma_eta/nsamples))
      
      cnt <- cnt+2
    }
  }
  
  return(out)
}
