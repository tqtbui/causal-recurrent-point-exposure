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

nuisance <- estimate(dat = dat_ben, t_fits = 1:60, kfolds = 5, tau = 64,
                     verbose = TRUE, 
                     pi.library = c("SL.glm"), 
                     event.library = c("survSL.km", "survSL.coxph"), 
                     cens.library = c("survSL.km", "survSL.coxph"), 
                     c.library = c("SL.glm"),
                     d.library = c("SL.glm"), 
                     covnames = c("l1", "l2"))

inputs <- inputvecs(nuisance)

#checking monotonicity
check_monotone <- function(vec) {
  ret <- TRUE
  tmp <- 0
  for (i in 1:(length(vec)-1)) {
    tmp.now <- ifelse(vec[i] == vec[i+1], 0, ifelse(vec[i] < vec[i+1], 1, -1))
    if (abs(tmp.now - tmp) == 2) {
      ret <- FALSE
      break
    } else {
      tmp <- tmp.now
      ret <- TRUE
    }
  }
  ret
}

check_k1 <- sapply(1:nrow(dat_ben), function(i) {
  check_monotone(nuisance$k.hat.1[i,])
})
sum(check_k1) != nrow(dat_ben)

check_k0 <- sapply(1:nrow(dat_ben), function(i) {
  check_monotone(nuisance$k.hat.0[i,])
})
sum(check_k0) != nrow(dat_ben)

check_h1 <- sapply(1:nrow(dat_ben), function(i) {
  check_monotone(nuisance$h.hat.1[i,])
})
sum(check_h1) != nrow(dat_ben)

check_h0 <- sapply(1:nrow(dat_ben), function(i) {
  check_monotone(nuisance$h.hat.0[i,])
})
sum(check_h0) != nrow(dat_ben)


check_b1 <- sapply(1:nrow(dat_ben), function(i) {
  sapply(1:dim(nuisance$b.hat.1)[2], function(j) {
    check_monotone(nuisance$b.hat.1[i,j,])
  })
})
check_b0 <- sapply(1:nrow(dat_ben), function(i) {
  sapply(1:dim(nuisance$b.hat.0)[2], function(j) {
    check_monotone(nuisance$b.hat.0[i,j,])
  })
})

res <- outputs(inputs)

tmp <- plotBB(res)

# res <- BBestimators(dat_ben, t_fits = 1:60, tau = 64, kfolds = 5, 
#                     verbose = TRUE, 
#                     pi.library = c("SL.glm"), 
#                     event.library = c("survSL.km", "survSL.coxph"), 
#                     cens.library = c("survSL.km", "survSL.coxph"), 
#                     c.library = c("SL.glm"),
#                     d.library = c("SL.glm"), 
#                     covnames = c("l1", "l2"))

#-----------------------------
#Visualization
#-----------------------------

# # visualizing the results
# library(dplyr)
# library(ggplot2)
# 
# #adding time=0
# obj <- res %>%
#   mutate(lb = est - 1.96*sd, 
#          ub = est + 1.96*sd)

# #transform confidence intervals
# eps <- 1e-09
# obj.1 <- obj %>%
#   filter(estimand %in% c("mu_1", "mu_0")) %>%
#   mutate(estp = pmax(est, eps),
#          g_mu = log(estp),
#          g_mu_sd = sqrt(((1/estp) * sd)^2),
#          lb = exp(g_mu - 1.96*g_mu_sd),
#          ub = exp(g_mu + 1.96*g_mu_sd), 
#          gest = exp(g_mu))
# obj.2 <- obj %>%
#   filter(estimand %in% c("eta_1", "eta_0")) %>%
#   mutate(estp = pmax(est, eps),
#          estp = pmin(estp, 1-eps),
#          g_mu = log(-log(estp)),
#          g_mu_sd = sqrt(((1/(estp*log(estp))) * sd)^2),
#          lb = exp(-exp(g_mu + 1.96*g_mu_sd)),
#          ub = exp(-exp(g_mu - 1.96*g_mu_sd)), 
#          gest = exp(-exp(g_mu)))
# obj <- rbind.data.frame(obj.1, obj.2)

# #adding time=0
# obj <- rbind(obj[1:4,], obj)
# obj
# obj[1:4,"time"] <- 0
# obj[1:4,"est"] <- c(0,0,1,1)
# obj[1:4,"sd"] <- 0
# obj[1:4,"lb"] <- c(0,0,1,1)
# obj[1:4,"ub"] <- c(0,0,1,1)
# 
# View(obj)
# 
# #plot
# png(filename = paste0("plot_bladder_result.png"), 
#     width = 900, height = 500)
# obj %>%
#   ggplot(aes(x = time, y = est)) +
#   geom_line() +
#   geom_line(aes(x = time, y = lb), linetype = "dashed") +
#   geom_line(aes(x = time, y = ub), linetype = "dashed") + 
#   facet_wrap(vars(estimand), scales = "free_y")
# dev.off()


#plot
png(filename = paste0("plot_bladder_result_untransformed.png"), 
    width = 900, height = 500)
tmp <- plotBB(fit, transf = FALSE)
dev.off()
