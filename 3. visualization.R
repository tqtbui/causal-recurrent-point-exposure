#################################
#Summarize results from the simulations into dataframes
#################################

result1 <- data.frame(size = double(),
                      simnum = integer(),
                      method = character(), 
                      cutoff = character(), 
                      time = double(),
                      estimand = character(), 
                      est = double(), 
                      sd = double())

result2 <- data.frame(size = double(),
                      simnum = integer(),
                      method = character(), 
                      cutoff = character(), 
                      time = double(),
                      estimand = character(), 
                      est = double(), 
                      sd = double())

result3 <- data.frame(size = double(),
                      simnum = integer(),
                      method = character(), 
                      cutoff = character(), 
                      time = double(),
                      estimand = character(), 
                      est = double(), 
                      sd = double())

for (i in 1:1000) {
  
  tmp <- try(load(paste0("scenario 1/sim_", i, ".RData")))
  if (!("try-error" %in% class(tmp))) {
    out$simnum <- i
    result1 <- rbind.data.frame(result1, out)
  }
  
  tmp <- try(load(paste0("scenario 2/sim_", i, ".RData")))
  if (!("try-error" %in% class(tmp))) {
    out$simnum <- i
    result2 <- rbind.data.frame(result2, out)
  }
  
  tmp <- try(load(paste0("scenario 3/sim_", i, ".RData")))
  if (!("try-error" %in% class(tmp))) {
    out$simnum <- i
    result3 <- rbind.data.frame(result3, out)
  }
}

result1$scenario <- 1
result2$scenario <- 2
result3$scenario <- 3
result <- rbind.data.frame(result1, result2, result3)

library(dplyr)
load("simulation_data.RData")
res.df <- result[-c(1:nrow(result)),,drop=FALSE]
for (i in 1:6) {
  tmp <- result %>%
    filter(time == i) %>%
    mutate(true = case_when(
      estimand == "mu_1" ~ mu_1[2*i], 
      estimand == "mu_0" ~ mu_0[2*i], 
      estimand == "eta_1" ~ eta_1[2*i], 
      estimand == "eta_0" ~ eta_0[2*i]
    ))
  res.df <- rbind.data.frame(res.df, tmp)
}

save(res.df, file = "results.RData")

#################################
#Plots
#################################

#calculate bias, variance ratio and coverage
load("results.RData")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# select the columns
res.df <- res.df %>% 
  mutate(estimator = method) %>%
  select(size, time, estimand, true, scenario, simnum, estimator, est, sd)

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
  select(size, time, estimand, scenario, simnum,
         estimator, true, est, sd, lb, ub)

# bias and coverage
res.df <- res.df %>% 
  mutate(scaled_bias = sqrt(size)*(est-true), 
         relative_bias = (est-true)/true*100, 
         sq_dist = (est-true)^2,
         coverage = as.numeric((true >= lb) & (true <= ub)))

# sumarize across simulation numbers
res.df.summarize <- res.df %>%
  group_by(size, time, estimand, scenario, estimator) %>%
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
  group_by(size, estimand, scenario, estimator) %>%
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

#------------------------------
# Reshaping datasets
#------------------------------

# long version of the df
res.df.summarize.long <- res.df.summarize.time %>%
  select(size, estimand, scenario, estimator, 
         scaled_bias__mean, relative_bias__mean, empirical_sd__mean, theoretical_sd__mean, coverage__mean, sd_ratio__mean, rmsq__mean,
         scaled_bias__sd, relative_bias__sd, empirical_sd__sd, theoretical_sd__sd, coverage__sd, sd_ratio__sd, rmsq__sd) %>%
  pivot_longer(cols = scaled_bias__mean:rmsq__sd, 
               names_to = c("summary", "mean_or_sd"), 
               names_sep = c("__")) %>%
  pivot_wider(names_from = mean_or_sd, 
              values_from = value)

# plots for a specific times
times <- 1:6
res.df.summarize.select <- res.df.summarize %>%
  filter(time %in% times) %>%
  select(size, estimand, scenario, estimator, 
         scaled_bias__mean, relative_bias__mean, empirical_sd__mean, theoretical_sd__mean, coverage__mean, sd_ratio__mean, rmsq__mean,
         scaled_bias__sd, relative_bias__sd, empirical_sd__sd, theoretical_sd__sd, coverage__sd, sd_ratio__sd, rmsq__sd) %>%
  pivot_longer(cols = scaled_bias__mean:rmsq__sd, 
               names_to = c("summary", "mean_or_sd"), 
               names_sep = c("__")) %>%
  pivot_wider(names_from = mean_or_sd, 
              values_from = value) %>%
  mutate(estimandp = case_when(
    estimand == "eta_1" ~ paste0("eta[1](", time, ")"),
    estimand == "eta_0" ~ paste0("eta[0](", time, ")"),
    estimand == "mu_1" ~ paste0("mu[1](", time, ")"),
    estimand == "mu_0" ~ paste0("mu[0](", time, ")")
  ))

#------------------------------
# Coverage plot
#------------------------------

# logit transformation
logit <- function(x) {1/(1+exp(-x))}

png(filename = paste0("plots/final/coverage_average.png"), 
    width = 1000, height = 700, res = 100)
res.df.summarize.long %>% 
  mutate(ref = case_when(
    summary == "coverage" ~ 0.95, 
    summary == "scaled_bias" ~ 0, 
    summary == "relative_bias" ~ 0, 
    summary == "empirical_sd" ~ 0,
    summary == "theoretical_sd" ~ 0,
    summary == "sd_ratio" ~ 1,
    summary == "rmsq" ~ 0
  )) %>%
  filter(summary %in% c("coverage")) %>%
  mutate(estimand = factor(estimand, 
                           levels = c("mu_1", "mu_0", "eta_1", "eta_0"), 
                           labels = c('mu[1](t)', 'mu[0](t)', 'eta[1](t)', 'eta[0](t)'),
                           ordered = TRUE), 
         scenario = factor(scenario, 
                           levels = c(1, 2, 3), 
                           labels = c("scenario~1", "scenario~2", "scenario~3"),
                           ordered = TRUE)) %>%
  mutate(ref = logit(ref), 
         ub = logit(mean+1.96*sd), 
         lb = logit(mean-1.96*sd), 
         mean = logit(mean)) %>%
  ggplot(aes(x = size, y = mean, color = estimator)) +
  geom_line() +
  geom_hline(aes(yintercept = ref), linetype = "solid", color = "darkgrey") +
  facet_grid(rows = vars(scenario), cols = vars(estimand), labeller = label_parsed) +
  scale_color_brewer(palette="Set1") +
  geom_errorbar(aes(ymin=lb, ymax=ub), width=200,
                position=position_dodge(0.05)) +
  labs(y = "coverage", x = "sample size") + 
  theme(legend.position="bottom")
dev.off()

png(filename = paste0("plots/final/coverage_points_mu.png"), 
    width = 1000, height = 700, res = 100)
res.df.summarize.select %>% 
  mutate(ref = case_when(
    summary == "coverage" ~ 0.95, 
    summary == "scaled_bias" ~ 0, 
    summary == "relative_bias" ~ 0, 
    summary == "empirical_sd" ~ 0,
    summary == "theoretical_sd" ~ 0,
    summary == "sd_ratio" ~ 1,
    summary == "rmsq" ~ 0
  )) %>%
  mutate(scenario = factor(scenario, 
                           levels = c(1, 2, 3), 
                           labels = c("scenario~1", "scenario~2", "scenario~3"),
                           ordered = TRUE)) %>%
  filter((summary %in% c("coverage")) & (estimand %in% c("mu_1", "mu_0"))) %>%
  ggplot(aes(x = size, y = mean, color = estimator)) +
  geom_line() +
  geom_hline(aes(yintercept = ref), linetype = "solid", color = "darkgrey") +
  facet_grid(rows = vars(scenario), cols = vars(estimandp), labeller = label_parsed) +
  scale_color_brewer(palette="Set1") +
  geom_errorbar(aes(ymin=mean-1.96*sd, ymax=mean+1.96*sd), width=200,
                position=position_dodge(0.05)) +
  labs(y = "coverage", x = "sample size") + 
  theme(legend.position="bottom") + 
  scale_x_continuous(breaks=c(500, 1000, 1500, 2000, 2500),
                     labels=c("", "1000", "", "2000", ""))
dev.off()

#------------------------------
# RMSE plot
#------------------------------

png(filename = paste0("plots/final/rmse_average.png"), 
    width = 1000, height = 700, res = 100)
res.df.summarize.long %>% 
  mutate(ref = case_when(
    summary == "coverage" ~ 0.95, 
    summary == "scaled_bias" ~ 0, 
    summary == "relative_bias" ~ 0, 
    summary == "empirical_sd" ~ 0,
    summary == "theoretical_sd" ~ 0,
    summary == "sd_ratio" ~ 1,
    summary == "rmsq" ~ 0
  )) %>%
  filter(summary %in% c("rmsq")) %>%
  mutate(estimand = factor(estimand, 
                           levels = c("mu_1", "mu_0", "eta_1", "eta_0"), 
                           labels = c('mu[1](t)', 'mu[0](t)', 'eta[1](t)', 'eta[0](t)'),
                           ordered = TRUE), 
         scenario = factor(scenario, 
                           levels = c(1, 2, 3), 
                           labels = c("scenario~1", "scenario~2", "scenario~3"),
                           ordered = TRUE)) %>%
  ggplot(aes(x = size, y = mean, color = estimator)) +
  geom_line() +
  geom_hline(aes(yintercept = ref), linetype = "solid", color = "darkgrey") +
  facet_grid(rows = vars(scenario), cols = vars(estimand), labeller = label_parsed) +
  scale_color_brewer(palette="Set1") +
  labs(y = "RMSE", x = "sample size") + 
  theme(legend.position="bottom")
dev.off()

png(filename = paste0("plots/final/rmse_points_eta.png"), 
    width = 1000, height = 700, res = 100)
res.df.summarize.select %>% 
  mutate(ref = case_when(
    summary == "coverage" ~ 0.95, 
    summary == "scaled_bias" ~ 0, 
    summary == "relative_bias" ~ 0, 
    summary == "empirical_sd" ~ 0,
    summary == "theoretical_sd" ~ 0,
    summary == "sd_ratio" ~ 1,
    summary == "rmsq" ~ 0
  )) %>%
  mutate(scenario = factor(scenario, 
                           levels = c(1, 2, 3), 
                           labels = c("scenario~1", "scenario~2", "scenario~3"),
                           ordered = TRUE)) %>%
  filter((summary %in% c("rmsq")) & (estimand %in% c("eta_1", "eta_0"))) %>%
  ggplot(aes(x = size, y = mean, color = estimator)) +
  geom_line() +
  geom_hline(aes(yintercept = ref), linetype = "solid", color = "darkgrey") +
  facet_grid(rows = vars(scenario), cols = vars(estimandp), labeller = label_parsed) +
  scale_color_brewer(palette="Set1") +
  labs(y = "RMSE", x = "sample size") + 
  theme(legend.position="bottom") + 
  scale_x_continuous(breaks=c(500, 1000, 1500, 2000, 2500),
                     labels=c("", "1000", "", "2000", ""))
dev.off()
