#individual project
#katherine hernandez

library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(reshape2)
library(arm)
library(rstan) #for extract function
library(rstanarm)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
theme_set(theme_classic())

#reading in file
weight = read.csv("Desktop/WeightTime.csv", as.is = TRUE)

# exploring the data structure
str(weight)
summary(weight)

# shaping the data set
weight = weight %>%
gather('four', 'six', 'eight', 'ten', 'twelve',key = "Time", value = "Weight")

#transforming weight gain on a log scale
weight$log_wt = log(weight$Weight)

# some exploratory plots
ggplot(weight, aes(Diet, weight$log_wt)) +
  geom_point(alpha = 0.5) + geom_smooth(method = "loess", span=99)

ggplot(weight, aes(Time, weight$log_wt)) +
  geom_point(alpha = 0.5) + geom_smooth(method = "loess", span=99)

weightplot = ggplot() + geom_boxplot(data = weight, position = position_dodge(width = 0.75), 
                                     aes(x = Diet, y = weight$log_wt, fill = Time), width = 0.5)
weightplot

# plotting abundance by diet
# there is some variation across diets
ggplot(weight, aes(Diet,weight$log_wt)) +
  geom_violin(alpha = 0.5) +
  geom_hline(yintercept = mean(weight$log_wt),
             color = "blue", linetype = "dashed") +
  coord_flip()


######################
#stan glmer model
# model looking at estimate of growth rates and how they vary across diets 
M1 = stan_glmer(log_wt ~ Time + Diet + Time:Diet + (1|Individual), 
           family=gaussian, weight)

summary(M1,probs=c(0.025,0.975),digits=2)


# assessment of model performance
# Rhats are at 1, so chains converged
# Neffs look relatively high (I think?) which is good
print(summary(M1)[,c(1,3,9,10)],digits=3)

# assessing the model performance with shinystan
# the histograms look good
# the trace plots converge
launch_shinystan(M1)

# model interpretation it looks like diet does have an effect on weight gain across individuals, with caterpillars
# that fed on lupine seeing an overall slower weight gain over time
plot(M1)

# not working for me, doesn't view log_wt as a parameter from my model
# not sure why
posterior_interval(M1, prob = 0.8, pars = "log_wt")

# Extract posterior samples
# tried to extract, I though that I was supposed to 
# extract stanfit from the model?
samples = extract(M1$stanfit)

names(samples)
str(samples$b) 
str(samples$beta) 

#some visualizations
library(bayesplot)
posterior = as.matrix(M1)
summary(posterior)
posterior2 = as.data.frame(posterior)
posterior2

# posterior predictive overlay
ppc_dens_overlay(y = M1$y,
                 yrep = posterior_predict(M1, draws = 50))

# posterior draws
M1 %>%
  posterior_predict(draws = 800) %>%
  ppc_stat_grouped(y = weight$log_wt,
                   group = weight$Diet,
                   stat = "median")

M1 %>%
  posterior_predict(draws = 800) %>%
  ppc_stat_grouped(y = weight$log_wt,
                   group = weight$Diet,
                   stat = "mean")

#posterior predictive intervals for my predictor
postpred_weight = ppc_intervals(y = weight$log_wt,
                             yrep = posterior_predict(M1),
                             x = weight$Diet, prob = 0.5) +
  labs(x = "Diet", y = "Log Weight",
       title = "50% Posterior Predictive Intervals") +
  panel_bg(fill = "gray", color = NA) +
  grid_lines(color = "white")


# predictions for every individual 
yrep = posterior_predict(M1)
meanyrep = colMeans(yrep) 
weight$yrep = meanyrep


# rsquared for regression model
rsq <- bayes_R2(M1)
print(median(rsq))

##############################################################################
#' simulate the model
#' setting parameters for the model simulation

# Parameters for simulated data
vy = 0.003   #Variance of y
va = 0.02   #Variance at a group level
mubar = 0.017 #mean among groups
ng = 7    #Number of groups in data
n = 151      #Number of data points within a group

# across group variance estimate
# at least my attempt
weight_est = weight %>%
  group_by(Diet) %>%
  summarise(mean(log_wt)) 
var(weight_est$`mean(log_wt)`)

w_sim = rnorm(ng, mubar, sqrt(vy)) 
w_sim
wfull = rnorm(n*ng,w_sim,sqrt(vy)) 
wfull
# compare distributions
hist(wfull)
hist(weight$log_wt)

# simulating the y and comparing models
mu = rnorm(ng,mubar,sqrt(va)) 
mu_v = rep(mu,each=n)        
y = rnorm(n*ng,mu_v,sqrt(vy)) 

weight_frame = data.frame(y=y,group=factor(rep(1:ng,each=n)))
weight_frame
weight_frame$wf = wfull 

# simulating model without predictors
# without predictors
bayesfit = stan_glmer(y ~ 1 + (1|group), data=weight_frame)

# simulating model with predictors
bayesfit_Pred = stan_glmer( y ~ wf + y:wf+ (1|group),
                             data=weight_frame)
summary(bayesfit)

# compare the simulated models to my model
# not that similar
simulation = print(summary(bayesfit)[,c(1,3,9,10)],digits=3)
simulation2 = print(summary(bayesfit_Pred)[,c(1,3,9,10)],digits=3)
originalmod = print(summary(M1)[,c(1,3,9,10)], digits=3)

##############
# build a model for melanization
#' Read in the data and explore the data structure
butterfly <- read.csv("Desktop/Butterfly_2018Data.csv", as.is = TRUE)
butterfly$Diet <- factor(butterfly$Diet)
str(butterfly)
summary(butterfly)

obsi <- 1:nrow(butterfly) #where df is the name of the data frame
# couldn't get the model to work, issue with the diet column
# the steps below is what I would do with the model, but just can't seem to figure this part out with the 
# diet, google didn't help, at least I didn't know what to search for...
M2 = stan_glmer(c(Mel,NoMel) ~ Diet + (1|obsi), family=binomial, data = butterfly)

# tried to remove diet and that did work
M3 = stan_glmer(c(Mel,NoMel) ~ Growth.Rate + (1|obsi), family=binomial, data = butterfly)

summary(M2,probs=c(0.025,0.975),digits=2)

print(summary(M2)[,c(1,3,9,10)],digits=3)

# assessing the model performance with shinystan
launch_shinystan(M2)

posterior_interval(M2, prob = 0.9, pars = "c(Mel, NoMel)")

# Extract posterior samples
samples = extract(M2$stanfit)

names(samples)
str(samples$b) 
str(samples$beta) 

# data visualtion of posterior
library(bayesplot)
posterior3 <- as.matrix(M2)
summary(posterior)
posterior4 <- as.data.frame(posterior3)
posterior4

# kind of disappointed i couldn't figure this out in time.
# I want to build a model where I could look at differences in melanization across diets
# eventually I would like to build a model where we can differences across time (generation)


