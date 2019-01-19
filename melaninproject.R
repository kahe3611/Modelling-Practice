library(ggplot2)
library(dplyr)
library(lme4)
library(reshape2)
library(arm)
library(rstan) #for extract function
library(rstanarm)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
theme_set(theme_classic())

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
posterior <- as.matrix(M2)
summary(posterior)
posterior2 <- as.data.frame(posterior)
posterior2