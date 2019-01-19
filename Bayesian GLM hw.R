# katherine hernandez
# week 7
# Frequentist GLM homework

#load in any packages needed
library(tidyverse)

# reading in data from my desktop
ant <- read.csv("Desktop/ants.csv", as.is = TRUE)

# check the dataset 
head(ant)
class(ant$site)
str(ant)
summary(ant)

# make factors for habitat and site
ant$habitat <- factor(ant$habitat)
ant$site <- factor(ant$site)

# plot data and look at the distribution
windows()
ggplot(data = ant) + 
  geom_point(mapping = aes(x = latitude, y = richness, color = habitat)) +
  geom_smooth(mapping = aes(x = latitude, y = richness, color = habitat)) +
  labs( x = "Latitude", y = "Ant Richness")

#work through frequentist glm code
ant_fit <- glm( richness ~ latitude + habitat + latitude:habitat, family = poisson(link = "log"), data = ant )

#looks at standard error and estimates
summary(ant_fit) 

#confidence intervals 
confint(ant_fit) 

#make correlation matrix for the parameters
cov2cor(vcov(ant_fit)) 

#The log likelihood
logLik(ant_fit) 

# create some new data to construct CI with
new_ant <- data.frame(latitude = rep(seq(min(ant$latitude), max(ant$latitude), length.out=100),2), 
                      habitat = factor(rep(c("forest","bog"), each=100)))

# use SE * 2 to create CI
preds <- predict(fit,newdata=new_ant,se.fit=TRUE)

antlp <- preds$fit 
#mean of the linear predictor

antse <- preds$se.fit     
#se of the linear predictor

lowerant <- antlp - 2 * antse  
#lower of 95% CI for linear predictor

upperant <- antlp + 2 * antse  
#upper 95% CI 

lowerant_exp <- exp(lowerant)  
#lower of 95% CI for response 

upperant_exp <- exp(upperant)       
#upper 95% CI for response

antlp_exp <- exp(mlp)           
#mean of response scale

# predictions
preds <- cbind(new_ant,preds,lowerant,upperant,antlp)
preds


# going to use ggplot to plot the data
windows()
ggplot(data = ant) + 
  geom_point(mapping = aes(x = latitude, y = richness, color = habitat)) +
  geom_smooth(mapping = aes(x = latitude, y = richness, color = habitat), method = "glm") +
  labs( x = "Latitude", y = "Ant Richness")


#Attempt at a normal model for the ant data
antfit2 <- lm( richness ~ latitude + habitat + latitude:habitat, data = ant )

#provides a summary of estimates and standard error, among other things....important things
summary(antfit2) 

#confidence intervals for the parameters
confint(antfit2)

#make a correlation matrix for the parameters.
cov2cor(vcov(antfit2)) 

#use log liklihood function on normal model
logLik(antfit2)     

# created new data frame for the ant data
new_ant2 <- data.frame(latitude = rep(seq(min(ant$latitude), max(ant$latitude), length.out=100),2),habitat = factor(rep(c("forest","bog"), each=100)))    
# predict new values
ant_pred <- predict(antfit2, newdata=new_ant2, interval = "confidence", se.fit=TRUE)

# model fit is okay, not sure if I got this part right at all....
windows()
plot(ant$latitude, ant$richness)
ant_lty <- c(2,1)
category <- c("forest","bog")
#making a for loop with the data
for ( i in 1:2 ) {
  lines(new_ant2$latitude, ant_pred$fit[,1], lty=ant_lty[i])
  lines(new_ant2$latitude, ant_pred$fit[,2], lty=ant_lty[i], col="blue")
  lines(new_ant2$latitude, ant_pred$fit[,3], lty=ant_lty[i], col="blue")
}


# made a plot of log transformed ant data
plot(log(ant$latitude), log(ant$richness))

#model with log transformed data, should see differences in distribution versus the backtransformed ant data
antfit3 <- lm( log(richness) ~ log(latitude) + habitat + log(latitude):habitat, data = ant )

# looked at the summary statistics
summary(antfit3)

#wanted to check out the diagnostics of the model
plot(antfit3,1:6,ask=FALSE) 

#####BAYESIAN ATTEMPT##########
#building off of the frequentist approach with the ant data and now will try to work through my attempts at a bayesian approach
# frequentist code is above, sorry for the need to scroll, needed it for reference.
# packages needed
library(ggplot2)
library(dplyr) #doesn't work for me, I keep getting an error with the R_Lang function
library(tidyr)
library(rstan) #to extract() samples (nb conflict with tidyr)
library(rstanarm) #nb function se over-rides rethinking version, still not sure what this does though...
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
theme_set(theme_grey()) #rstanarm overrides default ggplot theme: sets it back
source("source/hpdi.R") #calculating credible intervals

ants <- read.csv("Desktop/ants.csv", as.is = TRUE)
#explore the data and look at data frame
head(ants)
ants
str(ants)
class(ants$site)
unique(ants$site)

#make site and habitat into factors
ants$habitat <- factor(ants$habitat)
ants$site <- factor(ants$site)

# trying to construct a bayesian fit of the model, no transformation of the data
bayesmod <- stan_glm(richness ~ habitat + latitude + habitat:latitude,family=poisson, data=ant)
summary(bayesmod)
coefficients(bayesmod)
prior_summary(bayesmod)

# making draws of parameters from model
draws <- as.matrix(bayesmod)

# calculating a correlation matrix with vcov
vcov(bayesmod,correlation=TRUE) #Correlation matrix

# check out distributions for parameters of the model
hist(draws[,1], breaks = 50)
hist(draws[,2], breaks = 50)
hist(draws[,3], breaks = 50)
hist(draws[,4], breaks = 50)

# i'm going to try a different model and try my hand at centering and log transforming the data
bayesmod2 <- stan_glm(richness ~ habitat * scale(latitude), family = poisson(link="log"), data = ants)
prior_summary(bayesmod2)
coefficients(bayesmod2)
summary(bayesmod2)
#looks good I think

# central posterior intervals for both models
posterior_interval(bayesmod,prob=0.95) 
posterior_interval(bayesmod2,prob=0.95)

# compare model estimates
summary(bayesmod)
summary(bayesmod2)

# setting up priors for each model, not sure where to go from here though...
# do you try different combinations with the variables and create additional models?
# how do you know if 2, 3, 4, etc. models is enough for a dataset to make a model selection?
prior1 <- normal(location = 0, scale = .05, autoscale = FALSE)
prior2 <- normal(location = 0, scale = 4, autoscale = FALSE)


