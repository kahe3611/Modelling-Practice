# individual project
# katherine hernandez

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
butterfly$Diet <- as.factor(butterfly$Diet)
butterfly$Type = as.factor(butterfly$Type)
str(butterfly)
summary(butterfly)

#' Plot:
ggplot() + geom_point(mapping=aes(x=Diet,y=Growth.Rate,col=Type),data=butterfly)

# look for any correlation across data set

# EDA look at any correlations across variables
corrplot(butterfly) 

# organizing the data a little bit
butterfly$Diet <- factor(butterfly$Diet, c("Dandelion","Lupine","Mallow","Plantago","MD1"))
#exploring data some more and looking at the distribution of each variable
# hemocyte count
hemoplot = ggplot() + geom_boxplot(data = butterfly, position = position_dodge(width = 0.75), aes(x = Diet, y = butterfly$Hemocyte, fill = Type), width = 0.5)
hemoplot
hemoplot + theme_classic() + scale_y_continuous(name = "Hemocyte Count", breaks = c(200,400,600,800,1000,1200,1400)) +
  theme(axis.text = element_text( size = 18)) +
  theme(axis.title = element_text(size = 24)) +
  theme(legend.text = element_text(size = 16))

# melanin as a proportion (out of 255 pixels)
melaplot = ggplot() + geom_boxplot(data = butterfly, position = position_dodge(width = 1), aes(x = Diet, y = Melanization, fill = Type), width = 0.5)
melaplot
melaplot + theme_classic() + scale_y_continuous(name = "Melanization (%)", breaks = c(0,10,20,30,40,50,60)) +
  theme(axis.text = element_text( size = 18)) +
  theme(axis.title = element_text(size = 24)) +
  theme(legend.text = element_text(size = 16))

# pupal weights
pupaplot = ggplot() + geom_boxplot(data = pupa, position = position_dodge(width = .75), aes(x = Diet, y = Weight, fill = Type), width = 0.5)
pupaplot
pupaplot + theme_classic() + scale_y_continuous(name = "Pupal Weight (g)", breaks = c(0,.10,.20,.30,.40,.50,.60)) +
  theme(axis.text = element_text( size = 18)) +
  theme(axis.title = element_text(size = 24)) +
  theme(legend.text = element_text(size = 16))

# growth rates, the overall growth rate for each surviving individual
growthplot + theme_classic() + scale_y_continuous(name = "Growth Rates", breaks = c(0,.10,.20,.30,.40,.50,.60)) +
  theme(axis.text = element_text( size = 18)) +
  theme(axis.title = element_text(size = 24)) +
  theme(legend.text = element_text(size = 16))

### look at distributions for predictors included in model
## doesn't work, running into an error
ggplot(melt(butterfly), aes(butterfly$Growth.Rate.Log) + geom_histogram(bins = 25) +
  facet_wrap(~variable, scales = "free")

# transform growth rates and look at distribution again
butterfly$Growth.Rate.Log = log(butterfly$Growth.Rate)


#' My first try at a Bayesian fit with `rstanarm`:
bayesHxL <- stan_glmer(Growth.Rate ~ Diet + Growth.Rate + Growth.Rate:Diet + (1|Type), 
                       family=binomial, data=butterfly)
summary(bayesHxL,probs=c(0.025,0.975),digits=2)

#' Notice that the fit takes a while. This may be partly because of the high
#' correlation between the intercept and latitude parameters:
vcov(bayesHxL,correlation=TRUE) #Correlation matrix

#' Indeed, if we try a maximum likelihood fit the algorithm fails to converge
#' (although the failed fit is not too bad - compare to Bayesian above):
lmerHxL <- glmer(richness ~ habitat + latitude + habitat:latitude + (1|site),
                 family=poisson, data=ant)
summary(lmerHxL) #Corr matrix here also highlights the problem with latitude

#' Scaling and centering latitude could improve computational performance:
ant$latitude_s <- scale(ant$latitude)

#' Repeat Bayesian fit with scaled/centered latitude. The correlation is now
#' reduced markedly, the fit is about 20% faster, and we have more effectively
#' independent samples.
bayesHxL <- stan_glmer(richness ~ habitat + latitude_s + habitat:latitude_s + (1|site), 
                       family=poisson, data=ant)
summary(bayesHxL,probs=c(0.025,0.975),digits=2)
vcov(bayesHxL,correlation=TRUE) #Correlation matrix

#' Inspect diagnostics of the fit
#+ eval=FALSE
launch_shinystan(bayesHxL)
#' In particular, the posterior distributions for the linear coefficients
#' $\beta_i$ are all nicely symmetric, as expected for this type of model.
#' 

#' Scaling and centering also fixes convergence for the maximum likelihood fit
lmerHxL <- glmer(richness ~ habitat + latitude_s + habitat:latitude_s + (1|site),
                 family=poisson, data=ant) 
summary(lmerHxL)


#' ### Model selection approaches
#' 
#' Now that we've got a fairly decent basic model working, we're good to move on
#' to some more complex investigations. We'll now contrast five ways to approach
#' model/variable selection, asking the question: Should the model include an
#' interaction of habitat and latitude? \

#' 1. Frequentist NHST \
#' 2. Bayesian credible interval for the interaction coefficient $\beta_3$ \
#' 3. Cross-validation \
#' 4. AICc \
#' 5. LOOIC\

#' Finally, in the next section, we'll look at a bigger investigation and
#' consider the many possible models that could result from considering
#' combinations of the predictors (including a new one: elevation) and all their
#' higher order interactions.
#' 

#' #### 1. Frequentist NHST
#'
#' In the classical frequentist approach we could conduct a null-hypothesis
#' significance test (NHST), which essentially asks: Could the data plausibly be
#' explained by a model without the interaction? More precisely, what is the
#' probability, *p*, of a test statistic, given $\beta_3 = 0$? If *p* < 0.05
#' reject the null hypothesis and decide to add the interaction to the model.
#'
#' First fit the model without the interaction. This will be the null model in
#' the NHST.
lmerHL <- glmer(richness ~ habitat + latitude_s + (1|site),
                family=poisson, data=ant)

#' The usual NHST for GLMMs is the likelihood ratio test. Twice the log of the
#' likelihood ratio
#'
#' $$ 2 ln ( \frac{L(full)}{L(null)} ) $$
#'
#' will be our test statistic. Asymptotically, this has a Chi-squared
#' distribution. We can thus compute *p* directly:
chisq_stat <- as.numeric( 2 * ( logLik(lmerHxL) - logLik(lmerHL) ) )
chisq_stat
1 - pchisq(chisq_stat,df=1) #p-value.

#' Or we can use the `anova()` convenience function, where the deviance is also
#' reported. Minus twice the log-likelihood is also called the deviance, so the
#' test statistic is equivalently the difference in deviance of the models.
anova(lmerHxL,lmerHL)

#' Since *p* is greater than 0.05, we don't find the data implausible given the
#' model without the interaction. So, in this approach we would decide to leave
#' it out.
#' 

#' #### 2. Bayesian credible interval for the interaction coefficient $\beta_3$
#'
#' In the Bayesian setting, we could ask if 0 is a credible value for the
#' interaction coefficient. Inspecting the posterior summary, we find the
#' estimate is close to zero and the 95% credible interval widely spans zero.
#' Here, we are using the central posterior interval as a more numerically
#' stable approximation to the HPDI since the posterior is symmetric.
summary(bayesHxL,pars="beta",probs=c(0.025,0.975),digits=4)

#' Since zero is well within the credible interval, we would decide to leave the
#' interaction out of the model. Of course, we could come to the same decision
#' from the frequentist confidence interval for the interaction, which is almost
#' the same as the Bayesian credible interval:
confint(lmerHxL)

#' #### 3. Cross-validation
#'
#' The cross-validation (CV) approach (and the related information criteria in
#' the following sections) have a different underlying philosophy. Whereas the
#' above approaches focus on inference for the value of the interaction
#' coefficient, here, the focus is on the predictive performance of the model on
#' new data (in other words, out-of-sample performance). Since CV is a general
#' tool, we can apply it to either the likelihood or Bayesian analysis. Here,
#' we'll go with the likelihood fit, since it is faster. We'll do
#' leave-one-out-CV (LOOCV), which is a special case of k-fold CV where k = n.
#' LOOCV is a good choice for the multilevel model because we won't make any
#' fitted model too unbalanced by leaving out only one data point (recall that
#' data points are paired at sites, one each of forest and bog). As we are
#' fitting a model with a log link, it makes sense to calculate the mean square
#' error (MSE) on the log scale (the scale of the linear predictor).

n <- nrow(ant) #number of data points
# Model with interaction
errHxL <- rep(NA,n)
for ( i in 1:n ) {
  lmerHxL_loo <- glmer(richness ~ habitat + latitude_s + habitat:latitude_s + (1|site),
                       family=poisson, data=ant[-i,])
  errHxL[i] <- ( predict(lmerHxL_loo,newdata=ant[i,]) - log(ant[i,"richness"]) )^2 #MSE
  rm(lmerHxL_loo)
}
CV_HxL <- mean(errHxL)
# Model without interaction
errHL <- rep(NA,n)
for ( i in 1:n ) {
  lmerHL_loo <- glmer(richness ~ habitat + latitude_s + (1|site),
                      family=poisson, data=ant[-i,])
  errHL[i] <- ( predict(lmerHL_loo,newdata=ant[i,]) - log(ant[i,"richness"]) )^2 #MSE
  rm(lmerHL_loo)
}
CV_HL <- mean(errHL)
# Compare models
cbind(CV_HxL,CV_HL)

#' We see that the prediction error is basically the same for the two models,
#' with slightly better out-of-sample predictive performance (lower LOOCV) for
#' the model without the interaction. Thus, we should prefer the model without
#' the interaction.
#' 

#' #### 4. AICc
#'
#' AIC can be derived as a measure of out-of-sample predictive performance. It
#' is very popular in ecology. Comparing the AIC of the two models is
#' straightforward:
AIC(lmerHL)
AIC(lmerHxL)
#' The model with the lower AIC should be preferred, so we again see that the
#' model without the interaction is best.
#'
#' Burnham and Anderson recommend a small sample correction to AIC when n/k <
#' 40, where n is the number of data points and k is the number of estimated
#' parameters. The corrected value is
#'
#' $$ AICc = AIC + \frac{2k(k+1)}{n-k-1} $$
#'
#' We can make a simple helper function that takes the fitted model object and
#' calculates AICc
AICc <- function(fitmod) {
  ll <- logLik(fitmod)
  k <-  attr(ll, "df")
  n <- attr(ll,"nobs")
  return( -2 * as.numeric(ll) + 2 * k + 2 * k * (k + 1) / (n - k - 1) )
}

#' Now we can use it:
AICc(lmerHL)
AICc(lmerHxL)
#' We see that this doesn't change our decision and indeed favors the simpler
#' model even more than the uncorrected AIC.
#' 

#' #### 5. LOOIC
#'
#' Finally, the Bayesian approach that is the state of the art, the
#' leave-one-out information criterion. It would be very computationally
#' expensive to do LOOCV directly on Bayesian models (it would require refitting
#' the model leaving out each observation!). The LOOIC approach approximates
#' this for most points but uses a direct CV when the approximation is not good
#' (typically zero to a few points).
#' 
#' First fit the model without the interaction term:
bayesHL <- stan_glmer(richness ~ habitat + latitude_s + (1|site), 
                      family=poisson, data=ant)

#' Now calculate the LOOIC, which is very conveniently done automatically for
#' `rstanarm` objects using the `loo()` function. First for the model without
#' the interaction:
loo(bayesHL)
#' We'll follow the advice here, which is telling us that the approximation
#' might not be accurate for one data point, so a direct LOOCV simulation is
#' warranted.
loo(bayesHL,k_threshold = 0.7)
#' Now for the model with the interaction:
loo(bayesHxL)

#' Comparing LOOIC between the models, we find that the model without the
#' interaction has lower LOOIC, so we favor the simpler model without the
#' interaction.
#' 

#' ### Comparing multiple Bayesian models with LOOIC
#'
#' Often you want to compare many models. This is computationally expensive but
#' not too bad. Here we compare the many possible models that could result from
#' considering combinations of the predictors (including a new one: elevation)
#' and all their higher order interactions. Model fitting takes about 20 secs
#' per model (5 mins in all) and about the same for the LOOIC algorithm. If we
#' found this too arduous (say if each of our models takes a long time to fit),
#' then we can fall back to maximum likelihood and AICc.
#' 

#' Scale and center the new variable
ant$elevation_s <- scale(ant$elevation)

#' Single factor models
bayesH <- stan_glmer(richness ~ habitat + (1|site), 
                     family=poisson, data=ant)
bayesL <- stan_glmer(richness ~ latitude_s + (1|site), 
                     family=poisson, data=ant)
bayesE <- stan_glmer(richness ~ elevation_s + (1|site), 
                     family=poisson, data=ant)

#' Two factor models
# We have already fitted: bayesHL <- habitat + latitude_s
bayesHE <- stan_glmer(richness ~ habitat + elevation_s + (1|site), 
                      family=poisson, data=ant)
bayesLE <- stan_glmer(richness ~ latitude_s + elevation_s + (1|site), 
                      family=poisson, data=ant)

#' Two factor models with interactions
# We have already fitted: bayesHxL <- habitat + latitude_s + habitat:latitude_s
bayesHxE <- stan_glmer(richness ~ habitat + elevation_s + habitat:elevation_s + (1|site), 
                       family=poisson, data=ant)
bayesLxE <- stan_glmer(richness ~ latitude_s + elevation_s + latitude_s:elevation_s + (1|site), 
                       family=poisson, data=ant)

#' Three factor model
bayesHLE <- stan_glmer(richness ~ habitat + latitude_s + elevation_s + (1|site), 
                       family=poisson, data=ant)

#' Three factor model with single two-way interactions
bayesHLE_HxL <- stan_glmer(richness ~ habitat + latitude_s + elevation_s + 
                             habitat:latitude_s + (1|site), 
                           family=poisson, data=ant)
bayesHLE_HxE <- stan_glmer(richness ~ habitat + latitude_s + elevation_s + 
                             habitat:elevation_s + (1|site), 
                           family=poisson, data=ant)
bayesHLE_LxE <- stan_glmer(richness ~ habitat + latitude_s + elevation_s + 
                             latitude_s:elevation_s + (1|site), 
                           family=poisson, data=ant)

#' Three factor model with multiple two-way interactions
bayesHLE_HxL_HxE <- stan_glmer(richness ~ habitat + latitude_s + elevation_s + 
                                 habitat:latitude_s + habitat:elevation_s + (1|site), 
                               family=poisson, data=ant)
bayesHLE_HxL_LxE <- stan_glmer(richness ~ habitat + latitude_s + elevation_s + 
                                 habitat:latitude_s + latitude_s:elevation_s + (1|site), 
                               family=poisson, data=ant)
bayesHLE_HxE_LxE <- stan_glmer(richness ~ habitat + latitude_s + elevation_s + 
                                 habitat:elevation_s + latitude_s:elevation_s + (1|site), 
                               family=poisson, data=ant)
bayesHLE_HxL_HxE_LxE <- stan_glmer(richness ~ habitat + latitude_s + elevation_s + 
                                     habitat:latitude_s + habitat:elevation_s + 
                                     latitude_s:elevation_s + (1|site), 
                                   family=poisson, data=ant)

#' Three factor model with three-way interaction (full model)
bayesHxLxE <- stan_glmer(richness ~ habitat + latitude_s + elevation_s + 
                           habitat:latitude_s + habitat:elevation_s + 
                           latitude_s:elevation_s + 
                           habitat:latitude_s:elevation_s + (1|site), 
                         family=poisson, data=ant)

#' Calculate LOOICs. These simulations will take a while. We include
#' k_threshold=0.7 because many simulations had one or two observations where an
#' explicit LOOCV was needed. Because these are stochastic simulations,
#' sometimes it will be needed and sometimes not, so it's easiest to include the
#' option for all simulations.
loo_H <- loo(bayesH,k_threshold = 0.7)
loo_L <- loo(bayesL,k_threshold = 0.7)
loo_E <- loo(bayesE,k_threshold = 0.7)

loo_HL <- loo(bayesHL,k_threshold = 0.7)
loo_HE <- loo(bayesHE,k_threshold = 0.7)
loo_LE <- loo(bayesLE,k_threshold = 0.7)

loo_HxL <- loo(bayesHxL,k_threshold = 0.7)
loo_HxE <- loo(bayesHxE,k_threshold = 0.7)
loo_LxE <- loo(bayesLxE,k_threshold = 0.7)

loo_HLE <- loo(bayesHLE,k_threshold = 0.7)

loo_HLE_HxL <- loo(bayesHLE_HxL,k_threshold = 0.7)
loo_HLE_HxE <- loo(bayesHLE_HxE,k_threshold = 0.7)
loo_HLE_LxE <- loo(bayesHLE_LxE,k_threshold = 0.7)

loo_HLE_HxL_HxE <- loo(bayesHLE_HxL_HxE,k_threshold = 0.7)
loo_HLE_HxL_LxE <- loo(bayesHLE_HxL_LxE,k_threshold = 0.7)
loo_HLE_HxE_LxE <- loo(bayesHLE_HxE_LxE,k_threshold = 0.7)

loo_HLE_HxL_HxE_LxE <- loo(bayesHLE_HxL_HxE_LxE,k_threshold = 0.7)
loo_HxLxE <- loo(bayesHxLxE,k_threshold = 0.7)

#' Now we can compare the models, looking particularly at the LOOIC column.
compare_models(loo_H,loo_L,loo_E,loo_HL,loo_HE,loo_LE,loo_HxL,loo_HxE,loo_LxE,
               loo_HLE,loo_HLE_HxL,loo_HLE_HxE,loo_HLE_LxE,
               loo_HLE_HxL_HxE,loo_HLE_HxL_LxE,loo_HLE_HxE_LxE,
               loo_HLE_HxL_HxE_LxE,loo_HxLxE)
#' We see that the model with the best out-of-sample prediction error is
#' bayesHLE, the model with each of the variables but with no interactions.
#' Notice the standard errors of LOOIC - they are all large compared to the
#' LOOIC differences between models. Thus, we shouldn't get too fussed about
#' small differences in LOOIC between models. If we ran the simulations again,
#' we might get a slightly different ordering of models. A LOOIC difference of
#' about 2 indicates a fairly negligible difference between models. So, the top
#' four models all have about the same predictive performance.
#'
#' A maximum likelihood analysis followed by AICc gives very similar results,
#' identifying the same "best" model, then the next top three, with some
#' differences in the ordering of other models.
