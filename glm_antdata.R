#' ### Ant data Generalized Linear Model - Bayesian

#' Second in a series of scripts to analyze the ant data described in Ellison
#' (2004). This script includes plotting variations (ggplot & base) and Bayesian
#' inference from the GLM. Future scripts will consider multilevel models to
#' fully account for the design structure.
#' 
#' This script is `knitr` ready. Use `spin()` for reproducible report.\
#' This script will likely be updated a few times also.
#' 
#' Brett Melbourne\
#' 16 Oct 2018
#' 

#' Set up for Bayesian analysis (order is important):
#+ warning=FALSE

library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan) #to extract() samples (nb conflict with tidyr)
library(rstanarm) #nb function se over-rides rethinking version
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
theme_set(theme_grey()) #rstanarm overrides default ggplot theme: set it back
source("source/hpdi.R") #For calculating credible intervals

#' Read in the data:
ant <- read.csv("Desktop/ants.csv", as.is = TRUE)
ant$habitat <- factor(ant$habitat)
ant$site <- factor(ant$site)

#' Plot:
ggplot() +
  geom_point(mapping=aes(x=latitude,y=richness,col=habitat),data=ant)

#' Bayesian fit with `rstanarm`:
bysfitHxL <- stan_glm(richness ~ habitat + latitude + habitat:latitude, 
                      family=poisson, data=ant)
summary(bysfitHxL,digits=4)
vcov(bysfitHxL,correlation=TRUE) #Correlation matrix


#' #### Working with posterior samples

#' We'll first work with the samples directly (as we did in McElreath Ch 4).
#' There are also convenience functions to do many standard things but you will
#' often want to calculate new quantities from the samples directly to answer
#' science questions that aren't addressed by the standard output, so this is
#' important to know.

samples <- extract(bysfitHxL$stanfit)
class(samples)
str(samples)
names(samples)

#' #### Diagnostics

#' Histograms directly from the samples
hist(samples$alpha[,1],breaks=75,col="lightblue") #\beta_0 = intercept = alpha 
hist(samples$beta[,1],breaks=75,col="lightblue") #\beta_1
hist(samples$beta[,2],breaks=75,col="lightblue") #\beta_2
hist(samples$beta[,3],breaks=75,col="lightblue") #\beta_3

#' To do the same with `ggplot`, we first need the samples in a dataframe and
#' then in tidy format (using `gather`). One could make a function for this that
#' works generically on extracted samples. No doubt, there is a package that
#' does this.
samplesdf <- data.frame(samples$alpha,samples$beta)
names(samplesdf) <- c("alpha",paste(names(samples[2]),1:3,sep="_"))
samplesdf %>% 
  gather(key = "parameter", value = "sample") %>%
  ggplot() +
  geom_histogram(mapping = aes(x=sample,y=stat(density)),bins=75) +
  facet_wrap(facets = ~ parameter,scales="free")

#' Manual trace plot from the samples\
#' e.g for beta_0
intercept_trace <- extract(bysfitHxL$stanfit,pars="alpha",permuted=FALSE,inc_warmup=TRUE)
plot(NA,NA,type="n",ylim=range(intercept_trace),xlim=c(0,length(intercept_trace[,1,1])))
for ( i in 1:4 ) {
  lines(intercept_trace[,i,1],col=i)  
}
#' We see that the chains had converged very quickly.


#' #### Inference: 95% credible intervals (HPDI)

#' Directly from the samples:
hpdi(samples$alpha[,1], prob=0.95)
hpdi(samples$beta[,1], prob=0.95)
hpdi(samples$beta[,2], prob=0.95)
hpdi(samples$beta[,3], prob=0.95)

#' These are almost the same as the CPIs due to the symmetric posteriors. Here
#' is the convenience function for the parameter CPIs:
posterior_interval(bysfitHxL,prob=0.95)
#' There is a good argument for using the CPI here in lieu of the HPDI. We saw
#' from the histograms that the posteriors are quite symmetric but also that
#' there was still some noise in the tails of the samples. Thus, the CPI is
#' probably a more numerically stable estimate of the credible interval, even
#' though the CPI is not a credible interval itself.


#' #### Mean curves, regression intervals (HPDI), posterior predictive dist

#' Directly from the samples:\
#' The following is quite literal. We could make this more elegant but the steps
#' needed are clear this way.

# Initialize variables and storage 
latitude <- seq( from=41.9, to=45, length.out=50 ) #range for latitude
n <- length(latitude)
hpdi_bog <- matrix(NA,nrow=n,ncol=5) #to store hpdi values and mean
colnames(hpdi_bog) <- c("mnmu","mulo95","muhi95","ppdlo95","ppdhi95")
hpdi_forest <- matrix(NA,nrow=n,ncol=5)
colnames(hpdi_forest) <- c("mnmu","mulo95","muhi95","ppdlo95","ppdhi95")

# For each latitude, form the posterior
for ( i in 1:n ) {
  
  # First form samples for the linear predictor \eta
  eta_bog <- samples$alpha[,1] + 
    samples$beta[,2] * latitude[i]
  eta_forest <- samples$alpha[,1] + 
    samples$beta[,1] + 
    samples$beta[,2] * latitude[i] + 
    samples$beta[,3] * latitude[i]
  
  # Then use inverse link for samples of the posterior \mu
  mu_bog <- exp(eta_bog)
  mu_forest <- exp(eta_forest)
  
  # Sample from Poisson to get the posterior predictive distribution
  ppd_bog <- rpois(n=length(mu_bog),lambda=mu_bog)
  ppd_forest <- rpois(n=length(mu_forest),lambda=mu_forest)
  
  # Mean and intervals of these samples
  hpdi_bog[i,1] <- mean( mu_bog )
  hpdi_bog[i,2:3] <- hpdi( mu_bog, prob=0.95 )
  #hpdi_bog[i,4:5] <- hpdi( ppd_bog, prob=0.95 )
  hpdi_bog[i,4:5] <- quantile( ppd_bog, prob=c(0.025,0.975) ) #CPI
  hpdi_forest[i,1] <- mean( mu_forest )
  hpdi_forest[i,2:3] <- hpdi( mu_forest, prob=0.95 )
  #hpdi_forest[i,4:5] <- hpdi( ppd_forest, prob=0.95 )
  hpdi_forest[i,4:5] <- quantile( ppd_forest, prob=c(0.025,0.975) ) #CPI
  
}
#' Notice that we calculated expectations (means) and intervals directly on the
#' scale of the data (the "response" scale), not on the linear predictor scale.
#' If we calculated first on the linear predictor scale and then backtransformed
#' the intervals to the response scale they would be biased due to nonlinear
#' averaging. Also, the posterior predictive distribution (PPD) can, of course,
#' only be on the response scale. I used the CPI (`quantile()`) for the
#' posterior predictive distribution because the HPDI was not numerically
#' stable.
#'

#' Package in tidy format for plotting:
mcpreds_df <- data.frame(habitat=rep("bog",n),latitude,hpdi_bog)
mcpreds_df <- rbind(mcpreds_df,data.frame(habitat=rep("forest",n),latitude,hpdi_forest))
rm(latitude,n,hpdi_bog,hpdi_forest,eta_bog,eta_forest,mu_bog,mu_forest) #clean up

#' Plot the inference\
#' I have done a bit of customization for colors and labels. The credible
#' intervals for the means are the shaded regions while the dashed lines show
#' the posterior predictive interval.
ggplot() +
  geom_ribbon(mapping=aes(x=latitude,ymin=mulo95,ymax=muhi95,fill=habitat),
              alpha=0.2,show.legend=FALSE,data=mcpreds_df) +
  geom_point(mapping=aes(x=latitude,y=richness,col=habitat),
             show.legend=FALSE,data=ant) +
  geom_line(mapping=aes(x=latitude,y=mnmu,col=habitat),
            show.legend=FALSE,data=mcpreds_df) +
  geom_line(mapping=aes(x=latitude,y=ppdlo95,col=habitat),lty=2,
            show.legend=FALSE,data=mcpreds_df) +
  geom_line(mapping=aes(x=latitude,y=ppdhi95,col=habitat),lty=2,
            show.legend=FALSE,data=mcpreds_df) +
  geom_text(mapping=aes(x=42.9,y=3.3,label="Bog"),col="#d95f02") +
  geom_text(mapping=aes(x=43.85,y=9.5,label="Forest"),col="#1b9e77") +
  scale_fill_manual(values = c("#d95f02","#1b9e77")) +
  scale_color_manual(values = c("#d95f02","#1b9e77")) +
  ylim(0,20) +
  xlab("Latitude (degrees north)") +
  ylab("Ant species richness") 

#' Notice that the intervals for forest are wider than for bog. This is because
#' the uncertainty scales with the mean of the response. Also notice that the
#' intervals for the posterior predictive distribution have discrete steps. This
#' is because the data generating process is discrete (e.g. we cannot have 1.3
#' species). Also, there are a few blips in the intervals for the predictive
#' distribution and some wiggles in the mean intervals. This is due to Monte
#' Carlo error. Increase the number of iterations when training the model (e.g.
#' iter=10000) and these will largely go away.


#' #### Using convenience functions

#' Extracting everything manually from the samples is a fair bit of coding work.
#' The convenience functions in `rstanarm` make this easier for common tasks. As
#' in `glm()`, these functions take a `newdat` argument that simplifies coding.
#' To do what we just did manually above, first make a dataframe with the
#' desired values of the explanatory variables.
newd <- data.frame(latitude = rep(seq(min(ant$latitude), max(ant$latitude), 
                                      length.out=50),2),
                   habitat = factor(rep(c("bog","forest"),each=100)))

#' Then derive samples for the posterior distribution of the inverse link
#' function, i.e. Dist($\mu$), which we'll call `pmu`.
pmu <- posterior_linpred(bysfitHxL, transform = TRUE, newdata = newd)

#' This is a matrix with samples in rows and the variable combinations in
#' columns. The estimated means are then:
mnmu <- colMeans(pmu)
#' and the credible intervals for the mean are:
regression_intervals <- t(apply(pmu,2,hpdi))
colnames(regression_intervals) <- c("mulo95","muhi95")

#' For predictions, first derive samples for the posterior predictive
#' distribution, which we'll call ppd:
ppd <- posterior_predict(bysfitHxL, newdata = newd)
#' and the prediction intervals (here CPI) are then:
prediction_intervals <- t(apply(ppd,2,quantile,prob=c(0.025,0.975)))
colnames(prediction_intervals) <- c("ppdlo95","ppdhi95")

#' Plot
mcpreds_df <- cbind(newd,mnmu,regression_intervals,prediction_intervals)
ggplot() +
  geom_ribbon(mapping=aes(x=latitude,ymin=mulo95,ymax=muhi95,fill=habitat),
              alpha=0.2,show.legend=FALSE,data=mcpreds_df) +
  geom_point(mapping=aes(x=latitude,y=richness,col=habitat),
             show.legend=FALSE,data=ant) +
  geom_line(mapping=aes(x=latitude,y=mnmu,col=habitat),
            show.legend=FALSE,data=mcpreds_df) +
  geom_line(mapping=aes(x=latitude,y=ppdlo95,col=habitat),lty=2,
            show.legend=FALSE,data=mcpreds_df) +
  geom_line(mapping=aes(x=latitude,y=ppdhi95,col=habitat),lty=2,
            show.legend=FALSE,data=mcpreds_df) +
  geom_text(mapping=aes(x=42.9,y=3.3,label="Bog"),col="#d95f02") +
  geom_text(mapping=aes(x=43.85,y=9.5,label="Forest"),col="#1b9e77") +
  scale_fill_manual(values = c("#d95f02","#1b9e77")) +
  scale_color_manual(values = c("#d95f02","#1b9e77")) +
  ylim(0,20) +
  xlab("Latitude (degrees north)") +
  ylab("Ant species richness") 

plot(bysfitHxL)
bysfitHxL$coefficients
#' Hilariously, the `rstanarm` team were originally not even going to supply the
#' `posterior_linpred()` function and they [deliberately made it very hard to
#' find](https://github.com/stan-dev/rstanarm/issues/85), arguing for the
#' superiority of the `posterior_predict()` function. Thus, we can only use
#' `posterior_linpred()` because "we know what we are doing". However, that
#' stance is baffling. These two different functions serve completely different
#' needs: `posterior_predict()` allows us to calculate uncertainty for the
#' prediction of new observations, whereas `posterior_linpred()` allows us to
#' calculate uncertainty of the mean. Furthermore, `posterior_linpred()` should
#' be preferred to estimate means because it has less Monte Carlo error compared
#' to `posterior_predict()`. Comparing frequentist and Bayesian:

#' |        |Mean                     |Uncertainty of mean      |Uncertainty of prediction|
#' |:-------|:------------------------|:------------------------|:------------------------|
#' |lm      |predict()                |predict(int="confidence")|predict(int="prediction")|
#' |glm     |predict(type= "response")|predict(se.fit=TRUE)     |via bootstrap            |
#' |        |                         |or via bootstrap         |                         |
#' |stan_glm|mean(pmu)                |hpdi(pmu)                |hpdi(ppd)                |
#'
#' where:\
#' pmu <- posterior_linpred(transform = TRUE) \
#' ppd <- posterior_predict()
#'
#' In science, inference about the mean is often the goal rather than prediction
#' of new observations.