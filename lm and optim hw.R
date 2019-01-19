#caterpillar data from summer of 2018
# 1) Individual caterpillars identified by a number
#2) Diet
#3) Hemocyte count
#4) Growth Rates
#5) Melanization as a proportion (out of 255 pixels)
#6) pupal weights
## reading in the data

butterfly <- read.csv("Desktop/Butterfly_2018Data.csv", as.is = TRUE)
str(butterfly)
butterfly$logGrowth.Rate = log(butterfly$Growth.Rate)

## Explore variables

hist(butterfly$Growth.Rate)
hist(butterfly$Hemocyte)
hist(butterfly$logGrowth.Rate)
hist(butterfly$PupalWeight)
hist(butterfly$Hemocyte)

# Fitting model

Blm <- lm(Growth.Rate ~ Diet, butterfly)
summary(Blm)
x<- seq(min(butterfly$Growth.Rate), max(butterfly$Growth.Rate), length.out = 1000)



# Optim using SSQ 

# make a likelihood function
ssq.func <- function(params, x, y_obs){
  pred <- params[1] + params[2]*x
  ssq.sum <- sum((pred-y_obs)^2)
  return(ssq.sum)
}
# test it
ssq.func(params = c(.01, 0.1), x = butterfly$Hemocyte, y_obs = butterfly$Melanization)


# Optim model using log-likelihood

design_Full <- model.matrix(~Growth.Rate, data = butterfly)
design_Full
likeFunc <- function(params, X, y){
  k <- ncol(X)
  b <- params[1:k]
  sigma <- params[k+1]
  yhat <- X %*% b
  ll <- sum(dnorm(y, yhat, sigma, log=T))
  return(-ll)
}

#
fit_FullB <- optim(c(100,20,20), likeFunc, X=butterfly$Hemocyte, y=butterfly$Melanization, method = "BFGS", hessian=T)
fit_FullB$par
#standard error and t value calculations
SEs <- sqrt(diag(solve(fit_FullB$hessian)))
SEs
tval <- qt(0.025)
lower <- fit_FullN$par + tval*SEs[]
upper <- fit_FullN$par - tval*SEs[]
lower
upper
confint(Nlm)
params <- fit_FullN$par
Sigma <- solve(fit_FullN$hessian)
# attempts at making confidence interval
library(mvtnorm)
xPred = with(butterfly, seq(min(Hemocyte), max(Hemocyte), len=1000))
post_draw <- function(){
  post_draw <- rmvnorm(1, params, Sigma)
  yhat <- post_draw[1]+post_draw[2]*xPred
  return(yhat)
}
post_preds <- replicate(1000, post_draw())
lower <- apply(post_preds, 1, function(x) quantile(x, 0.025))
upper <- apply(post_preds, 1, function(x) quantile(x, 0.975))
{plot(Hemocyte~Melanization, data=butterfly, pch=14)
  polygon(c(xPred, rev(xPred)), c(lower, rev(upper)),
          border=NA, col="lightgray")
  yhat = params[1] + params[2]*xPred
  lines(xPred, yhat, col="red", lw=2)
  points(Hemocyte~Melanization, butterfly, pch=16)}



