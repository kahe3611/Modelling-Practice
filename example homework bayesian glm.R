ants <- read.csv("ants.csv", as.is = TRUE)
head(ants)
ants
class(ants$site)
unique(ants$site)
ants$habitat <- factor(ants$habitat)
ants$site <- factor(ants$site)
theme_set(theme_grey()) # re set ggplot defualt since rstan will change 
```


# Running the model
Model 1: default priors, no centering
The priors are: intercept is normally distributed around 0, sd = 10
Slopes are centered around 0 with sd of 2.5 but are adjusted then based on each predictor's scale
Not sure what the prior is for SD (or does it not exist??)
```{r, echo = FALSE}
fit_b <- stan_glm(richness ~ habitat + latitude + habitat:latitude,
family= poisson(link = "log"), data = ants)
prior_summary(fit_b)
coefficients(fit_b)
summary(fit_b)
```
Model 2: center variables
```{r, echo = FALSE}
fit_b2 <- stan_glm(richness ~ habitat * scale(latitude), family = poisson(link="log"), data = ants)
prior_summary(fit_b2)
```

Model 3: center variables, add more informative priors
```{r}
my_prior <- c(normal(location = c(0,0,0), scale = c(5,1,.1), autoscale = FALSE))
fit_b3 <- stan_glm(richness ~ habitat * scale(latitude), family = poisson(link="log"), data = ants,
prior = my_prior)
prior_summary(fit_b3)
summary(fit_b3)
my_prior1 <- normal(location = 0, scale = .01, autoscale = FALSE)
my_prior2 <- normal(location = 0, scale = 5, autoscale = FALSE)
my_prior3 <- normal(location = 0, scale = 5, autoscale = FALSE)
fit_b4 <- stan_glm(richness ~ habitat * scale(latitude), family = poisson(link="log"), data = ants,
prior =c( my_prior1, my_prior2, my_prior3))
summary(fit_b4)
summary(fit_b)
```

# Model output exploration
Compare model estimates
```{r}
summary(fit_b)
summary(fit_b2)
summary(fit_b4)
```

Sample the posterior.
4 chains for each parameter
X axis is iteration (2000 iters)
Y is the parameter value
The first half of the iteration are dumped (gets rid of a lot of good data so it's pretty conservative)
The samples we get from posterior sampling are shuffled (so if we take 10 we don't get the first 10)
                                                         
                                                         We have 4000 samples from each.
                                                         One row = 1 draw of each parameter (they all belong together; it's not shuffling a beta1 with a beta2 from a different iteration)
Think of the posterior as being multivariate

```{r}
samples_1 <- extract(fit_b$stanfit)
median(samples_1$beta[,1])
fit_b$coefficients # matches
# we have four chains with 
samples_2 <- extract(fit_b2$stanfit)
hist(samples_1$beta[,1]) #e.g. histogram of \beta_1
hist(samples_1$beta[,2]) #e.g. histogram of \beta_1
hist(samples_1$beta[,3]) #e.g. histogram of \beta_1
# Convenience functions that estimate things from the samples
methods(class="stanreg") 
posterior_interval(fit_b,prob=0.95) #central posterior intervals, nb default=0.9
posterior_interval(fit_b2,prob=0.95) #central posterior intervals, nb default=0.9
posterior_interval(fit_b3,prob=0.95) #central posterior intervals, nb default=0.9
vcov(fit_b,correlation=TRUE) #Correlation matrix
vcov(fit_b2,correlation=TRUE) #Correlation matrix
```
# Plotting (from class)
Problem: have posteriors for params
But, we want a posterior for the LINE (need posterior for every point on the line)

Make a dataframe that has forest AND bog in it and a range of different latitudes
With each of the 4000 samples, calculate mu i for each row of the dataframe
We have 4000 mu's (one for every value of latitude). Have a posterior of mu for every latitude and habitat
Now we can take a interval

Make the new dataframe
```{r}
newdata1 <- data.frame(
habitat = c(rep("bog", 1000), rep("forest",1000)),
latitude = rep(seq(min(ants$latitude), max(ants$latitude), len = 1000),2)
)
```
Get expected values from the newdataframe
```{r}
posterior.fits <- posterior_linpred(fit_b,newdata = newdata1, transform = TRUE)
hdpi.fits <- data.frame(lo = rep(NA,2000), hi = rep(NA,2000))
for(i in 1:2000){
hdpi.fits$lo[i] <-HPDI(posterior.fits[,i], 0.89)[1]
hdpi.fits$hi[i] <-HPDI(posterior.fits[,i], 0.89)[2]
hdpi.fits$med[i] <-median(posterior.fits[,i])
}
newdata$cilo <- hdpi.fits$lo
newdata$cihi <- hdpi.fits$hi
newdata$med <- hdpi.fits$med
?posterior_predict
postpred.1 <- predictive_interval(fit_b, newdata, prob = c(0.89))
newdata$hi <- postpred.1[,1]
newdata$lo <- postpred.1[,2]
mycols <- c("saddlebrown", "forestgreen")
ggplot()+
geom_ribbon(data = newdata, aes(x=latitude, ymin =lo, ymax =hi, fill = habitat), alpha = 0.2)+
geom_ribbon(data = newdata, aes(x = latitude, ymin = cilo, ymax =cihi, fill = habitat), alpha = 0.4)+
geom_line(data = newdata, aes(x=latitude, y = med, col =habitat))+
geom_point(data = ants, aes(x=latitude, y = richness, col = habitat))+
scale_fill_manual(values = mycols)+
scale_color_manual(values=mycols)+
labs(x="latitude", y = "species richness")
```

Now try plotting with the centered model
```{r}
# go from -2 to +2 sds in latitude
#fit_b2 <- stan_glm(richness ~ habitat * scale(latitude), family = poisson(link="log"), data = ants)
newdata2 <- data.frame(
habitat = c(rep("bog", 1000), rep("forest",1000)),
latitude = rep(seq(min(ants$latitude), max(ants$latitude), len = 1000),2)
)
posterior.fits2 <- posterior_linpred(fit_b2,newdata = newdata1, transform = TRUE)
hdpi.fits2 <- data.frame(lo = rep(NA,2000), hi = rep(NA,2000))
for(i in 1:2000){
hdpi.fits2$lo[i] <-HPDI(posterior.fits[,i], 0.89)[1]
hdpi.fits2$hi[i] <-HPDI(posterior.fits[,i], 0.89)[2]
hdpi.fits2$med[i] <-median(posterior.fits[,i])
}
newdata2$cilo <- hdpi.fits2$lo
newdata2$cihi <- hdpi.fits2$hi
newdata2$med <- hdpi.fits2$med
postpred.2 <- predictive_interval(fit_b2, newdata2, prob = c(0.89))
newdata2$hi <- postpred.2[,1]
newdata2$lo <- postpred.2[,2]
ggplot()+
geom_ribbon(data = newdata, aes(x=scale(latitude), ymin =lo, ymax =hi, fill = habitat), alpha = 0.2)+
geom_ribbon(data = newdata, aes(x = scale(latitude), ymin = cilo, ymax =cihi, fill = habitat), alpha = 0.4)+
geom_line(data = newdata, aes(x=scale(latitude), y = med, col =habitat))+
geom_point(data = ants, aes(x=scale(latitude), y = richness, col = habitat))+
scale_fill_manual(values = mycols)+
scale_color_manual(values=mycols)+
labs(x="latitude", y = "species richness")
```

# Plotting
See previous file (ants_data) for the unscaled plots

```{r}
# make a dataframe of the output
fits <- fit_b2 %>% 
as_data_frame %>%
rename(intercept = `(Intercept)`) 
# get the HDPI for the parameters.
hdpi.1 <- posterior_interval(fit_b2, prob = 0.89)
# get predictor dataframes for forest and bog, we will predict richness for each draw from the posterior
newd.f <- data.frame(
intercept = rep(1, 100),
habitat = rep(1,100),
latitude = seq(min(ants$latitude), max(ants$latitude), len=100)
)
newd.f$for.lat <- newd.f$habitat*newd.f$latitude
newd.b <- data.frame(
intercept = rep(1, 100),
habitat = rep(0,100),
latitude = seq(min(ants$latitude), max(ants$latitude), len=100)
)
newd.b$b.lat <- newd.b$habitat*newd.b$latitude
newd.b$scalelat <- scale(newd.b$latitude)
# make a vector of the "best" parameters
params.b <- fit_b2$coefficients %>% data.frame()
# predict richness using the "best" parameters
# using matrix math here...is there a way to actually use a predict function to simplify this?
newd.b$best <- as.vector(exp(as.matrix(newd.b[1:4]) %*% as.matrix(params.b)))
newd.f$best <- as.vector(exp(as.matrix(newd.f[1:4]) %*% as.matrix(params.b)))
# now we want to predict a richness using each draw from the posterior
# vectors to store these in (1 for forest, 1 for bog)
mus.f <- matrix(NA, nrow=100, ncol = 1000)
mus.b <- matrix(NA, nrow = 100, ncol = 1000)
lats = seq(min(ants$latitude), max(ants$latitude), len=100)
# vectors to store the confidence interval (I think that's what I'm making)
hdpi.mat.b <- matrix(NA, nrow = 100, ncol = 2)
hdpi.mat.f <- matrix(NA, nrow = 100, ncol = 2)
# sample 1000 draws from the posterior because I don't want to loop through all 4000
sample.post <- sample_n(fits, size = 1000) %>% data.frame()
# for each draw of 1000, predict mean richness (if I wanted prediction interval I'd rpois it)
for (i in 1:1000){
  mus.f[,i] <- exp(sample.post[i,1]+sample.post[i,2]+(sample.post[i,3]+sample.post[i,4])*scale(lats))
  mus.b[,i] <- exp(sample.post[i,1]+sample.post[i,3]*scale(lats))
}
# now I have a vector giving predicted richness for each value of latitude. 1000 for each parameter input.
# take the HDPI of this, for each latitude
for(i in 1:100){
  hdpi.mat.b[i,] <- HPDI(mus.b[i,])
  hdpi.mat.f[i,] <- HPDI(mus.f[i,])
}
# now plot it 
# something is weird because I scaled it
{plot(richness~latitude, data = ants, pch = 16, col = as.factor(habitat), ylim = c(0,21), type = "n")
  for(i in 1:1000){
    params <- sample_n(fits,1) %>% data.frame()
    newd.f$Pred.r <- exp(as.matrix(newd.f[1:4]) %*% t(params))
    newd.b$Pred.r <- exp(as.matrix(newd.b[1:4]) %*% t(params))
    lines(y = newd.f$Pred.r, x = newd.f$latitude, lwd = .2, col = alpha("forestgreen",.1))
    lines(y = newd.b$Pred.r, x = newd.b$latitude, lwd = .2, col = alpha("brown",.1))}
  lines(y = hdpi.mat.b[,1], x=lats, lwd =2, lty =2, col = "saddlebrown")
  lines(y = hdpi.mat.b[,2], x=lats, lwd =2, lty =2, col = "saddlebrown")
  lines(y = newd.b$best, x=lats, lwd =2,  col = "saddlebrown")
  lines(y = hdpi.mat.f[,1], x=lats, lwd =2, lty =2, col = "darkgreen")
  lines(y = hdpi.mat.f[,2], x=lats, lwd =2, lty =2, col = "darkgreen")
  lines(y = newd.f$best, x=lats, lwd =2, col = "darkgreen")
  points(richness~latitude, pch=21, bg = as.factor(habitat), col = "black", ants)
}
# try with ggplot
hdpi.f <- hdpi.mat.f %>% data.frame() %>% mutate(lats = lats)
hdpi.b <- hdpi.mat.b %>% data.frame() %>% mutate(lats = lats)
mycolors = c("brown", "forestgreen")
ggplot()+
  geom_ribbon(mapping=aes(x=lats,ymin=X1,ymax=X2),
              alpha=0.2,fill = "forestgreen", data=hdpi.f)+
  geom_ribbon(mapping=aes(x=lats,ymin=X1,ymax=X2),
              alpha=0.2,fill = "brown", data=hdpi.b)+
  geom_line(data = newd.f, mapping = aes(x=latitude, y = best), color="forestgreen")+
  geom_line(data = newd.b, mapping = aes(x=latitude, y = best), color="brown")+
  geom_point(data = ants, aes(x=latitude, y = richness, col = habitat))+
  scale_color_manual(values =mycolors)+
  guides(color=guide_legend(title="Habitat Type"))+
  xlab("Latitude") + ylab("Species richness")
```

Scaling the latitude predictor makes this kind of plotting tough...not sure quite how to troubleshoot it.