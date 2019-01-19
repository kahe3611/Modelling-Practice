# homework from chapter 3 in McElreath

# Homework problem 3M1: w = water and n = #of tosses

w <- 8

n <- 15

grid_approx <- seq(from = 0, to = 1, length.out = 1000)

prior <- rep(x = 1, length(grid_approx))

p_likelihood <- dbinom(x = w, size = n, prob = grid_approx)

unstand.posterior <- p_likelihood * prior

posterior <- unstand.posterior / sum(unstand.posterior)

plot(posterior ~ grid_approx, type = "l")




# Homework 3M2: need to draw 10000 samples from previous code

trials <- 1e4

samples <- sample(x = grid_approx, size = trials, prob = posterior, replace = TRUE)

#installed HPDinterval package, otherwise ran into errors; can't get hdi to work probably
#library(HDInterval)
#hdi(samples = samples, prob = .9)

??HPDI


# Homework 3M3
# got this code to work
posterior_distribution <- rbinom(n = trials, size = n, prob = samples)

mean(posterior_distribution == 8)




# Homework 3M4

predictive_dist <- rbinom(n = trials, size = 9, prob = samples)

mean(predictive_dist == 6)



# I skipped Homework 3M5


## 3H1: Tried loading "rethinking" and "homeworkch3" but doesn't work for me.

install.packages(c("devtools","mvtnorm","loo","coda"), repos="https://cloud.r-project.org/",dependencies=TRUE) 
library(devtools)
install_github("rmcelreath/rethinking")

library(rethinking)

data(homeworkch3)




total_births <- length(birth1) + length(birth2)

boys_born <- sum(birth1 + birth2)

girls_born <- total_births - boys_born




p_grid <- seq(from = 0, to = 1, length.out = 1000)

prior <- rep(x = 1, length(p_grid))

likelihood <- dbinom(x = boys_born, size = total_births, prob = p_grid)

unstandard_posterior <- likelihood * prior

posterior <- unstandard_posterior / sum(unstandard_posterior)

plot(posterior ~ p_grid, type = "l")




p_grid[which.max(posterior)]




# homework 3H2

trials <- 1e4

samples <- sample(x = p_grid, size = trials, prob = posterior, replace = TRUE)

HPDI(samples = samples, prob = c(.5, .89, .97))




# homework 3H3

n <- total_births

posterior_predictive <- rbinom(n = trials, size = n, prob = samples)

dens(posterior_predictive)

abline(v = boys.born, col = "red")




# homework 3H4

n <- 100

sum(birth1)

posterior_pred <- rbinom(n = trials, size = n, prob = samples)

dens(posterior_pred, adj = .1)

abline(v = sum(birth1), col = "red" )




# homework 3H5

boysborn_aftergirls <- birth2[birth1 == 0]
boysborn_aftergirls
posterior_pred2 <- rbinom(n = trials, size = length(boysborn_aftergirls), prob = samples)
posterior_pred2
dens(posterior_pred2)

abline(v = sum(boysborn_aftergirls), col = "red")
