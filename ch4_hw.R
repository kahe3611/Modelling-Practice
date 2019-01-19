# Chapter 4 Homework 4H1 - 4H3

#4H1 
library(rethinking)

sampost = function(flist, data, n=10000){
  quadapprox = map(flist, data)
  posterior_sample = extract.samples(quadapprox,n)
  return(posterior_sample)
}

data(Howell1)
d <- Howell1

d2 <- d[ d$age >= 18 , ]
min(d2$height)
max(d2$height)
mean(d2$height)
weights <-c(46.95,43.72,64.78,32.59,54.63)
d2$weight.center <- d2$weight - mean(d2$weight) 
?sampost
#fit a model
# couldn't get sampost function to work, not sure what to try here
#nevermind, I missed the code at the beginning of ch 4 to create "sampost"
model <- sampost(
  alist(
    height ~ dnorm(mu, sigma) ,
    mu <- a + b*weight.center,
    a ~ dnorm(154, 100) , #average height
    b ~ dnorm(0, 10) ,
    sigma ~ dunif(0, 50)
  ),
  data=d2 
  )

precis(model)

#can't get this to work because the model isn't working

weights.center <- weights - mean(d2$weight)
#post<- extract.samples(model)
sampost

sim.height <- sapply(weights.center, function(weight) {
  rnorm(
    n = nrow(post),
    mean = post$a + post$b*weight,
    sd = post$sigma
  )
})

height.PI <- apply(sim.height, 2, PI, prob=0.89)
height.HPDI <- apply(sim.height, 2, HPDI, prob=0.89)
height.HPDI
height.mean <- apply(sim.height, 2, mean)
height.mean
pred_df <- data.frame("Individual"=1:5, "Weight"=weights, "Exptected_height"=height.mean, 
                      "PI_89_lower"=height.PI[1,], "PI_89_upper"=height.PI[2,])
pred_df

#4H2: fitting a linear regression 

d18 <- d[ d$age < 18, ]
d18$weight.c <- d18$weight - mean(d18$weight)   # centering the data

# fit the model

model_18 <- map(
  alist(
    height ~ dnorm( mu, sigma) ,
    mu <- a + b*weight.c ,
    a ~ dnorm( 156, 100) ,
    b ~ dnorm( 0, 10) ,
    sigma ~ dunif(0, 50)
  ),
  data=d18
)
precis(model_18)

weight.seq <- seq(from=-15, to=30, length.out = 30)             
#get posterior sample
post <- extract.samples(model_18) 

# computing mu
mu.link <- function(weight.c) post$a + post$b*weight.c           
mu <- sapply(weight.seq, mu.link)

mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob=0.89)

# generating some predicted heigh with this model
sim.height <- sapply( weight.seq, function(weight) {
  rnorm(
    n = nrow(post),
    mean = post$a + post$b*weight,
    sd = post$sigma
  )
})

#heigh hpdi and mean
height.HPDI <- apply(sim.height, 2, HPDI, prob=0.89)
height.mean <- apply(sim.height, 2, mean)

# plotting the data and putting a regression line
plot(height ~ weight.c, data=d18, col=col.alpha(rangi2, 0.9), ylim=c(50, 180))   
lines(weight.seq, mu.mean)                                      
shade( mu.HPDI, weight.seq)                                     
shade( height.HPDI, weight.seq)                                 

#4H3: use entire data set with the model provided

d <- Howell1
# fit the model

model.l <- map(
  alist(
    height ~ dnorm( mu, sigma) ,
    mu <- a + b*log(weight) ,
    a ~ dnorm( 178, 100) ,
    b ~ dnorm( 0, 100) ,                
    sigma ~ dunif(0, 50)
  ),
  data=d
)
precis(model.l)

weight.seq <- seq(from=2, to=65, length.out = 70) 
min(d$weight)
max(d$weight)
# min(d$weight) = 4.25, max(d$weight) = 62.99

post <- extract.samples(model.l)                             
# generate mu
mu.link <- function(weight) post$a + post$b*log(weight)       
mu <- sapply(weight.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob=0.89)

#model height
sim.height <- sapply( weight.seq, function(weight) {
  rnorm(
    n = nrow(post),
    mean = post$a + post$b*log(weight),
    sd = post$sigma
  )
})

height.HPDI <- apply(sim.height, 2, HPDI, prob=0.89)
height.mean <- apply(sim.height, 2, mean)

# plotting data and getting map regression line
#trying to draw hpdi area around the regression line
plot(height ~ log(weight), data=d, col=col.alpha(rangi2, 0.6))
lines(log(weight.seq), mu.mean) 

#this did work and added the shading
shade( height.HPDI, log(weight.seq))                                 


#plotting without log scale and added the shading around the regression line
plot(height ~ weight, data=d, col=col.alpha(rangi2, 0.6))
lines(weight.seq, mu.mean)

shade( height.HPDI, weight.seq)
