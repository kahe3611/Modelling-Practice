hemo2 <- read.csv("~/Desktop/hemo_gen1.csv")
hemo2
set.seed(4.6) #make example reproducible
n <- 30 #size of dataset
b0 <- 50
b1 <- 10
s <- 40
x <- hemo2$Diet
y <- hemo2$Hemocytes
plot(x,y)

set.seed(4.6) #make example reproducible
n <- 30 #size of dataset
b0 <- 50
b1 <- 10
s <- 40
x <- runif(n,min=0,max=25)
y <- b0 + b1 * x + rnorm(n,sd=s)
plot(x,y)

lmod <- function(b0,b1,x){
  return( b0 + b1*x )
}

lmod(b0=300,b1=-9,x)

ysim <- function(mu,sigma) {
  return(rnorm(n=length(mu),mean=mu,sd=sigma))
}

par(mfrow=c(2,2),mar=c(5, 4, 0, 2) + 0.1)
ylim <- c(0,400)
for ( run in 1:4 ) {
  yout <- ysim( mu=lmod(b0=300,b1=-9,x), sigma=30 )
  plot(x,yout,ylim=ylim,ylab="y")
}

lm_nll <- function(p,y,x) {
  mu <- lmod(b0=p[1],b1=p[2],x) #call the linear model
  nll <- -sum(dnorm(y,mean=mu,sd=p[3],log=TRUE)) #-1 * sum of log-likelihoods 
  return(nll)
}

p <- c(70,8,30)
lm_nll(p,y,x)

fitlm <- optim(p=c(50,200/25,75/2),lm_nll,y=y,x=x)
fitlm

lm_slopefixed_nll <- function(p,b1,y,x) {
  mu <- lmod(b0=p[1],b1,x) #call the linear model
  nll <- -sum(dnorm(y,mean=mu,sd=p[2],log=TRUE)) #-1 * sum of log-likelihoods 
  return(nll)
}

fitlm_slopenull <- optim(p=c(200,75),lm_slopefixed_nll,b1=0,y=y,x=x)
fitlm_slopenull

exp(-fitlm$value)/exp(-fitlm_slopenull$value)

nll_bl <- rep(NA, 50)
b1_range <- seq(6,12, length.out = (nll_b1))
par <- c(74,35)
i <- 1
for(b1 in b1_range){
  nll_b1[1] <- optim(p = par, lm_slopefixed_nll, b1 = b1, y = y, x = x)$value
  i <- i + 1
}

likprof_b1 <- exp(-nll_bl)
likeratio_b1 <- exp(fitlm$value-nll_b1)

par(mfrow = c(1,3))
plot(b1_range, nll_b1)