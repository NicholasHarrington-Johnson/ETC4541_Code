################
##            ##
## Question 2 ##
##            ##   
## Part (d)   ##
##            ##
################

T=100

# Initialise values
x.1 <- rnorm(1,0,100)
x <- rep(0,T)

# True value of sigma.eta
sigma.eta <- 1

eta.t <- c(0,rnorm((T-1),mean=0,sd=(sigma.eta^2)))

x <- x.1+cumsum(eta.t)

# Number of observations
T <- length(x)

# Define initial values for nu_n and sigma.hat
# These are actually known so I am not going to change these
# Starting value for prior parameter nu_bar
nu_n <- 1
# Starting value for prior parameter 
# Just guess this
sigma.2.hat <- 1

# Function that describes the posterior of sigma

s.posterior <- function(x,nu_n,sigma.2.hat){
  # Figure out parameters to go into inverted gamma
  # Start with x vector
  s  <- x[2:T]-x[1:(T-1)]
  sq <- s^2
  ssq <- sum(sq)
  
  #  Define input parameters
  
  vbar <- T+nu_n-1
  
  sigma.bar <- (nu_n*sigma.2.hat+ssq)/vbar
    
  # Calculate shape and rate
  sigma <- 1/sqrt(rgamma(1,shape=vbar/2,rate=vbar*(sigma.bar)/2))
  return(sigma)
}

