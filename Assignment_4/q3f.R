################
##            ##
## Question 3 ##
##            ##   
## Part (f)   ##
##            ##
################

## Sampling sigma eta

s.sigma.eta <- function(T,nu_eta,alpha.t,alpha.tp,sigma.2.hat_eta){
  # Defining relevant input parameters
  
  nu_bar <- T+nu_eta
  # alpha vector summed
  alp2 <- sum((alpha.t-alpha.tp)^2)
  sig.bar <- (nu_bar*sigma.2.hat_eta+alp2)/(nu_bar)
  
  # Sampling
  
  sigma <- 1/sqrt(rgamma(1,shape=nu_bar/2,rate=nu_bar*(sig.bar)/2))
  
  return(sigma)
  
}

## Sampling sigma epsilon

s.sigma.eps <- function(T,nu_eps,y,alpha.t,sigma.2.hat_eps){
  # Defining relevant input parameters
  
  nu_bar <- T+nu_eps
  # y and alpha vector
  y2 <- sum((y-alpha.t)^2)
  
  sig.bar <- (nu_bar*sigma.2.hat_eps+y2)/(nu_bar)
  
  # Sampling
  
  sigma <- 1/sqrt(rgamma(1,shape=nu_bar/2,rate=nu_bar*(sig.bar)/2))

  return(sigma)

}