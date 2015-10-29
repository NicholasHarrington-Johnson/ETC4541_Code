################
##            ##
## Question 3 ##
##            ##   
## Part (f)   ##
##            ##
################

## Sampling sigma eta

s.sigma_eta <- function(T,nu_eta,alpha.vec,sigma.2.hat_eta){
  # Defining relevant input parameters
  
  nu_bar <- T+nu_eta
  # alpha vector summed
  alp2 <- sum((alpha.vec[2:(T+1)]-alpha.vec[1:T])^2)
  sig.bar <- (nu_bar*sigma.2.hat_eta+alp2)/(nu_bar)
  
  # Sampling
  
  sigma <- 1/sqrt(rgamma(1,shape=nu_bar/2,rate=nu_bar*(sig.bar^2)/2))
  
  return(sigma)
  
}

## Sampling sigma epsilon

s.sigma.eps <- function(T,nu_eps,y,alpha.vec,sigma.2.hat_eps){
  # Defining relevant input parameters
  
  nu_bar <- T+nu_eps
  # y and alpha vector
  # Note that alpha vector element 0 is ignored
  # Alpha vector is of length (T+1)
  
  y2 <- sum((y-alpha.vec[2:length(alpha.vec)])^2)
  
  sig.bar <- (nu_bar*sigma.2.hat_eps+y2)/(nu_bar)
  
  # Sampling
  
  sigma <- 1/sqrt(rgamma(1,shape=nu_bar/2,rate=nu_bar*(sig.bar^2)/2))

  return(sigma)

}