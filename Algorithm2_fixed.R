#####################################################################
#
#  The program implements a Random Walk Metropolis Hastings algorithm
#  for the mean of a student-t data model with the 'standard' 
#  noninformative prior proporitonal to 1/sigma
#
#  Code written by Catherine Forbes (let me know if you spot any errors!)
#
#######################################################################

# Generate some iid data y = mu + sig*eps where eps is standardised student-t(5)
set.seed(23395710)

mutrue <- 5
sigtrue <- 1
psitrue <- 5
lam<-sqrt((psitrue-2)/psitrue)
nobs<-50
eps<-rt(n=nobs,df=psitrue)*lam
y<-mutrue+sigtrue*eps

# take a look at the data just generated
par(mfrow=c(1,1))
plot(y,type="l")
plot(density(y))


# set up all functions needed for Algorithm 2

logpkernel.f = function(y_in, mu_in, sig_in, psi_in){
  
  # this function produces the log of equation (4) on slide 32
  
  T = length(y_in) 
  logkernel = -(T+1)*log(sig_in) # a scalar
  terms = -0.5*(psi_in+1)*log(1+(((y_in - mu_in)/sig_in)^2)/(psi_in-2))   # a vector
  
 out = logkernel+sum(terms)   # a scalar
 return(out)
}

IG.f = function(sigin, vin, sigsqhatin){
  
  # this is essential the "Gael_dens_IG.f" function
  # it computes the value of the pdf of the IG(vin,sigsqhatin) distribution
  # at each value of the variate 'sigin'
  
  out = (2/gamma(vin/2))*((vin*sigsqhatin/2)^(vin/2))*exp(-vin*sigsqhatin/(2*sigin^2))/(sigin^(vin+1))
  return(out)
}

mu_RWMH.f = function(mu_last, sig_in, y_in, psi_in, sig_tune){
  
  # this function completes an RWMH step for the full conditional of mu
  # see notes slide 43
  
  mu_c = rnorm(1,mean=mu_last,sd=sig_tune) #candidate draw
  
  # calculate the MH ratio, using logs for stability
  # note do not need to evaluate the candidate density at all!
  logalpha = logpkernel.f(y_in, mu_c, sig_in, psi_in) 
  logalpha = logalpha - logpkernel.f(y_in, mu_last, sig_in, psi_in)
  
  alpha = min(1, exp(logalpha))
   
  # choose between the candidate mu and the mu from the previous iteration
  if(alpha==1){
    out = mu_c
  } else if (runif(1,0,1)<=alpha){
    out = mu_c
  } else
    out = mu_last
  
return(out)

}


logpzetakernel.f = function(y_in, mu_in, zeta_in, psi_in){
  
  # this function produces the log of target p*(zeta)
  # see slide 44
  
  T = length(y_in) 
  logkernel = -zeta_in*(T+1) # a scalar
  terms = -0.5*(psi_in+1)*log(1+(((y_in - mu_in)/exp(zeta_in))^2)/(psi_in-2))   # a vector
  
  out = logkernel+sum(terms)   # a scalar
  return(out)
}


sig_RWMH.f = function(mu_in, sig_last, y_in, psi_in, sig_tune){
  
  # this function completes an RWMH step for the full conditional of sig
  # see slide 44
  
  # transform sig_last so no longer constrained
  zeta_last = log(sig_last)
  
  zeta_c = rnorm(1,mean=zeta_last,sd=sig_tune) #candidate draw
  
  # calculate the MH ratio, using logs for stability
  
  # note do not need to evaluate the candidate density at all!
  logalpha = logpzetakernel.f(y_in, mu_in, zeta_c, psi_in) 
  logalpha = logalpha - logpzetakernel.f(y_in, mu_in, zeta_last, psi_in)
  
  alpha = min(1, exp(logalpha))
  
  # choose between the candidate sig and the sig from the previous iteration
  if(alpha==1){
    out = zeta_c
    
  } else if (runif(1,0,1)<=alpha){
    out = zeta_c
    
  } else
    out = zeta_last
    
  #transform back to return sigma, not zeta
  out = exp(out) 
  return(out)
  
}


### Algorithm 2 (Random Walk MH)

# define number of desired iterations
B = 100 # burn-in draws
NG = 1000 # retained draws

# fix psi to its true value
psi = psitrue

# set up storage for MCMC draws
RWMHGibbs = matrix(c(0),nrow=(B+NG),ncol=2)

# initialise MCMC
RWMHGibbs[1,] = c(mean(y),sd(y))

# PLAY AROUND WITH THESE VALUES
tuner<-c(0.1,0.1) #tuning parameters for mu and sig, respectively

# run the MH within Gibbs scheme!
for(iter in 2:(B+NG)){
  # sample mu (FIXED HERE)
  RWMHGibbs[iter,1] = mu_RWMH.f(RWMHGibbs[(iter-1),1], RWMHGibbs[(iter-1),2], y, psi,tuner[1])
  
  # sample sigma (FIXED HERE)
  RWMHGibbs[iter,2] = sig_RWMH.f(RWMHGibbs[iter,1], RWMHGibbs[(iter-1),2], y, psi,tuner[2])
}
  
# Quick plot
quickgraph(RWMHGibbs,B,2)


## GGPLOT
mu_draws <- RWMHGibbs[,1]
sigma_draws <- RWMHGibbs[,2]
draw.num <- c(1:nrow(MHGibbs))

plot1 <- ggplot(data.frame(mu_draws,sigma_draws,draw.num),aes(x=draw.num,y=mu_draws))+geom_line(color="red")+labs(y=expression(paste(mu, " iterates")),x="Iterations")+ggtitle("Iterates from Algorithm 2 for Student t data problem")+geom_vline(xintercept=B,colour="green")
plot2 <- ggplot(data.frame(mu_draws,sigma_draws,draw.num),aes(x=draw.num,y=sigma_draws))+geom_line(color="blue")+labs(x="Iteratations",y=expression(paste(sigma, " iterates")))+geom_vline(xintercept=B,colour="green")

savepdf("Algorithm_2")
grid.arrange(plot1,plot2)
dev.off()