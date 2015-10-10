#####################################################################
#
#  The program implements Algorithm 1: a Metropolis Hastings within Gibbs algorithm
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

# Plot a kernel density estimate
savepdf("Data_K_Est")
m <- ggplot(as.data.frame(y),aes(x=y))
m + geom_density(aes(colour="y",fill="y"),alpha=0.4) + ggtitle("Estimated Density of Data")+theme(legend.position="none")
dev.off()

# Print summary statistics as latex code
y.out <- lm(y~1)
summary(y.out)

# set up all functions needed for Algorithm 1

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

mu_MH.f = function(mu_last, sig_in, y_in, psi_in){
  
  # this function completes an MH step for the full conditional of mu
  # see notes slide 32
  
  # create data summaries needed for normal candidate distribution
  ybar = mean(y_in)
  sy = sd(y_in)
  T = length(y_in)
  
  mu_c = rnorm(1,mean=ybar,sd=(sy/sqrt(T))) #candidate draw
  
  # calculate the MH ratio, using logs for stability
  logalpha = logpkernel.f(y_in, mu_c, sig_in, psi_in) 
  logalpha = logalpha - logpkernel.f(y_in, mu_last, sig_in, psi_in)
  logalpha = logalpha - log(dnorm(mu_c, mean = ybar, sd = (sy/sqrt(T)))) 
  logalpha = logalpha + log(dnorm(mu_last, mean = ybar, sd = (sy/sqrt(T))))
  
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


sig_MH.f = function(mu_in, sig_last, y_in, psi_in){
  
  # this function completes an MH step for the full conditional of sig
  # see slide 32
  
  # create data summaries needed for IG candidate distribution
  # see, for example, slide 17
  
  T = length(y_in)
  sighatsq = sum((y_in - mu_in)^2)/T
  
  sig_shape = T/2
  sig_scale = T*sighatsq/2
    
  sig_c = 1/sqrt(rgamma(1,shape=sig_shape,rate = sig_scale))  #candidate draw
  
  # calculate the MH ratio, using logs for stability
  logalpha = logpkernel.f(y_in, mu_in, sig_c, psi_in) 
  logalpha = logalpha - logpkernel.f(y_in, mu_in, sig_last, psi_in)
  logalpha = logalpha - log(IG.f(sig_c, T, sighatsq))
  logalpha = logalpha + log(IG.f(sig_last, T, sighatsq))
  
  alpha = min(1, exp(logalpha))
  
  # choose between the candidate sig and the sig from the previous iteration
  if(alpha==1){
    out = sig_c
    
  } else if (runif(1,0,1)<=alpha){
    out = sig_c
    
  } else
    out = sig_last
    
  return(out)
  
}


### Algorithm 1 (MH using normal data full conditionals as candidate)

# define number of desired iterations
B = 100 # burn-in draws
NG = 1000 # retained draws

# fix psi to its true value
psi = psitrue

# set up storage for MCMC draws
MHGibbs = matrix(c(0),nrow=(B+NG),ncol=2)

# initialise MCMC
MHGibbs[1,] = c(mean(y),sd(y))

# run the MH within Gibbs scheme!
for(iter in 2:(B+NG)){
  # sample mu
  MHGibbs[iter,1] = mu_MH.f(MHGibbs[(iter-1),1], MHGibbs[(iter-1),2], y, psi)

  # sample sigma
  MHGibbs[iter,2] = sig_MH.f(MHGibbs[iter,1], MHGibbs[(iter-1),2], y, psi)
}

# Quick plot
quickgraph(MHGibbs,B,1)

## ggplot of MCMC
mu_draws <- MHGibbs[,1]
sigma_draws <- MHGibbs[,2]
draw.num <- c(1:nrow(MHGibbs))

plot1 <- ggplot(data.frame(mu_draws,sigma_draws,draw.num),aes(x=draw.num,y=mu_draws))+geom_line(color="red")+labs(y=expression(paste(mu, " iterates")),x="Iterations")+ggtitle("Iterates from Algorithm 1 for Student t data problem")+geom_vline(xintercept=B,colour="green")
plot2 <- ggplot(data.frame(mu_draws,sigma_draws,draw.num),aes(x=draw.num,y=sigma_draws))+geom_line(color="blue")+labs(x="Iteratations",y=expression(paste(sigma, " iterates")))+geom_vline(xintercept=B,colour="green")

savepdf("Algorithm_1")
grid.arrange(plot1,plot2)
dev.off()






