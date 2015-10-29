########################################################################
#
# Program to complete an MCMC for the local level model (states only)
# Generating state variables one at a time (i.e. with Gibbs sampling)
# conditional on the true variance parameters (simulated data)
#
# Code written by Catherine Forbes 2015
#
########################################################################

## set parameters
sigeta.true = 1
sigeps.true = 2

# Generate data values
set.seed(2339570)
T=100

y = rep(0,T) # initialise observed values
alpha.true = rep(0,T) # initialise true alpha values for periods 1:T

# Generate states at time t=0 and t=1, and observation at time t=1
alpha0.true = 25 # true alpha at period t=0
alpha.true[1] = rnorm(1, alpha0.true, sigeta.true) 
y[1] = rnorm(1,alpha.true[1],sigeps.true)

# Generate states and observations for times t=2 to t=T
for(i in 2:T){
  alpha.true[i] = rnorm(1, alpha.true[(i-1)], sigeta.true)
  y[i] = rnorm(1, alpha.true[i], sigeps.true)
}

# Plot the data
par(mfrow=c(1,1))
ry = range(y,alpha.true, alpha0.true) # ensure range is sufficient
plot(0:T,c(NA,y[1:T]),col=4,ylim=ry,xlab="t",ylab=expression(y[t]),pch=16)
title("Simulated local level data")
lines(0:T,c(alpha0.true,alpha.true),col="red")

# Now think about inference for the state variables, given the data and the variance parameters

# functions 

gen_alpha0.f=function(state1,sigsq0,sigsqstate){
  # this function generates alpha_0 from its full conditional distribution
  
  # condition on alpha_(t+1)
  condmean = state1*sigsq0/(sigsq0+sigsqstate)
#  condvar = sigsq0 - (sigsq0^2)/(sigsq0+sigsqstate) 
  condvar = sigsq0*sigsqstate/(sigsq0+sigsqstate) # more stable calculation
  
  out = rnorm(1,condmean,sqrt(condvar))
  return(out)
}

gen_alphat.f=function(stateminus1,stateplus1,obs_y,sigsqobs,sigsqstate){
  # this function generates alpha_t from its full conditional distribution
  
  # first condition on alpha_t on alpha_(t-1) and alpha(t+1)
  m = (stateminus1+stateplus1)/2 
  V = sigsqstate/2 
  
  # then condition on the observation too
  condmean = m + (obs_y - m)*V/(V+sigsqobs)
#  condvar = V - V^2/(V+sigsqobs)
  condvar = sigsqobs*V/(V+sigsqobs) # more stable calculation
  
  out = rnorm(1,condmean,sqrt(condvar))
  return(out)
}

gen_alphaT.f=function(stateminus1,data,sigsqobs,sigsqstate){
  # this function generates alpha_T from its full conditional distribution
  
  # first condition on alpha_t on alpha_(t-1)
  m = stateminus1
  V = sigsqstate
  
  # then condition on the observation too
  condmean = m + (data - m)*V/(V+sigsqobs)
# condvar = V - V^2/(V+sigsqobs)
 condvar = V*sigsqobs/(V+sigsqobs) # more stable calculation
  
  out = rnorm(1,condmean,sqrt(condvar))
  return(out)
}

# condition on the true sigmas
P0l0 = 10^2 # note a0l0=1 so has dropped out of the expressions
sigeta = sigeta.true
sigeps = sigeps.true

BGibbs = 100
MGibbs = 1000 
Mreps = BGibbs + MGibbs # total number of replications of the state vector

State0 = rep(0,Mreps) # storage for generated states at time t=0
States = matrix(c(0),Mreps,T) # storage generated states at times t=1,2,..,T

# arbitrary starting values
State0[1] = mean(y)
States[1,] = y

# now generate those states, one at a time, conditional on all other states 
# (and conditional on the data, and the other parameters)

for(iter in 2:Mreps){
  
  State0[iter]=gen_alpha0.f(States[(iter-1),1],P0l0,sigeta^2)
  
  for(t in 1:(T-1)){  
    if(t==1){
      States[iter,t]=gen_alphat.f(State0[iter],States[(iter-1),(t+1)],y[t],sigeps^2,sigeta^2)
    } else {
      States[iter,t]=gen_alphat.f(States[iter,(t-1)],States[(iter-1),(t+1)],y[t],sigeps^2,sigeta^2)
      }
  }
  States[iter,T]=gen_alphaT.f(States[iter,(T-1)],y[T],sigeps^2,sigeta^2)
  
  # Add generation of Gibbs draws for two variance parameters here
  
}

# Take a look at the generated states
plot(0:T,c(State0[1],States[1,]),col="grey",ylim=range(State0,States),xlab="t",ylab="generated states")
legend(x=70,y=25,legend=c("true states","observations","generated states"),cex=0.7,col=c("red","blue","grey"),lty=c(1,NA,NA),pch=c(NA,1,1))
title("States generated one at a time")

for(iter in 1:Mreps){
points(0:T,c(State0[iter],States[iter,]),col="grey",pch=1)
}
lines(0:T,c(alpha0.true,alpha.true),col="red")
points(1:T,y,col="blue")

library(coda)
my.mcmc = mcmc(States,start=BGibbs+1)
summary(my.mcmc)

effectiveSize(my.mcmc)
