# Set priors for my functions
nu_eps <- 1
nu_eta <- 1
sigma.2.hat_eps <- 1
sigma.2.hat_eta <- 1
## set parameters
sigeta.true = 1
sigeps.true = 2

# Generate data values
set.seed(23395710)
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

# condition on the true values of the parameters

# Now think about inference for the state variables, given the data and the variance parameters
# First run FF (only once) and save all output

# apred represents the one step ahead predicted state mean, a_t|(t-1)
# Ppred represents the one step ahead predicted state variance, P_t|(t-1)
# afilt represents the filtered state mean, a_t|t
# Pfilt represents the filtered state variance, P_t|t
# v represents the one step ahead prediction error, or innovation
# M represents the Kalman gain 
# F represents the variance of the innovation
# sigsqeps is the variance of the measurement error (eps)
# sigsqeta is the variance of the state error (eta)

# Forward filter
FFfunction <- function(sigeta,sigeps){
  # initialise variables
  apred = rep(0,T)
  Ppred = rep(0,T)
  afilt = rep(0,T)
  Pfilt = rep(0,T)
  v = rep(0,T)
  M = rep(0,T)
  F = rep(0,T)
  
  # setting initial state distribution
  a0l0 = 0
  P0l0 = 100
for(i in 1:T){
  if(i==1){
    apred[i] = a0l0 
    Ppred[i] = P0l0 + sigeta^2
  } else {
    apred[i] = afilt[(i-1)]
    Ppred[i] = Pfilt[(i-1)] + sigeta^2
  }
  
  v[i] = y[i]-apred[i]
  F[i] = Ppred[i] + sigeps^2
  M[i] = Ppred[i]
  
  afilt[i] = apred[i]+M[i]*v[i]/F[i]
  Pfilt[i] = Ppred[i]-M[i]^2/F[i]
}
  aP <- data.frame(afilt,Pfilt)
  return(aP)
}

init <- FFfunction(sigeta.true,sigeps.true)

lines(1:T,afilt[1:T],col="green")

# Backward sampling - want to replicate this many times, so put into a function
BSt.f = function(filtat,filtPt,nextstatestar,sigsqstate){
  
  vtstar = nextstatestar - filtat
  Ftstar = filtPt + sigsqstate
  Mtstar = filtPt

  asmt = filtat + Mtstar*vtstar/Ftstar
  Psmt = sqrt(filtPt - Mtstar^2/Ftstar)
  if(Psmt>1e9){
    print("Psmt is too big")
    Psmt <- 1e9
  }
  if(is.nan(Psmt)){
    print("Psmt is negative square root")
    Psmt <- abs(filtPt - Mtstar^2/Ftstar)
    Psmt <- sqrt(Psmt)
  }
  out = rnorm(1,mean = asmt,sd=Psmt)

  return(out)
}
# Draw iteration 1 of theta from its full conditional based on the first forward filter
s.eps.1 <- s.sigma.eps(T,nu_eps,y,afilt,sigma.2.hat_eps)
s.eta.1 <- s.sigma.eta(T,nu_eta,alpha.t=afilt,alpha.tp=c(a0l0,afilt[1:(T-1)]),sigma.2.hat_eta)

Bl_BGibbs = 100
Bl_MGibbs = 10000
Bl_Mreps = Bl_BGibbs + Bl_MGibbs # total number of replications of the state vector

BState0 = rep(0,Bl_Mreps) # storage for generated states for time t=0
BStates = matrix(c(0),Bl_Mreps,T) # storage generated states for times t=1,2,..,T

# same arbitrary starting values
BState0[1] = mean(y)
BStates[1,] = y

eta.vec <- rep(s.eta.1,Bl_Mreps)
eps.vec <- rep(s.eps.1,Bl_Mreps)
for(iter in 1:Bl_Mreps){
  
  BStates[iter,T] = rnorm(1,mean=init[T,1],sqrt(init[T,2]))
  for(t in (T-1):0){
    if(t==0){
      sigeta <- 1
      BState0[iter]=BSt.f(a0l0,P0l0,BStates[iter,1],sigeta^2)
    } else {
      alph <- FFfunction(eta.vec[iter],eps.vec[iter]) 
      BStates[iter,t]=BSt.f(alph[t,1],alph[t,2],BStates[iter,(t+1)],eta.vec[iter]^2)
      eta.vec[iter+1] <- s.sigma.eta(T,nu_eta,alpha.t=alph[,1],alpha.tp=c(a0l0,alph[1:(T-1),1]),sigma.2.hat_eta)
      eps.vec[iter+1] <- s.sigma.eps(T,nu_eps,y=y,alpha.t=alph[,1],sigma.2.hat_eps)
    }
  }
}

plot(0:T,c(BState0[1],BStates[1,]),col="darkgrey",ylim=range(BState0,BStates))
legend(x=70,y=25,legend=c("true states","observations","generated states"),cex=0.7,col=c("red","blue","darkgrey"),lty=c(1,NA,NA),pch=c(NA,1,1))
title("States generated in one block")


for(iter in 1:Bl_Mreps){
  points(0:T,c(BState0[iter],BStates[iter,]),col="darkgrey",pch=1)
}
lines(0:T,c(alpha0.true,alpha.true),col="red")
points(1:T,y,col="blue")


library(coda)
my.Bl_mcmc = mcmc(BStates,start=(Bl_BGibbs+1))
summary(my.Bl_mcmc)

effectiveSize(my.Bl_mcmc)


