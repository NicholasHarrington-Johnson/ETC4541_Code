#####################################################################
#
#  The program implements Algorithm 3 for the student-t data model, which is a
#  full Gibbs sampling algorithm for the augmented model that represents
#  the Student-t variables as an infinite mixture of normals, where the mixing is
#  according to an IG distribution. The model here uses the 'standard' 
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


# set up all functions needed for Algorithm 3

mu_fc.f = function(sig_in, w_in, y_in, psi_in){
  
  # this function generates from the full conditional of mu
  # see notes slide 52
  
  denom = sum(1/(w_in)^2)
  # (FIXED HERE)
  mean_mu = sum(y_in/w_in^2)/denom
  sd_mu = sig_in*sqrt(((psi_in-2)/psi_in)/denom)
  out = rnorm(1,mean=mean_mu,sd=sd_mu) #mu draw
  
  return(out)
  
}


sig_fc.f = function(mu_in, w_in, y_in, psi_in){
  
  # this function generates from the full conditional of mu
  # see notes slide 52
  
  T = length(y_in)
  sse = sum(((y_in - mu_in)/w_in)^2)
  sighatsq = psi_in*sse/(T*(psi_in - 2))
 
  out = 1/sqrt(rgamma(1,shape=(T/2),rate=(T*sighatsq/2)))
  
  return(out)
  
}


wblock_fc.f = function(mu_in, sig_in, y_in, psi_in){
  
  # this function generates from the "blocked" full conditional of 
  # the T-dimensional vector w
  # see notes slide 53
  
  T = length(y_in)
  
  # set rate vector (FIXED HERE)
  sighatsq_vec = psi_in*(1 + (((y_in - mu_in)/sig_in)^2)/(psi_in - 2))/(psi_in + 1)
  
  # note (FIXED HERE)
  v_param = psi_in+1
  out = 1/sqrt(rgamma(T,shape=(v_param/2),rate=(v_param*sighatsq_vec/2)))
  
  return(out)
  
}

### Algorithm 3 (Augmented Gibbs Sampler)

# define number of desired iterations
B = 450 # burn-in draws
NG = 2200 # retained draws

# fix psi to its true value
psi = psitrue

# set up storage for MCMC draws

AGibbs = matrix(c(0),nrow=(B+NG),ncol=nobs+2)

# initialise MCMC
AGibbs[1,] = c(mean(y),sd(y),rep(1,nobs))

#AGibbs[1,] = c(start_mean,start_sd,rep(1,nobs))

# run the augmented Gibbs scheme!
for(iter in 2:(B+NG)){
  # sample mu
  AGibbs[iter,1] = mu_fc.f(AGibbs[(iter-1),2], AGibbs[(iter-1),3:(nobs+2)], y, psi)
  
  # sample sigma
  AGibbs[iter,2] = sig_fc.f(AGibbs[iter,1], AGibbs[(iter-1),3:(nobs+2)], y, psi)

  AGibbs[iter,3:(nobs+2)] = wblock_fc.f(AGibbs[iter,1], AGibbs[iter,2], y, psi)
}

# Quick plot
quickgraph(AGibbs,B,3)

## GGPLOT

mu_draws <- AGibbs[,1]
sigma_draws <- AGibbs[,2]
draw.num <- c(1:nrow(AGibbs))

plot1 <- ggplot(data.frame(mu_draws,sigma_draws,draw.num),aes(x=draw.num,y=mu_draws))+geom_line(color="red")+labs(y=expression(paste(mu, " iterates")),x="Iterations")+ggtitle("Iterates from Algorithm 3 for Student t data problem")+geom_vline(xintercept=B,colour="green")
plot2 <- ggplot(data.frame(mu_draws,sigma_draws,draw.num),aes(x=draw.num,y=sigma_draws))+geom_line(color="blue")+labs(x="Iteratations",y=expression(paste(sigma, " iterates")))+geom_vline(xintercept=B,colour="green")

savepdf("Algorithm_3")
grid.arrange(plot1,plot2)
dev.off()

