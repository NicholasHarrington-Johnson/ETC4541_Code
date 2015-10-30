#####################################################################
#
#  Use the "coda" package to get some MCMC summaries and diagnostics
#  See the "coda.pdf" file on Moodle for details of all functions used here
#  
#  This file uses input from Algorithm 1 (MHGibbs) and Algorithm 2 (RWMHGibbs)
#  associated with the Student-t data problem
#
#  Code written by Catherine Forbes  
#
#######################################################################

# Add Algorithm 3

# Be sure you have installed the coda package before you start
library(coda)
library(ggmcmc)
# Create MCMC objects
B1 = 700 # burnin may be redefined here
M1 = 2000 # number of retained draws redefined here too
B2 <- 600
M2 <- 2500
B3 <- 450
M3 <- 2200
## Reedefine for part c
B1 = 700 # burnin may be redefined here
#M1 = 10000 # number of retained draws redefined here too
B2 <- 600
#M2 <- 10000
B3 <- 450
#M3 <- 10000

MH.mcmc = mcmc(MHGibbs[(B1+1):(B1+M1),])
RWMH.mcmc = mcmc(RWMHGibbs[(B2+1):(B2+M2),])
# For A3 just ignore ws
A3.mcmc = mcmc(AGibbs[(B3+1):(B3+M3),c(1:2)])

# Numerical distributional summaries
summary(MH.mcmc)
summary(RWMH.mcmc)
summary(A3.mcmc)

####################################################################################

# GGPLOT stuff

# A1
A_1 <- ggs(MH.mcmc)
ggmcmc(A_1,file="Alg_1_Density.pdf",plot="density")
ggmcmc(A_1,file="Alg_1_running.pdf",plot="running")
ggmcmc(A_1,file="Alg_1_ACF.pdf",plot="autocorrelation")

# A2
A_2 <- ggs(RWMH.mcmc)
ggmcmc(A_2,file="Alg_2_Density.pdf",plot="density")
ggmcmc(A_2,file="Alg_2_running.pdf",plot="running")
ggmcmc(A_2,file="Alg_2_ACF.pdf",plot="autocorrelation")

# A3
A_3 <- ggs(A3.mcmc)
ggmcmc(A_3,file="Alg_3_Density.pdf",plot="density")
ggmcmc(A_3,file="Alg_3_running.pdf",plot="running")
ggmcmc(A_3,file="Alg_3_ACF.pdf",plot="autocorrelation")


####################################################################################


# Traceplots
traceplot(MH.mcmc,main="Traceplot for Algorithm 1")
traceplot(RWMH.mcmc,main="Traceplot for Algorithm 2")
traceplot(A3.mcmc,main="Traceplot for Algorithm 3")
# Put all traceplots on one page
par(mfrow=c(3,2))
traceplot(MH.mcmc,main="Traceplot for Algorithm 1")
traceplot(RWMH.mcmc,main="Traceplot for Algorithm 2")
traceplot(A3.mcmc,main="Traceplot for Algorithm 3")

# Autocorrelation in MCMC output
# numerical summaries

autocorr(MH.mcmc,lags=c(seq(1,10,1),seq(20,50,10)))
autocorr(RWMH.mcmc,lags=c(seq(1,10,1),seq(20,50,10)))
autocorr(A3.mcmc,lags=c(seq(1,10,1),seq(20,50,10)))

# compare with...
autocorr(MH.mcmc)
autocorr(RWMH.mcmc)
autocorr(A3.mcmc)

autocorr.diag(MH.mcmc,lags=c(seq(1,10,1),seq(20,50,10)))
autocorr.diag(RWMH.mcmc,lags=c(seq(1,10,1),seq(20,50,10)))
autocorr.diag(A3.mcmc,lags=c(seq(1,10,1),seq(20,50,10)))

# ACF plots
autocorr.plot(MH.mcmc,main="ACF for Algorithm 1")
autocorr.plot(RWMH.mcmc,main="ACF for Algorithm 2")
autocorr.plot(A3.mcmc,main="ACF for Algorithm 3")
plog <- ggplot(as.data.frame(MHGibbs[,1]),aes(x=c(seq(-5,10,1))))
# Density plots 
par(mfrow=c(3,2))
densplot(MH.mcmc,main="Density from Algorithm 1")
densplot(RWMH.mcmc,main="Density from Algorithm 2")
densplot(A3.mcmc,main="Density from Algorithm 3")

# Cumuplots (Cumulative quantile plots)
# I have modified the default "auto.layout=TRUE" to fit all on one page
cumuplot(MH.mcmc,main="Cumuplot from Algorithm 1",auto.layout=FALSE)
cumuplot(RWMH.mcmc,main="Cumuplot from Algorithm 2",auto.layout=FALSE)
cumuplot(A3.mcmc,main="Cumuplot from Algorithm 3",auto.layout=FALSE)

# HPD intervals
HPDinterval(MH.mcmc)
HPDinterval(RWMH.mcmc)
HPDinterval(A3.mcmc)

# Effective sample size for each variable in each Markov chain
effectiveSize(MH.mcmc)
effectiveSize(RWMH.mcmc)
effectiveSize(A3.mcmc)

# Inefficiency Factors for each variable in each Markov chain
M <-10000
M/effectiveSize(MH.mcmc)
M/effectiveSize(RWMH.mcmc)
M/effectiveSize(A3.mcmc)

# Rejection rates for each variable in each Markov chain
rejectionRate(MH.mcmc)
rejectionRate(RWMH.mcmc)
rejectionRate(A3.mcmc)

# Gelman
#gelman.diag(MH.mcmc,confidence=0.95,transform=FALSE,autoburnin = FALSE,multivariate = TRUE)

################################################################
#Do an MCMC object for theta
MH_theta <- data.frame(AGibbs[(B3+1):(B3+M3),1]/AGibbs[(B3+1):(B3+M3),2])
colnames(MH_theta) <- c("theta")
MH_theta.mcmc <- mcmc(MH_theta)
A_T <- ggs(MH_theta.mcmc)
# I don't actually like the way this density plot comes out
#ggmcmc(A_T,file="theta_density.pdf",plot="density")
m <- ggplot(MH_theta,aes(x=theta))+geom_density(aes(fill="theta",colour="theta"),alpha=0.4)+
  ggtitle("Estimated Density of theta")+theme(legend.position="none")
savepdf("theta_density")
print(m)
dev.off()
summary(MH_theta.mcmc)

################################################################

# Alt stuff

A3_alt.mcmc <- mcmc(alt_AGibbs[(B3+1):(B3+M3),c(1:2)])
S_A3_alt <- ggs(A3_alt.mcmc)

S_2_chains <- mcmc.list(A3.mcmc,A3_alt.mcmc)
vis_2_chain <- ggs(S_2_chains)
ggmcmc(vis_2_chain,file="2_chain_running.pdf",plot="running")
ggmcmc(vis_2_chain,file="2_chain_acf.pdf",plot="autocorrelation")
ggmcmc(vis_2_chain,file="2_chain_density.pdf",plot="density")
ggmcmc(vis_2_chain,file="2_chain_tp.pdf",plot="traceplot")
