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

# Create MCMC objects
B1 = 100 # burnin may be redefined here
M1 = 1000 # number of retained draws redefined here too
B2 <- 100
M2 <- 1000
B3 <- 100
M3 <- 1000

MH.mcmc = mcmc(MHGibbs[(B1+1):(B1+M1),])
RWMH.mcmc = mcmc(RWMHGibbs[(B2+1):(B2+M2),])
# For A3 just ignore ws
A3.mcmc = mcmc(AGibbs[(B3+1):(B3+M3),1:2])

# Numerical distributional summaries
summary(MH.mcmc)
summary(RWMH.mcmc)
summary(A3.mcmc)
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
M1/effectiveSize(MH.mcmc)
M2/effectiveSize(RWMH.mcmc)
M3/effectiveSize(A3.mcmc)

# Rejection rates for each variable in each Markov chain
rejectionRate(MH.mcmc)
rejectionRate(RWMH.mcmc)
rejectionRate(A3.mcmc)