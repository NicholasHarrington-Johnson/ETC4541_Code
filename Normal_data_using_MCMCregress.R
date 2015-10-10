#####################################################################
#
#  The program implements a Random Walk Metropolis Hastings algorithm
#  for a normal data regression problem, using the MCMCpack function "MCMCregress"
#
#  Code written by Catherine Forbes (let me know if you spot any errors!)
#
#######################################################################
library(MCMCpack)

#rm(list=ls())

# Generate some iid data y = mu + sig*eps where eps~iid N(0,1)
set.seed(2015)

betatrue = c(21,2.5)
sigtrue = 6

nobs = 100
x = runif(nobs,min=5,max=15)
y = betatrue[1]+x*betatrue[2]+sigtrue*rnorm(nobs,0,1)

xy.df = data.frame(x,y)

par(mfrow=c(1,1))
plot(x,y,ylim=range(y-1,y+1),xlim=range(x+1,x-1),col="blue")
title("Scatterplot of normal regression data")
abline(lsfit(x,y),col="blue",lwd=2)


classical = lm(y~x)
summary(classical)
plot(classical)

bayesian = MCMCregress(y~x, data=xy.df)
summary(bayesian)
plot(bayesian)
