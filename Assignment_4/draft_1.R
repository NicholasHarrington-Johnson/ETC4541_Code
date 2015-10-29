source("functions.R")
#############
## Question 2
## (b)

# Generate data values
set.seed(2339570)
T=100

# Initialise values
x.1 <- rnorm(1,0,100)
x <- rep(0,T)

# True value of sigma.eta
sigma.eta <- 1

eta.t <- c(0,rnorm((T-1),mean=0,sd=(sigma.eta^2)))

x <- x.1+cumsum(eta.t)

dat <- data.frame(seq_len(T),x)
colnames(dat) <- c("t","x")
p <- ggplot(dat,aes(x=t,y=x))+
  geom_line()+
  ggtitle(latex2exp("Simulated $x$ data"))
savepdf("Simulated_x_dat")
print(p)
dev.off()
dif.dat <- dat$x[2:nrow(dat)]-dat$x[1:(nrow(dat)-1)]
dif.dat <- data.frame(seq_len(T-1),dif.dat)
colnames(dif.dat)<-c("t","Difference")
d <- ggplot(dif.dat)+
  geom_density(aes(x=Difference),fill="red",alpha=0.4)+
  ggtitle(latex2exp("Density Plot Simulated $x$ data"))
savepdf("Density_x_dat")
print(d)
dev.off()

# Log likelihood function
# I've just assumed x as predefined
loglikelihood <- function(sigma){
  y <- (x[2:length(x)]-x[1:(length(x)-1)])
  y.2 <- y^2
  const.1 <- (100*2*pi)^-0.5
  const.2 <- 1/(sigma *( (2*pi)^0.5))
  l <- log(const.1)+(T-1)*log(const.2)-((x[1]^2)/200)-(sum(y.2)/(2*sigma^2))
  return(l)
}

potential.sigma <- seq(0.6,1.4,by=0.01)
log.l <- unlist(lapply(potential.sigma,loglikelihood))
# Print actual best estimate
print(potential.sigma[which(log.l==max(log.l))])

# Plot Log Likelihood
max.l <- potential.sigma[which(log.l==max(log.l))]

poi <- data.frame(sigmas =c(sigma.eta,max.l),sigma=c("True_Sigma","Estimated_Sigma"))

ldat <- data.frame(potential.sigma,log.l)

colnames(ldat) <- c("Sigma","Likelihood")
lp <- ggplot(ldat,aes(x=Sigma))+
  geom_line(aes(y=Likelihood))+
  geom_vline(aes(xintercept=sigmas,color=sigma),show_guide=TRUE,data=poi)+
  labs(x=latex2exp("$\\sigma"),y="Log Likelihood")+
  ggtitle(latex2exp("Log Likelihood for $\\sigma_{\\eta}$"))

savepdf("Log_Likelihood")
print(lp)
dev.off()

###


