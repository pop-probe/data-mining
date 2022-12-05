library(nimble)
library(coda)

#1-(a)
n=1000
p=4

Xmat <- rep(0,n*p)
dim(Xmat) <- c(n,p)

for (i in 1:n){
  for (j in 1:p) {
    Xmat[i,j] <- rnorm(1,0,1)
  }
}

head(Xmat)
dim(Xmat)
class(Xmat)

#1-(b)
trueBeta = c(0.5,-0.5,0,1)
dim(trueBeta) = c(4,1)

Mux = function(Xi,beta){
  return(exp(Xi%*%beta))
}

Y <- rep(0,n*1)
dim(Y) <- c(n,1)



for (i in 1:n) {
  Y[i] = rpois(1,Mux(Xmat[i,],trueBeta))
}

head(Y)

#2 MCMC in hand
#setting
niter = 10000
beta_count = 0
prior_sd = sqrt(10)
prop_sd = 0.02

#Initialize beta
beta <- rep(0,niter*p)
dim(beta) <- c(niter,p)

head(beta)

#mcmc sampling
for (i in 2:niter) {
  
  #updating beta
  beta0_star = rnorm(1,beta[i-1,1],prop_sd)
  beta1_star = rnorm(1,beta[i-1,2],prop_sd)
  beta2_star = rnorm(1,beta[i-1,3],prop_sd)
  beta3_star = rnorm(1,beta[i-1,4],prop_sd)
  
  beta_star = c(beta0_star, beta1_star, beta2_star, beta3_star)
  dim(beta_star) = c(4,1)

  dnorm_star = log(dnorm(beta0_star,0,prior_sd)) + log(dnorm(beta1_star,0,prior_sd)) + log(dnorm(beta2_star,0,prior_sd)) + log(dnorm(beta3_star,0,prior_sd))
  dnorm_pre = log(dnorm(beta[i-1,1],0,prior_sd)) + log(dnorm(beta[i-1,2],0,prior_sd)) + log(dnorm(beta[i-1,3],0,prior_sd)) + log(dnorm(beta[i-1,4],0,prior_sd))
  
  beta_alpha = sum(log(dpois(Y,Mux(Xmat,beta_star))),base=exp(1)) + dnorm_star - sum(log(dpois(Y,Mux(Xmat,beta[i-1,]))),base=exp(1)) - dnorm_pre
  
  if (log(runif(1), base = exp(1)) < beta_alpha) {
    beta[i,] <-beta_star
    beta_count <- beta_count + 1
  }
  else beta[i,] <- beta[i-1,]
  
}

#trace plot
par(mfrow=c(2,2))
plot(beta[,1],type="l")
plot(beta[,2],type="l")
plot(beta[,3],type="l")
plot(beta[,4],type="l")

#acceptance probabilty
accept_prob = beta_count/niter
accept_prob

#density plot
par(mfrow=c(2,2))
plot(density(beta[,1]))
plot(density(beta[,2]))
plot(density(beta[,3]))
plot(density(beta[,4]))

#95% HPD-interval
library(HDInterval)
hdi(beta[,1], credMass = 0.95)
hdi(beta[,2], credMass = 0.95)
hdi(beta[,3], credMass = 0.95)
hdi(beta[,4], credMass = 0.95)

#effective sample size
effectiveSize(beta[,1])
effectiveSize(beta[,2])
effectiveSize(beta[,3])
effectiveSize(beta[,4])

#posterior mean
mean(beta[,1])
mean(beta[,2])
mean(beta[,3])
mean(beta[,4])

#3 Nimble package implementation 
library(mvtnorm)
n=1000
p=4

Xmat = rmvnorm(n,rep(0,4),diag(4))
trueBeta = c(0.5,-0.5,0,1)
Y = rpois(n,exp(Xmat%*%trueBeta))

model_glm <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    Mu[i] <- exp(XB[i])
    Y[i] ~ dpois(Mu[i])
  }
  XB[1:n] <- X[1:n,1:p]%*%beta[1:p]
  
  # Parameter Model
  beta[1:p] ~ dmnorm(mean=M[1:p], cov=Cov[1:p,1:p])
  
})

niter <- 1e4
consts <- list(n=dim(Xmat)[1], p=dim(Xmat)[2], X=Xmat, M=rep(0,dim(Xmat)[2]), Cov=diag(dim(Xmat)[2]))
dat <- list(Y=Y)
inits <- list(beta=rep(0,dim(Xmat)[2]))

# Run MCMC
pt<-proc.time()
samples_glm  <- nimbleMCMC(model_glm, data = dat, inits = inits,
                           constants=consts,
                           monitors = c("beta"),
                           samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE,
                           niter = niter, nburnin = 0, thin=1, nchains = 1)   
ptFinal_glm<-proc.time()-pt
ptFinal_glm

dim(samples_glm)

par(mfrow=c(2,2))
for(i in 1:dim(samples_glm)[2]){ ts.plot(samples_glm[,i]) }

par(mfrow=c(2,2))
for(i in 1:dim(samples_glm)[2]){ plot(density(samples_glm[,i])) }  

effectiveSize(as.mcmc(samples_glm))
HPDinterval(as.mcmc(samples_glm))

#posterior mean
mean(samples_glm[,1])
mean(samples_glm[,2])
mean(samples_glm[,3])
mean(samples_glm[,4])

#acceptance probabilty
beta_count_nimble <- 0

for(i in (2:niter)){
  if(samples_glm[i,1] != samples_glm[i-1,1]){
    beta_count_nimble = beta_count_nimble + 1
  }
}

accept_prob_nimble = beta_count_nimble/niter
accept_prob_nimble

#4 adaptMCMC package implementation 
library(adaptMCMC)

init.pars = c(beta0=0, beta1=0, beta2=0, beta3=0)

logPost = function(pars) {
  with(as.list(pars), {
    beta = pars
    Mu <- exp(Xmat%*%beta)
    # likelihood + prior 
    sum(dpois(Y, Mu, log=TRUE)) + sum(dnorm(beta,mean=0,sd=1,log=TRUE))
  })
}

out.mcmc <- MCMC(p=logPost, n=1e4, init=init.pars, scale=c(0.1,0.1,0.1,0.1), adapt=TRUE, acc.rate=.3)

par(mfrow=c(2,2))
for(i in 1:dim(out.mcmc$samples)[2]){ ts.plot(out.mcmc$samples[,i]) }

par(mfrow=c(2,2))
for(i in 1:dim(out.mcmc$samples)[2]){ plot(density(out.mcmc$samples[,i])) }

HPDinterval(as.mcmc(out.mcmc$samples))
effectiveSize(as.mcmc(out.mcmc$samples))

#posterior mean
mean(out.mcmc$samples[,1])
mean(out.mcmc$samples[,2])
mean(out.mcmc$samples[,3])
mean(out.mcmc$samples[,4])

#acceptance probabilty
beta_count_adapt <- 0

for(i in (2:niter)){
  if(out.mcmc$samples[i,1] != out.mcmc$samples[i-1,1]){
    beta_count_adapt = beta_count_adapt + 1
  }
}

accept_prob_adapt = beta_count_adapt/niter
accept_prob_adapt

#5 Comparing beta
par(mfrow=c(1,1))
#beta0
plot(density(beta[,1]), main="Beta0")
lines(density(samples_glm[,1]),col=3)
lines(density(out.mcmc$samples[,1]),col=4)
abline(v=trueBeta[1], col="red")
legend("topleft", legend=c("handMCMC","nimble","adaptMCMC"),pch=1, col=c(1,3,4))

#beta1
plot(density(beta[,2]), main="Beta1")
lines(density(samples_glm[,2]),col=3)
lines(density(out.mcmc$samples[,2]),col=4)
abline(v=trueBeta[2], col="red")
legend("topright", legend=c("handMCMC","nimble","adaptMCMC"),pch=1, col=c(1,3,4))

#beta2
plot(density(beta[,3]), main="Beta2")
lines(density(samples_glm[,3]),col=3)
lines(density(out.mcmc$samples[,3]),col=4)
abline(v=trueBeta[3], col="red")
legend("topright", legend=c("handMCMC","nimble","adaptMCMC"),pch=1, col=c(1,3,4))

#beta3
plot(density(beta[,4]), main="Beta3")
lines(density(samples_glm[,4]),col=3)
lines(density(out.mcmc$samples[,4]),col=4)
abline(v=trueBeta[4], col="red")
legend("topleft", legend=c("handMCMC","nimble","adaptMCMC"),pch=1, col=c(1,3,4))

