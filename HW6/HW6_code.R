library(nimble)
library(coda)
library(mvtnorm)
library(invgamma)
library(Rlab)


#1-(b)
true_beta = c(0.5, -0.5, 1, -1, 0.7, 0, 0, 0, 0, 0)
dim(true_beta) = c(p,1)

Y <- rep(0,n)
dim(Y) <- c(n,1)

Mx = function(X,B){
  return(exp(X%*%B)/(1+exp(X%*%B)))
}

true_Mx = Mx(Xmat,true_beta)

for (i in 1:n){
  Y[i] <- rbinom(1,1,true_Mx[i])
}

head(Y)

#2 MCMC in hand
#setting
niter = 10000
beta_count = 0
lambda_count = 0
sigma1_count = 0
sigma2_count = 0
prop_std = 0.2

#Initialize parameters
beta = matrix(rep(0,niter*p),niter,p)
lambda = matrix(rep(0,niter*p),niter,p)
sigma1 = matrix(rep(0,niter*p),niter,1)
sigma2 = matrix(rep(0,niter*p),niter,1)

beta_star = matrix(rep(0,p),1,p)
lambda_star = matrix(rep(0,p),1,p)

for (k in 1:p){
  lambda[1,k] = rbinom(1,1,1/2)
}

sigma1[1] = rinvgamma(1, 1, rate = 20)
sigma2[1] = rgamma(1, 1, rate = 20)

for (k in 1:p){
  beta[1,k] = lambda[1,k]*rnorm(1,0,sigma1[1]) + (1-lambda[1,k])*rnorm(1,0,sigma2[1])
}


#mcmc sampling
for (i in 2:niter) {
  
  #proposing betas
  for (k in 1:p){
    beta_star[k] = rnorm(1,0,prop_std)
  }
  
  #calculating alpha for proposed beta
  y_part = 0
  
  #Replacing NAN and Inf to 0
  logMx_betastar = log(Mx(Xmat,t(beta_star)))
  logMx_betastar[is.na(logMx_betastar)] <- 0
  logMx_betastar[is.infinite(logMx_betastar)] <- 0
  
  logMx_1betastar = log(1-Mx(Xmat,t(beta_star)))
  logMx_1betastar[is.na(logMx_1betastar)] <- 0
  logMx_1betastar[is.infinite(logMx_1betastar)] <- 0
  
  logMx_beta = log(Mx(Xmat,beta[i-1,]))
  logMx_beta[is.na(logMx_beta)] <- 0
  logMx_beta[is.infinite(logMx_beta)] <- 0
  
  logMx_1beta = log(1-Mx(Xmat,beta[i-1,]))
  logMx_1beta[is.na(logMx_1beta)] <- 0
  logMx_1beta[is.infinite(logMx_1beta)] <- 0
  
  y_part = y_part + t(Y)%*%logMx_betastar + t(1-Y)%*%logMx_1betastar
  y_part = y_part - t(Y)%*%logMx_beta - t(1-Y)%*%logMx_1beta
    
  beta_part = 0
  for (k in 1:p){
    beta_part = beta_part + log( dnorm(beta_star[k], 0, lambda[i-1,k]*sigma1[i-1] + (1-lambda[i-1,k])*sigma2[i-1])+ 1e-10 )
    beta_part = beta_part - log( dnorm(beta[i-1,k], 0, lambda[i-1,k]*sigma1[i-1] + (1-lambda[i-1,k])*sigma2[i-1])+ 1e-10 )
  }
  
  alpha = y_part + beta_part
  
  #updating beta
  if ( log(runif(1), base = exp(1)) < alpha ) {
    beta[i,] <- beta_star
    beta_count <- beta_count + 1
  }else {
    beta[i,] <- beta[i-1,]
  }
  
  #proposing lambdas
  for (k in 1:p){
    lambda_star[k] = rbinom(1,1,1/2)
  }
  
  #influenced betas
  for (k in 1:p){
    beta_star[k] = lambda_star[k]*rnorm(1,beta[i,k],sigma1[i-1]) + (1-lambda_star[k])*rnorm(1,beta[i,k],sigma2[i-1])
  }
  
  #calculating alpha for lambdas
  y_part = 0
  
  #Replacing NAN and Inf to 0
  logMx_betastar = log(Mx(Xmat,t(beta_star)))
  logMx_betastar[is.na(logMx_betastar)] <- 0
  logMx_betastar[is.infinite(logMx_betastar)] <- 0
  
  logMx_1betastar = log(1-Mx(Xmat,t(beta_star)))
  logMx_1betastar[is.na(logMx_1betastar)] <- 0
  logMx_1betastar[is.infinite(logMx_1betastar)] <- 0
  
  logMx_beta = log(Mx(Xmat,beta[i,]))
  logMx_beta[is.na(logMx_beta)] <- 0
  logMx_beta[is.infinite(logMx_beta)] <- 0
  
  logMx_1beta = log(1-Mx(Xmat,beta[i,]))
  logMx_1beta[is.na(logMx_1beta)] <- 0
  logMx_1beta[is.infinite(logMx_1beta)] <- 0
  
  y_part = y_part + t(Y)%*%logMx_betastar + t(1-Y)%*%logMx_1betastar
  y_part = y_part - t(Y)%*%logMx_beta - t(1-Y)%*%logMx_1beta
  
  beta_part = 0
  for (k in 1:p){
    beta_part = beta_part + log( dnorm(beta_star[k], 0, lambda_star[k]*sigma1[i-1] + (1-lambda_star[k])*sigma2[i-1])+ 1e-10 )
    beta_part = beta_part - log( dnorm(beta[i,k], 0, lambda[i-1,k]*sigma1[i-1] + (1-lambda[i-1,k])*sigma2[i-1])+ 1e-10 )
  }
  
  lambda_part = 0
  for (k in 1:p){
    lambda_part = lambda_part + log( dbinom(lambda_star[k],1,1/2) )
    lambda_part = lambda_part - log( dbinom(lambda[i-1,k],1,1/2) )
  }
  
  alpha = 0
  alpha = y_part + beta_part + lambda_part
  
  #updating lambda
  if ( log(runif(1), base = exp(1)) < alpha ) {
    lambda[i,] <- lambda_star
    lambda_count <- lambda_count + 1
  }else {
    lambda[i,] <- lambda[i-1,]
  }
  
  #proposing sigma1
  sigma1_star = rnorm(1,0,prop_std)
  
  #influenced betas
  for (k in 1:p){
    beta_star[k] = lambda[i,k]*rnorm(1,beta[i,k],sigma1_star) + (1-lambda[i,k])*rnorm(1,beta[i,k],sigma2[i-1])
  }
  
  #calculating alpha for sigma1
  y_part = 0
  
  #Replacing NAN and Inf to 0
  logMx_betastar = log(Mx(Xmat,t(beta_star)))
  logMx_betastar[is.na(logMx_betastar)] <- 0
  logMx_betastar[is.infinite(logMx_betastar)] <- 0
  
  logMx_1betastar = log(1-Mx(Xmat,t(beta_star)))
  logMx_1betastar[is.na(logMx_1betastar)] <- 0
  logMx_1betastar[is.infinite(logMx_1betastar)] <- 0
  
  logMx_beta = log(Mx(Xmat,beta[i,]) + 1e+10)
  logMx_beta[is.na(logMx_beta)] <- 0
  logMx_beta[is.infinite(logMx_beta)] <- 0
  
  logMx_1beta = log(1-Mx(Xmat,beta[i,]))
  logMx_1beta[is.na(logMx_1beta)] <- 0
  logMx_1beta[is.infinite(logMx_1beta)] <- 0
  
  y_part = y_part + t(Y)%*%logMx_betastar + t(1-Y)%*%logMx_1betastar
  y_part = y_part - t(Y)%*%logMx_beta - t(1-Y)%*%logMx_1beta
  
  beta_part = 0
  tmp_pre = 0
  tmp_star = 0
  for (k in 1:p){
    tmp_star = log( dnorm(beta_star[k], 0, lambda[i,k]*sigma1_star + (1-lambda[i,k])*sigma2[i-1])+ 1e-10 )
    tmp_star[is.na(tmp_star)] <- 0
    tmp_star[is.infinite(tmp_star)] <- 0
    
    tmp_pre = log( dnorm(beta[i,k], 0, lambda[i,k]*sigma1[i-1] + (1-lambda[i,k])*sigma2[i-1])+ 1e-10 )
    tmp_pre[is.na(tmp_star)] <- 0
    tmp_pre[is.infinite(tmp_star)] <- 0
    
    beta_part = beta_part + tmp_star
    beta_part = beta_part - tmp_pre
  }
  
  sigma1_star = log( dinvgamma(sigma1_star, 1, scale = 20) )
  sigma1_pre = log( dinvgamma(sigma1[i-1], 1, scale = 20) )
  sigma1_star[is.na(sigma1_star)] <- 0
  sigma1_star[is.infinite(sigma1_star)] <- 0
  sigma1_pre[is.na(sigma1_pre)] <- 0
  sigma1_pre[is.infinite(sigma1_pre)] <- 0
  
  sigma1_part = sigma1_star - sigma1_pre
  
  alpha = 0
  alpha = y_part + beta_part + sigma1_part
  
  #updating sigma1
  if ( log(runif(1), base = exp(1)) < alpha ) {
    sigma1[i,] <- sigma1_star
    sigma1_count <- sigma1_count + 1
  }else {
    sigma1[i,] <- sigma1[i-1,]
  }
  
  #proposing sigma2
  sigma2_star = rnorm(1,0,prop_std)
  
  #influenced betas
  for (k in 1:p){
    beta_star[k] = lambda[i,k]*rnorm(1,beta[i,k],sigma1[i]) + (1-lambda[i,k])*rnorm(1,beta[i,k],sigma2_star)
  }
  
  #calculating alpha for sigma2
  y_part = 0
  
  #Replacing NAN and Inf to 0
  logMx_betastar = log(Mx(Xmat,t(beta_star)))
  logMx_betastar[is.na(logMx_betastar)] <- 0
  logMx_betastar[is.infinite(logMx_betastar)] <- 0
  
  logMx_1betastar = log(1-Mx(Xmat,t(beta_star)))
  logMx_1betastar[is.na(logMx_1betastar)] <- 0
  logMx_1betastar[is.infinite(logMx_1betastar)] <- 0
  
  logMx_beta = log(Mx(Xmat,beta[i,]))
  logMx_beta[is.na(logMx_beta)] <- 0
  logMx_beta[is.infinite(logMx_beta)] <- 0
  
  logMx_1beta = log(1-Mx(Xmat,beta[i,]))
  logMx_1beta[is.na(logMx_1beta)] <- 0
  logMx_1beta[is.infinite(logMx_1beta)] <- 0
  
  y_part = y_part + t(Y)%*%logMx_betastar + t(1-Y)%*%logMx_1betastar
  y_part = y_part - t(Y)%*%logMx_beta - t(1-Y)%*%logMx_1beta
  
  beta_part = 0
  tmp_pre = 0
  tmp_star = 0
  for (k in 1:p){
    tmp_star = log( dnorm(beta_star[k], 0, lambda[i,k]*sigma1[i] + (1-lambda[i,k])*sigma2_star)+ 1e-10 )
    tmp_star[is.na(tmp_star)] <- 0
    tmp_star[is.infinite(tmp_star)] <- 0
    
    tmp_pre = log( dnorm(beta[i,k], 0, lambda[i,k]*sigma1[i] + (1-lambda[i,k])*sigma2[i-1])+ 1e-10 )
    tmp_pre[is.na(tmp_star)] <- 0
    tmp_pre[is.infinite(tmp_star)] <- 0
    
    beta_part = beta_part + tmp_star
    beta_part = beta_part - tmp_pre
  }
  
  sigma2_star = log( dgamma(sigma2_star, 1, rate = 20) )
  sigma2_pre = log( dgamma(sigma2[i-1], 1, rate = 20) )
  sigma2_star[is.na(sigma2_star)] <- 0
  sigma2_star[is.infinite(sigma2_star)] <- 0
  sigma2_pre[is.na(sigma2_pre)] <- 0
  sigma2_pre[is.infinite(sigma2_pre)] <- 0
  
  sigma2_part = sigma2_star - sigma2_pre
  
  alpha=0
  alpha= y_part + beta_part + sigma2_part
  #updating sigma2
  if ( log(runif(1), base = exp(1)) < alpha ) {
    sigma2[i,] <- sigma2_star
    sigma2_count <- sigma2_count + 1
  }else {
    sigma2[i,] <- sigma2[i-1,]
  }
}

#End of hand-made MCMC

#trace plot
par(mfrow=c(1,1))
plot(beta[,1],type="l",col=1, ann=FALSE)
plot(beta[,2],type="l",col=1, ann=FALSE)
plot(beta[,3],type="l",col=1, ann=FALSE)
plot(beta[,4],type="l",col=1, ann=FALSE)
plot(beta[,5],type="l",col=1, ann=FALSE)
plot(beta[,6],type="l",col=1, ann=FALSE)
plot(beta[,7],type="l",col=1, ann=FALSE)
plot(beta[,8],type="l",col=1, ann=FALSE)
plot(beta[,9],type="l",col=1, ann=FALSE)
plot(beta[,10],type="l",col=1, ann=FALSE)

legend("bottomright",legend=c("beta1(0.5)","beta2(-0.5)","beta3(1)","beta4(-1)","beta5(0.7)","beta6(0)","beta7(0)","beta8(0)","beta9(0)","beta10(0)"),fill=c(1,2,3,4,5,6,7,8,9,10))

par(mfrow=c(1,2))
plot(sigma1,type="l",col=1, ann=FALSE)
plot(sigma2,type="l",col=2, ann=FALSE) 

#density plot
par(mfrow=c(1,1))
plot(density(beta[,8]),type="l",col=8, ann=FALSE)
points(density(beta[,1]),type="l",col=1) 
points(density(beta[,2]),type="l",col=2) 
points(density(beta[,3]),type="l",col=3) 
points(density(beta[,4]),type="l",col=4) 
points(density(beta[,5]),type="l",col=5) 
points(density(beta[,6]),type="l",col=6) 
points(density(beta[,7]),type="l",col=7) 
points(density(beta[,9]),type="l",col=9)
points(density(beta[,10]),type="l",col=10) 

legend("topleft",legend=c("beta1(0.5)","beta2(-0.5)","beta3(1)","beta4(-1)","beta5(0.7)","beta6(0)","beta7(0)","beta8(0)","beta9(0)","beta10(0)"),fill=c(1,2,3,4,5,6,7,8,9,10))


par(mfrow=c(1,2))
plot(density(sigma1),type="l",col=1, ann=FALSE)
plot(density(sigma2),type="l",col=2, ann=FALSE) 

#95% HPD-interval
library(HDInterval)
hdi(beta[,1], credMass = 0.95)
hdi(beta[,2], credMass = 0.95)
hdi(beta[,3], credMass = 0.95)
hdi(beta[,4], credMass = 0.95)
hdi(beta[,5], credMass = 0.95)
hdi(beta[,6], credMass = 0.95)
hdi(beta[,7], credMass = 0.95)
hdi(beta[,8], credMass = 0.95)
hdi(beta[,9], credMass = 0.95)
hdi(beta[,10], credMass = 0.95)

hdi(sigma1, credMass = 0.95)
hdi(sigma2, credMass = 0.95)


#posterior mean
mean(beta[,1])
mean(beta[,2])
mean(beta[,3])
mean(beta[,4])
mean(beta[,5])
mean(beta[,6])
mean(beta[,7])
mean(beta[,8])
mean(beta[,9])
mean(beta[,10])

mean(sigma1)
mean(sigma2)

#acceptance probabilty
beta_count/n
lambda_count/n
sigma1_count/n
sigma2_count/n

(beta_count+lambda_count+sigma1_count+sigma2_count)/n

library(coda)
#effective sample size
effectiveSize(beta[,1])
effectiveSize(beta[,2])
effectiveSize(beta[,3])
effectiveSize(beta[,4])
effectiveSize(beta[,5])
effectiveSize(beta[,6])
effectiveSize(beta[,7])
effectiveSize(beta[,8])
effectiveSize(beta[,9])
effectiveSize(beta[,10])

effectiveSize(sigma1)
effectiveSize(sigma2)

#Posterior inclusion probability
sum(lambda[,1])/n
sum(lambda[,2])/n
sum(lambda[,3])/n
sum(lambda[,4])/n
sum(lambda[,5])/n
sum(lambda[,6])/n
sum(lambda[,7])/n
sum(lambda[,8])/n
sum(lambda[,9])/n
sum(lambda[,10])/n

#3 Nimble package implementation 

dim(Y) = NULL

model_glm <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    Mu[i] <- exp(XB[i])/(1+exp(XB[i]))
    Y[i] ~ dbern(prob = Mu[i])
  }
  XB[1:n] <- X[1:n,1:p]%*%IBeta[1:p]
  
  # Parameter Model
  Sigma1 ~ dinvgamma(shape =1, rate =20)
  Sigma2 ~ dgamma(shape =1, rate = 20)
  
  for(i in 1:p){
    Lambda[i] ~ dbern(prob= 1/2)
    Beta1[i] ~ dnorm(mean=0, sd=Sigma1)
    Beta2[i] ~ dnorm(mean=0, sd=Sigma2)
    IBeta[i] <- Lambda[i]*Beta1[i] + (1-Lambda[i])*Beta2[i]
  }
})

niter <- 1e4
consts <- list(n=dim(Xmat)[1], p=dim(Xmat)[2], X=Xmat)
dat <- list(Y=Y)
inits <- list(IBeta=rep(0,dim(Xmat)[2]),Sigma1=1e-10,Sigma2=1e-10,Lambda=rep(1,dim(Xmat)[2]))

# Run MCMC
pt<-proc.time()
samples_glm  <- nimbleMCMC(model_glm, data = dat, inits = inits,
                           constants=consts,
                           monitors = c("IBeta","Sigma1","Sigma2","Lambda"),
                           samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE,
                           niter = niter, nburnin = 0, thin=1, nchains = 1)   
ptFinal_glm<-proc.time()-pt
ptFinal_glm

dim(samples_glm)

#Trace plot and Density plot
par(mfrow=c(1,1))
plot(samples_glm[,1], ann=FALSE)
plot(samples_glm[,2], ann=FALSE)
plot(samples_glm[,3], ann=FALSE)
plot(samples_glm[,4], ann=FALSE)
plot(samples_glm[,5], ann=FALSE)
plot(samples_glm[,6], ann=FALSE)
plot(samples_glm[,7], ann=FALSE)
plot(samples_glm[,8], ann=FALSE)
plot(samples_glm[,9], ann=FALSE)
plot(samples_glm[,10], ann=FALSE)
plot(samples_glm[,11], ann=FALSE)
plot(samples_glm[,12], ann=FALSE)

effectiveSize(as.mcmc(samples_glm))

HPDinterval(as.mcmc(samples_glm))

#posterior mean
for(k in (1:12)){
  if(k<11){
    cat("Posterior mean of beta",k,": ",mean(samples_glm[,k]),"\n")
  } else{
    cat("Posterior mean of Sigma",k-10,": ",mean(samples_glm[,10+k]),"\n")
  }
}

#acceptance probabilty

for(k in (1:12)){
  count_nimble <- 0
  for(i in (2:niter)){
    if(samples_glm[i,k] != samples_glm[i-1,k]){
      count_nimble = count_nimble + 1
    }
  }
  if(k < 11){
    cat("Acceptance rate for beta",k,": ",count_nimble/niter,"\n")
  }
  else{
    cat("Acceptance rate for sigma",k-10,": ",count_nimble/niter,"\n")
  }
}

#Posterior inclusion probability
for(k in (11:20)){
  count <- sum(samples_glm[,k])/niter
  cat("Posterior inclusion probability by lambda1",k-10,": ",count,"\n")
}