setwd("C:/Users/Ben/Documents/Pred_Obs/Bayesian")

sink("BPhylo.jags")

cat("model{
  for(i in 1:length(y)){
      y[i] ~ dbern(p[i])
      logit(p[i]) <- alpha[Species[i]] + beta[Species[i]] * dist[i] + gamma[Species[i]] * pow(dist[i],2)
    }

    for (j in 1:s){
    beta[j] ~ dnorm(slope,tauSlope)
    alpha[j] ~ dnorm(intercept,tauIntercept)
    gamma[j] ~ dnorm(slope2,tauSlope2)
    } 

    slope ~ dnorm(0,0.001)
    tauSlope ~ dgamma(0.001,0.001)

    slope2 ~ dnorm(0,0.001)
    tauSlope2 ~ dgamma(0.001,0.001)
    
    intercept ~ dnorm(0,0.001)
    tauIntercept ~ dgamma(0.001,0.001)

    sigmaSlope <- pow(1/tauSlope,.5)
    sigmaSlope2 <- pow(1/tauSlope2,.5)
    sigmaIntercept<- pow(1/tauIntercept,.5)

    }",fill = TRUE)

sink()
