setwd("C:/Users/Ben/Documents/Pred_Obs/Bayesian")

sink("NullD.jags")

cat("model{
    for(i in 1:length(y)){
    y[i] ~ dbern(p[i])
    logit(p[i]) <- alpha[Species[i]] 
    }
    
    for (j in 1:s){
    alpha[j] ~ dnorm(intercept,tauIntercept)
    } 
    
    intercept ~ dnorm(0,0.001)
    tauIntercept ~ dgamma(0.001,0.001)
    sigmaIntercept<- pow(1/tauIntercept,.5)
    
    }",fill = TRUE)

sink()
