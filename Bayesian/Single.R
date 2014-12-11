setwd("C:/Users/Ben/Documents/Pred_Obs/Bayesian")

sink("Single.jags")

cat("model{
    for(i in 1:length(y)){
    y[i] ~ dbern(p[i])
    logit(p[i]) <- intercept + linear * dist[i] + polynomial * pow(dist[i],2)
    }
    
    intercept ~ dnorm(0,0.001)
    linear ~ dnorm(0,0.001)
    polynomial ~ dnorm(0,0.001)
    
    }",fill = TRUE)

sink()
