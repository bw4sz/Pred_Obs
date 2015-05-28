setwd("C:/Users/Ben/Documents/Pred_Obs/Bayesian")

sink("Single.jags")

cat("model{
    for(i in 1:length(y)){
    y[i] ~ dbern(p[i])
    logit(p[i]) <- alpha + beta * dist[i] + gamma * pow(dist[i],2)
    }
    
    alpha ~ dnorm(0,0.001)
    beta ~ dnorm(0,0.001)
    gamma ~ dnorm(0,0.001)
    
    }",fill = TRUE)

sink()
