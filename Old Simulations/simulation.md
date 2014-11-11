Simulated phylogenetic community matrices to test the response of logistic regression
========================================================

I have modified the simulation data from ives and helmus (2011) and fit our hummingbird data. The goal of the simulation will be to test the response of a logistic model under simultaneous abiotic filtering and competition. Immense credit must go to Matt Helmus for the very elegant structure and code. 

I've edited the code to constrain the richness of the assemblages. Unconstrained, the model would create unrealistically large assemblages, which would dilute the effect of competition by allowing ample spots for rare co-occurence to happen. Instead of drawing randomly against the environmental value, the assemblages are first drawn a richness value from a poisson with lambda = 9, which is the mean richness of hummingbird assemblages in the area (Graham 2009, Weinstein accepted 2014)

Read in data
------------


```r

opts_chunk$set(warning = FALSE, message = FALSE, echo = TRUE)

require(picante)
require(ggplot2)
require(reshape)
# Set dropbox path
droppath <- "C:/Users/Ben/Dropbox/"

# Set git path
gitpath <- "C:/Users/Ben/Documents/Pred_Obs/"

# source input functions
source(paste(gitpath, "pglmmrichess.R", sep = ""))

# Bring in phylogenetic and trait info Bring in Phylogenetic Data
trx <- read.nexus(paste(gitpath, "InputData\\ColombiaPhylogenyUM.tre", sep = ""))
spnames <- read.table(paste(gitpath, "InputData\\SpNameTree.txt", sep = ""), 
    sep = "\t", header = TRUE)

# Replace tip.label with Spnames# replace the tiplabels with periods, which
# is the biomod default
trx$tip.label <- gsub("_", ".", as.character(spnames$SpName))
ctrx <- cophenetic(trx)

# Bring in the assemblages, cleaned from JP
Sites <- read.csv(paste(gitpath, "InputData//Sites.csv", sep = ""), row.names = 1)
plot(trx, cex = 0.3)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 


Our original data has 170 species and 229 localities. After data cleaning steps in the manuscript, we end up with 100 species and 200 localities.

__Model Testing__

I test the logistic model against the Ives and Helmus (2011) Model II and III simulations. Text from the original paper is below each model.


Model II
---

"Phylogenetically closely related species are often assumed to be ecologically similar, and therefore, if there is a set of strong environmental drivers determining the distribution of species, phylogenetically related species that share common responses to the environmental drivers should be more likely to co-occur."


```r
modelflag = 2
sim.dat <- pglmm.sim(trx, nsites = 10, modelflag = modelflag, second.env = FALSE, 
    compscale = 2)
```




```r
require(ggplot2)
require(reshape)

siteXsppm <- melt(sim.dat$Y)
```



We should briefly test that the simulation is working - there should be phylogenetic signal (p > .05) for Blomberg K following:

Blomberg, S. P., T. Garland Jr., A. R. Ives (2003) Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. Evolution, 57, 717–745.



```r
require(phytools)
phylosig(trx, sim.dat$bspp1, test = TRUE)
```

```
## $K
## [1] 1.121
## 
## $P
## [1] 0.001
```


In addition, for each model, take two sister species and inspect the correlation in probability of co-occurence. For Model I this should be very high.

__Aglaiocercus kingi__

![](http://farm4.static.flickr.com/3010/3093879517_c33fb350d1.jpg)

__Aglaiocercus coelestis__

![](http://upload.wikimedia.org/wikipedia/commons/7/79/Violet-tailed_Sylph_%28Aglaiocercus_coelestis%29.jpg)

Correlation among probability of occurence between _Aglaiocercus kingi_ and _Aglaiocercus coelestis_ is 0.963 

_Compute closest phylogenetic neighbor for each assemblage_


```r
sp.lists <- apply(sim.dat$Y, 1, function(x) {
    names(x[x == 1])
})
names(sp.lists) <- 1:nrow(sim.dat$Y)

colnames(siteXsppm) <- c("Locality", "Species", "P_A")

# Seperate out presence and absence, see
PA_phylo <- rbind.fill(apply(siteXsppm, 1, function(f) {
    # The function which is begin repeated for all communities, in both present
    # and absent, for all species get the related species (mean patristic
    # distance,median, and sum of distances)
    other <- unlist(sp.lists[names(sp.lists) %in% as.numeric(f[["Locality"]])])
    other.list <- other[!other %in% f[["Species"]]]
    
    # closest cophenetic neighbor
    hum.min <- min(ctrx[rownames(ctrx) %in% f[["Species"]], colnames(ctrx) %in% 
        other.list], na.rm = TRUE)
    
    return(data.frame(Species = f[["Species"]], Locality = f[["Locality"]], 
        Phylo.Relatedness = hum.min, P_A = f[["P_A"]]))
}))
```


_Fit a logistic model_


```r

# make numeric, not states
PA_phylo$P_A <- as.numeric(as.character(PA_phylo$P_A))

modII <- ggplot(PA_phylo, aes(x = Phylo.Relatedness, y = P_A)) + geom_point() + 
    geom_smooth(method = "glm", family = "binomial", formula = (y ~ (poly(x, 
        2))))
modII
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



Model III
---
The third model we consider assumes that phylogenetically
closely related species might show similar responses to
an environmental factor yet after factoring
out this attraction phylogenetically closely related
species tend not to co-occur. This scenario could
arise if there were simultaneous environmental ﬁltering,
with related species responding similarly to the ﬁlter,
and competition in which related species are more likely
to exclude each other


```r
modelflag = 3
sim.dat <- pglmm.sim(trx, nsites = 10, modelflag = modelflag, second.env = FALSE, 
    compscale = 1)
```




```r
require(ggplot2)
require(reshape)

siteXsppm <- melt(sim.dat$Y)
```


For each model we will take two siste species and inspect there correlation of probability of co-occurence. For Model III we think this should be reduced compared to model I and II.

Correlation among probability of occurence between _Aglaiocercus kingi_ and _Aglaiocercus coelestis_ is 0.9073 

_Compute closest phylogenetic neighbor for each assemblage_


```r
sp.lists <- apply(sim.dat$Y, 1, function(x) {
    names(x[x == 1])
})
names(sp.lists) <- 1:nrow(sim.dat$Y)

colnames(siteXsppm) <- c("Locality", "Species", "P_A")

# Seperate out presence and absence, see
PA_phylo <- rbind.fill(apply(siteXsppm, 1, function(f) {
    # The function which is begin repeated for all communities, in both present
    # and absent, for all species get the related species (mean patristic
    # distance,median, and sum of distances)
    other <- unlist(sp.lists[names(sp.lists) %in% as.numeric(f[["Locality"]])])
    other.list <- other[!other %in% f[["Species"]]]
    
    # closest cophenetic neighbor
    hum.min <- min(ctrx[rownames(ctrx) %in% f[["Species"]], colnames(ctrx) %in% 
        other.list], na.rm = TRUE)
    
    return(data.frame(Species = f[["Species"]], Locality = f[["Locality"]], 
        Phylo.Relatedness = hum.min, P_A = f[["P_A"]]))
}))
```


_Fit a logistic model_


```r

# make numeric, not states
PA_phylo$P_A <- as.numeric(as.character(PA_phylo$P_A))

modIII <- ggplot(PA_phylo, aes(x = Phylo.Relatedness, y = P_A)) + geom_point() + 
    geom_smooth(method = "glm", family = "binomial", formula = (y ~ (poly(x, 
        2))))
modIII
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 


Compare models
======


```r
require(gridExtra)
multiplot(modII, modIII, cols = 3)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 


To Do
--------

1. Reduce assemblage size.
2. Non monotonic repulsion matrix? Its unclear the form of the function. 


Model II - fifty runs
=====================


```r

require(doSNOW)
require(foreach)
require(ggplot2)

cl <- makeCluster(3, "SOCK")
registerDoSNOW(cl)

modLoop <- foreach(i = 1:5, .errorhandling = "pass") %dopar% {
    require(picante)
    require(reshape)
    require(ggplot2)
    out <- modRun(2, nsites = 10, compscale = 1)
    return(list(out))
}

names(modLoop) <- 1:length(modLoop)
moddf <- melt(modLoop, id.vars = colnames(modLoop[[1]][[1]]))

modIIloop <- ggplot(moddf, aes(x = Phylo.Relatedness, y = P_A, col = L1)) + 
    geom_point() + geom_smooth(method = "glm", family = "binomial", formula = (y ~ 
    (poly(x, 2))), se = FALSE) + geom_smooth(method = "glm", family = "binomial", 
    formula = (y ~ (poly(x, 2))), se = FALSE, aes(group = 1), linetype = "dashed", 
    col = "black", size = 0.9) + labs("Iteration")
modIIloop
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 


Model III - fifty runs
=======================


```r

require(doSNOW)
require(foreach)
require(ggplot2)

cl <- makeCluster(3, "SOCK")
registerDoSNOW(cl)
modLoop <- foreach(i = 1:5, .errorhandling = "pass") %dopar% {
    require(picante)
    require(reshape)
    require(ggplot2)
    out <- modRun(3, nsites = 10, compscale = 1)
    return(list(out))
}

names(modLoop) <- 1:length(modLoop)
moddf <- melt(modLoop, id.vars = colnames(modLoop[[1]][[1]]))

modIIIloop <- ggplot(moddf, aes(x = Phylo.Relatedness, y = P_A, col = L1)) + 
    geom_point() + geom_smooth(method = "glm", family = "binomial", formula = (y ~ 
    (poly(x, 2))), se = FALSE) + geom_smooth(method = "glm", family = "binomial", 
    formula = (y ~ (poly(x, 2))), se = FALSE, aes(group = 1), linetype = "dashed", 
    col = "black", size = 0.9) + labs(col = "Iteration")
modIIIloop
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 


This shows a very similiar pattern to our observed data. However, only a few runs show the reduced co-occurence among very closely related species. My hunch is this is because of the very large assemblage size in the simulations. On average assemblages have ~40-50 species, where in our observed data they have ~9-15. Interestingly, if we ratchet up the competition scaling, ie, the amount of repulsion, we do see curves that match the observed pattern.

__Model III, fifty runs, increasing competition with each run__


```r
require(doSNOW)
require(foreach)
require(ggplot2)

cl <- makeCluster(3, "SOCK")
registerDoSNOW(cl)
modLoop <- foreach(i = 1:5, .errorhandling = "pass") %dopar% {
    require(picante)
    require(reshape)
    require(ggplot2)
    out <- modRun(3, nsites = 10, compscale = i * 10)
    return(list(out))
}

names(modLoop) <- 1:length(modLoop) * 10
moddf <- melt(modLoop, id.vars = colnames(modLoop[[1]][[1]]))

modIIIloop <- ggplot(moddf, aes(x = Phylo.Relatedness, y = P_A, col = as.numeric(L1))) + 
    geom_point() + geom_smooth(method = "glm", family = "binomial", formula = (y ~ 
    (poly(x, 2))), se = FALSE, aes(group = L1)) + geom_smooth(method = "glm", 
    family = "binomial", formula = (y ~ (poly(x, 2))), se = FALSE, aes(group = 1), 
    linetype = "dashed", col = "black", size = 0.9) + labs("Competition Scaling")

modIIIloop
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 


__Model III, fifty runs, increasing number of sites with each run__


```r
require(doSNOW)
require(foreach)
require(ggplot2)

cl <- makeCluster(3, "SOCK")
registerDoSNOW(cl)
modLoop <- foreach(i = 1:10, .errorhandling = "pass") %dopar% {
    require(picante)
    require(reshape)
    require(ggplot2)
    out <- modRun(3, nsites = i * 5, compscale = 1)
    return(list(out))
}
stopCluster(cl)

names(modLoopa) <- 1:length(modLoopa) * 5
```

```
## Error: object 'modLoopa' not found
```

```r
moddfa <- melt(modLoopa, id.vars = colnames(modLoopa[[1]][[1]]))
```

```
## Error: object 'modLoopa' not found
```

```r

modIIIloopa <- ggplot(moddfa, aes(x = Phylo.Relatedness, y = P_A, col = as.numeric(L1))) + 
    geom_point() + geom_smooth(method = "glm", family = "binomial", formula = (y ~ 
    (poly(x, 2))), se = FALSE, aes(group = L1)) + geom_smooth(method = "glm", 
    family = "binomial", formula = (y ~ (poly(x, 2))), se = FALSE, aes(group = 1), 
    linetype = "dashed", col = "black", size = 0.9) + scale_color_continuous("Number of Sites", 
    low = "blue", high = "red")
```

```
## Error: object 'moddfa' not found
```

```r

modIIIloopa
```

```
## Error: object 'modIIIloopa' not found
```



Polynomial function for competition
-----

Achieved through squaring the Vcomp covariance matrix of trait values for species. This leads to increased competition among closely related species as compared to distance related species.

_One run_



```r
modelflag = 3
sim.dat <- pglmm.simM2(trx, nsites = 10, modelflag = modelflag, second.env = FALSE, 
    compscale = 1)
```

```
## Error: the leading minor of order 2 is not positive definite
```




```r
require(ggplot2)
require(reshape)

siteXsppm <- melt(sim.dat$Y)
```


For each model we will take two siste species and inspect there correlation of probability of co-occurence. For Model III we think this should be reduced compared to model I and II.

Correlation among probability of occurence between _Aglaiocercus kingi_ and _Aglaiocercus coelestis_ is 0.9073 

_Compute closest phylogenetic neighbor for each assemblage_


```r
sp.lists <- apply(sim.dat$Y, 1, function(x) {
    names(x[x == 1])
})
names(sp.lists) <- 1:nrow(sim.dat$Y)

colnames(siteXsppm) <- c("Locality", "Species", "P_A")

# Seperate out presence and absence, see
PA_phylo <- rbind.fill(apply(siteXsppm, 1, function(f) {
    # The function which is begin repeated for all communities, in both present
    # and absent, for all species get the related species (mean patristic
    # distance,median, and sum of distances)
    other <- unlist(sp.lists[names(sp.lists) %in% as.numeric(f[["Locality"]])])
    other.list <- other[!other %in% f[["Species"]]]
    
    # closest cophenetic neighbor
    hum.min <- min(ctrx[rownames(ctrx) %in% f[["Species"]], colnames(ctrx) %in% 
        other.list], na.rm = TRUE)
    
    return(data.frame(Species = f[["Species"]], Locality = f[["Locality"]], 
        Phylo.Relatedness = hum.min, P_A = f[["P_A"]]))
}))
```


_Fit a logistic model_


```r

# make numeric, not states
PA_phylo$P_A <- as.numeric(as.character(PA_phylo$P_A))

modIII <- ggplot(PA_phylo, aes(x = Phylo.Relatedness, y = P_A)) + geom_point() + 
    geom_smooth(method = "glm", family = "binomial", formula = (y ~ (poly(x, 
        2))))
modIII
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19.png) 


__fifty runs__ 


```r
require(doSNOW)
require(foreach)
require(ggplot2)

cl <- makeCluster(3, "SOCK")
registerDoSNOW(cl)
modLoopS <- foreach(i = 1:5, .errorhandling = "pass") %dopar% {
    require(picante)
    require(reshape)
    require(ggplot2)
    out <- modRun2(3, nsites = 10, compscale = i)
    return(list(out))
}
stopCluster(cl)

names(modLoopb) <- 1:length(modLoopS)
```

```
## Error: object 'modLoopb' not found
```

```r
moddfS <- melt(modLoopS, id.vars = colnames(modLoopS[[1]][[1]]))
```

```
## Error: cannot coerce class "c("simpleError", "error", "condition")" to a
## data.frame
```

```r

modIIIloopS <- ggplot(moddfS, aes(x = Phylo.Relatedness, y = P_A, col = as.numeric(L1))) + 
    geom_point() + geom_smooth(method = "glm", family = "binomial", formula = (y ~ 
    (poly(x, 2))), se = FALSE, aes(group = L1)) + geom_smooth(method = "glm", 
    family = "binomial", formula = (y ~ (poly(x, 2))), se = FALSE, aes(group = 1), 
    linetype = "dashed", col = "black", size = 0.9) + scale_color_continuous("Number of sites", 
    low = "blue", high = "red")
```

```
## Error: object 'moddfS' not found
```

```r

modIIIloopS
```

```
## Error: object 'modIIIloopS' not found
```



Compare models of multiple runs
======


```r
require(gridExtra)
multiplot(modIIloop, modIIIloop, modIIIloopa, modIIIloopS, cols = 3)
```

```
## Error: object 'modIIIloopa' not found
```



Save workspace


```r
save.image(paste(droppath, "Simulationworkspace.RData", sep = ""))
```

