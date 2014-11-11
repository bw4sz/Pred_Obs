Non-linear competition coefficients for phylogenetic community structure analysis
========================================================

Using the pglmm.sim code to create assemblages where adaptation to an environmental gradient is phylogenetically correlated, but competition increases with relatedness. This is model III in the Ives and Helmus 2011 paper. 

Aim
---
To create a polynomial function to describe the competition coefficient, such that distantly related species compete exponentially less than closely related. To me, this makes biological sense, since a small difference when you are _closely_ related should make more of a difference in terms of niche space as the same distance when you are very _distantly related_.

Approach
---
Reviewing the code for the pglmm.sim function in the picante package, it appears this code be accomplished in one of two areas. Either in the code chunk for model III here (from picante function):


```r
Xscale <- 1
Mscale <- 0
Vscale1 <- 1
Vscale2 <- 1
compscale <- compscale
b0scale <- 0
Vcomp <- solve(V, diag(nspp))
Vcomp <- Vcomp/max(Vcomp)
Vcomp <- compscale * Vcomp
iDcomp <- t(chol(Vcomp))
colnames(Vcomp) <- rownames(Vcomp)
```


or under the repulsion flag:


```r
if (repulseflag == 1) {
    bcomp <- NULL
    for (i in 1:m) {
        bcomp <- cbind(bcomp, iDcomp %*% rnorm(nspp))
    }
    bcomp0 <- 0
    Xcomp <- exp(bcomp0 + bcomp)/(1 + exp(bcomp0 + bcomp))
    X <- X * t(Xcomp)
}
```


Bring in data
----

```r
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
plot(trx, cex = 0.3)
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 



Testing the relationship between 'Vcomp', probability of occurrence and cophenetic distance 
----


```r

V <- vcv.phylo(trx, corr = TRUE)
nspp <- dim(V)[1]
compscale = 100

Vcomp <- solve(V, diag(nspp))
Vcomp <- Vcomp/max(Vcomp)
Vcomp <- compscale * Vcomp
iDcomp <- t(chol(Vcomp))
colnames(Vcomp) <- rownames(Vcomp)


Vm <- melt(Vcomp)

colnames(Vm) <- c("To", "From", "Vcomp")


pco <- cophenetic(trx)
# remove within species relatedness
diag(pco) <- NA

co <- melt(pco)

colnames(co) <- c("To", "From", "cophenetic")

dat <- merge(Vm, co)

head(dat)
```

```
##                      To                   From      Vcomp cophenetic
## 1 Adelomyia.melanogenys  Adelomyia.melanogenys  2.6898487         NA
## 2 Adelomyia.melanogenys Aglaeactis.cupripennis -0.0078686     0.4446
## 3 Adelomyia.melanogenys      Aglaeactis.pamela -0.0078686     0.4446
## 4 Adelomyia.melanogenys Aglaiocercus.coelestis -0.3588613     0.2789
## 5 Adelomyia.melanogenys     Aglaiocercus.kingi -0.3588613     0.2789
## 6 Adelomyia.melanogenys      Amazilia.amabilis -0.0002276     0.4649
```

```r

ggplot(dat, aes(x = cophenetic, Vcomp)) + geom_point()
```

```
## Warning: Removed 170 rows containing missing values (geom_point).
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


First attempt to 
----
