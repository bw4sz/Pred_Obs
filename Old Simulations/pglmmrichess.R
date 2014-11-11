##new source function for pglmm sim to modify

#loop function

#define function
modRun<-function(modelflag,nsites,compscale=1){
  sim.dat<-pglmm.sim(trx,nsites=10,modelflag=modelflag,second.env=FALSE,compscale=1)
  
  siteXsppm<-melt(sim.dat$Y)
  sp.lists<-apply(sim.dat$Y,1,function(x){
    names(x[x==1])
  })
  names(sp.lists)<-1:nrow(sim.dat$Y)
  
  colnames(siteXsppm)<-c("Locality","Species","P_A")
  
  #Seperate out presence and absence, see 
  PA_phylo<-rbind.fill(apply(siteXsppm,1,function(f){ 
    
    #The function which is begin repeated for all communities, in both present and absent, for all species
    
    #get the related species (mean patristic distance,median, and sum of distances) 
    
    other<-unlist(sp.lists[names(sp.lists) %in% as.numeric(f[["Locality"]])])
    other.list<-other[!other %in% f[["Species"]]]
    
    #closest cophenetic neighbor
    hum.min<-min(ctrx[rownames(ctrx) %in% f[["Species"]],colnames(ctrx) %in% other.list],na.rm=TRUE)    
    
    return(data.frame(Species=f[["Species"]],Locality=f[["Locality"]],Phylo.Relatedness=hum.min,P_A=f[["P_A"]]))
  }))
  
  PA_phylo$P_A<-as.numeric(as.character(PA_phylo$P_A))
  
  modI<-ggplot(PA_phylo,aes(x=Phylo.Relatedness,y=P_A)) + geom_point() + geom_smooth(method="glm",family="binomial",formula=(y~(poly(x,2)))) 
  
  print(modI)
  
  return(PA_phylo)}


#now modloop for the covariance squared function

#define function
modRun2<-function(modelflag,nsites,compscale=1){
  sim.dat<-pglmm.simM2(trx,nsites=10,modelflag=modelflag,second.env=FALSE,compscale=100)
  
  siteXsppm<-melt(sim.dat$Y)
  sp.lists<-apply(sim.dat$Y,1,function(x){
    names(x[x==1])
  })
  names(sp.lists)<-1:nrow(sim.dat$Y)
  
  colnames(siteXsppm)<-c("Locality","Species","P_A")
  
  #Seperate out presence and absence, see 
  PA_phylo<-rbind.fill(apply(siteXsppm,1,function(f){ 
    
    #The function which is begin repeated for all communities, in both present and absent, for all species
    
    #get the related species (mean patristic distance,median, and sum of distances) 
    
    other<-unlist(sp.lists[names(sp.lists) %in% as.numeric(f[["Locality"]])])
    other.list<-other[!other %in% f[["Species"]]]
    
    #closest cophenetic neighbor
    hum.min<-min(ctrx[rownames(ctrx) %in% f[["Species"]],colnames(ctrx) %in% other.list],na.rm=TRUE)    
    
    return(data.frame(Species=f[["Species"]],Locality=f[["Locality"]],Phylo.Relatedness=hum.min,P_A=f[["P_A"]]))
  }))
  
  PA_phylo$P_A<-as.numeric(as.character(PA_phylo$P_A))
  
  modI<-ggplot(PA_phylo,aes(x=Phylo.Relatedness,y=P_A)) + geom_point() + geom_smooth(method="glm",family="binomial",formula=(y~(poly(x,2)))) 
  
  print(modI)
  
  return(PA_phylo)}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#covariance squared function

pglmm.simM2<-function (tree=trx, nsites = 30, modelflag = 3, second.env = FALSE, 
          compscale = 1) 
{
  bspp2 <- NULL
  Vcomp <- NULL
  envirogradflag2 <- 0
  if (is(tree)[1] == "phylo") {
    if (is.null(tree$edge.length)) {
      tree <- compute.brlen(tree, 1)
    }
    V <- vcv.phylo(tree, corr = TRUE)
  }
  else {
    V <- tree
  }
  nspp <- dim(V)[1]
  
  if (modelflag == 3) {
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
    envirogradflag1 <- 1
    if (second.env) {
      envirogradflag2 <- 1
    }
    elimsitesflag <- 0
    repulseflag <- 1
  }
  
  Vforsim <- V
  iD <- t(chol(Vforsim))
  mx <- t(as.matrix((-(nsites)/2):(nsites/2)))
  m <- length(mx)
  if (envirogradflag1 == 1) {
    e <- iD %*% rnorm(nspp)
    e <- Vscale1 * (e - mean(e))/apply(e, 2, sd)
    bspp1 <- e
    X <- 1/(1 + exp(-(b0scale * array(1, c(m, 1)) %*% rnorm(nspp) + 
                        t(mx) %*% t(e))))
    X <- Xscale * X
    Xsmooth <- X
    X1 <- diag(1 - Mscale * runif(m)) %*% X
  }
  if (envirogradflag2 == 1) {
    e <- iD %*% rnorm(nspp)
    e <- Vscale2 * (e - mean(e))/apply(e, 2, sd)
    bspp2 <- e
    mx2 <- as.matrix(mx[sample(m)])
    X <- 1/(1 + exp(-(b0scale * array(1, c(m, 1)) %*% rnorm(nspp) + 
                        mx2 %*% t(e))))
    X <- Xscale * X
    Xsmooth <- Xsmooth * X
    X2 <- diag(1 - Mscale * runif(m)) %*% X
    X <- X1 * X2
  }
  else {
    X <- X1
  }
  if (repulseflag == 1) {
    bcomp <- NULL
    for (i in 1:m) {
      bcomp <- cbind(bcomp, iDcomp %*% rnorm(nspp))
    }
    bcomp0 <- 0
    Xcomp <- exp(bcomp0 + bcomp)/(1 + exp(bcomp0 + bcomp))
    X <- X * t(Xcomp)
  }
  Y <- matrix(0, ncol = nspp, nrow = m)
  
  #use species probabilities to pick species
  richnessS<-rpois(nrow(X),9)
  
  #draw species
  richC<-melt(sapply(1:nrow(X),function(g){
    sample(colnames(X),prob=X[g,],size=richnessS[g])    
  }))
  
  #create siteXspp 
  Y<-t(as.matrix(table(richC)))
  #Y[matrix(runif(nspp * m), ncol = nspp) < X] <- 1
  #colnames(Y) <- colnames(X)
  pick <- (rowSums(Y) > 0)
  Y <- Y[pick, ]
  mx <- mx[pick]
  m <- length(mx)
  if (elimsitesflag == 1) {
    pick <- (rowSums(Y) > 1)
    Y <- Y[pick, ]
    mx <- mx[pick]
    m <- length(mx)
  }
  return(list(Vphylo = V, Vcomp = Vcomp, Y = Y, X = X, u = mx, 
              bspp1 = bspp1, bspp2 = bspp2))
}
