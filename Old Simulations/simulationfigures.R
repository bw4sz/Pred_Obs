require(picante)
require(ggplot2)
require(reshape)
#Set dropbox path
#droppath<-"C:/Users/Ben/Dropbox/"

#Set git path
gitpath<-"C:/Users/Ben/Documents/Pred_Obs/"

#load from cluster, optionally

load(paste(gitpath,"Cluster/simulation.RData",sep=""))
#source input functions
source(paste(gitpath,"pglmmrichess.R",sep=""))

#Bring in phylogenetic and trait info
#Bring in Phylogenetic Data
trx<-read.nexus(paste(gitpath,"InputData\\ColombiaPhylogenyUM.tre",sep=""))
spnames<-read.table(paste(gitpath,"InputData\\SpNameTree.txt",sep=""), sep = "\t", header = TRUE)

#Replace tip.label with Spnames#
#replace the tiplabels with periods, which is the biomod default
trx$tip.label<-gsub("_",".",as.character(spnames$SpName))
ctrx<-cophenetic(trx)

#Bring in the assemblages, cleaned from JP
Sites<-read.csv(paste(gitpath,"InputData//Sites.csv",sep=""),row.names=1)
plot(trx,cex=.3)

##Model 2
require(doSNOW)
require(foreach)
require(ggplot2)

cl<-makeCluster(3,"SOCK")
registerDoSNOW(cl)

modLoop<-foreach(i=1:2,.errorhandling="pass") %do% {
  require(picante)
  require(reshape)
  require(ggplot2)
  out<-modRun(2,nsites=10,compscale=1)
  return(list(out))}

stopCluster(cl)

names(modLoop)<-1:length(modLoop)
moddf<-melt(modLoop,id.vars=colnames(modLoop[[1]][[1]]))

modIIloop<-ggplot(moddf,aes(x=Phylo.Relatedness,y=P_A,col=L1)) + geom_smooth(method="glm",family="binomial",formula=(y~(poly(x,2))),se=FALSE) + geom_smooth(method="glm",family="binomial",formula=(y~(poly(x,2))),se=FALSE,aes(group=1),linetype="dashed",col="black",size=.9) + labs("Iteration")
modIIloop

#model 3

cl<-makeCluster(3,"SOCK")
registerDoSNOW(cl)

modLoop<-foreach(i=1:2,.errorhandling="pass") %do% {
  require(picante)
  require(reshape)
  require(ggplot2)
  out<-modRun(3,nsites=10,compscale=1)
  return(list(out))}

names(modLoop)<-1:length(modLoop)
moddfb<-melt(modLoop,id.vars=colnames(modLoop[[1]][[1]]))

modIIIloop<-ggplot(moddfb,aes(x=Phylo.Relatedness,y=P_A,col=L1)) + geom_smooth(method="glm",family="binomial",formula=(y~(poly(x,2))),se=FALSE) + geom_smooth(method="glm",family="binomial",formula=(y~(poly(x,2))),se=FALSE,aes(group=1),linetype="dashed",col="black",size=.9) + labs("Iteration")
modIIIloop

p<-multiplot(modIIloop,modIIIloop)
p

#combine runs
dat<-list(moddf,moddfb)
names(dat)<-c("Model_II","Model_III")

dat<-melt(dat,id.var=colnames(moddf))

Bothmod<-ggplot(dat,aes(x=Phylo.Relatedness,y=P_A,col=L1)) + geom_smooth(method="glm",family="binomial",formula=(y~(poly(x,2))),se=TRUE) 
Bothmod

save.image(paste(gitpath,"simulation.RData",sep=""))

