y<-sdat[[x]][1:10000,]
out<-Bayes(presence=y$P_A,distance=y$Phylo.Relatedness,species=y$Species,runs = runsl[x],burn=((runsl[x]/25)-2000),Iteration=y$Iteration,hier=T)
pars<-melt(out$BUGSoutput$sims.array)

colnames(pars)<-c("Draw","Chain","parameter","estimate")


spec<-pars[str_detect(pars$parameter,"alpha|beta|gamma"), ]
parmelts<-pars[!str_detect(pars$parameter,"alpha|beta|gamma|deviance"), ]

ggplot(droplevels(parmelts),aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_wrap(~parameter,scales="free") + theme_bw() + labs(col="Chain")

gdf<-parmelts

trajsplit<-dcast(gdf,...~parameter,value.var="estimate")

#trajsplit<-gdcast[sample(1:nrow(gdcast),500,replace=F),]

#define trajectory function
trajF<-function(intercept,linear,polynomial,x){
  p<-inv.logit(intercept + linear * x  + polynomial * x^2)
  return(p)
}

specCast<-list()

for(x in 1:length(trajsplit)){
  a<-trajsplit[[x]]
  dat<-sdat[[which(names(sdat) %in% names(trajsplit)[x])]]
  
  est<-list()
  for(p in 1:nrow(a)){
    
    est[[p]]<-data.frame(Relatedness=dat$Phylo.Relatedness,y=trajF(a[p,]$intercept,a[p,]$linear,a[p,]$polynomial,dat$Phylo.Relatedness))
    
  }
  est<-rbind_all(est)
  mtraj<-group_by(est,Relatedness) %>% summarise(mean=mean(y),upper=quantile(y,0.975),lower=quantile(y,0.025))
  
  specCast[[x]]<-mtraj
}

names(specCast)<-names(trajsplit)

#name and reduce
for(x in 1:length(specCast)){
  specCast[[x]]$Assemblage<-names(specCast)[x]
}

specdf<-rbind_all(specCast)

#define another facet
#add in observed predicted and eventually simulated labels

specdf[specdf$Assemblage %in% c("Env","Env + Dispersal","Observed"),"Type"]<-"Observed and Predicted"

#overlay together
ggplot(data=specdf[specdf$Assemblage %in% c("Env","Env + Dispersal","Observed"),],aes(x=Relatedness,y=mean,ymin=lower,ymax=upper,fill=Assemblage)) + geom_ribbon(alpha=.5) + geom_line() + theme_bw()+ labs(x="Phylogenetic distance to closest related species",y="Probability of Occurrence") + scale_fill_manual(values=c("grey80","grey60","black"),labels=c("Environment","Environment + Dispersal","Observed"))  
