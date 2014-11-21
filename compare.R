library(dplyr)

droppath<-"C:/Users/Ben/Dropbox/"

#Load image if desired
load(paste(droppath,"Thesis/Pred_Realized/Assemblages/Threshold0.2/Results/PredictedRealized.RData",sep=""))
new<-PA_m2

load("C:/Users/Ben/Desktop/Run.Rdata")
old<-PA_m2
str(old)
str(new)

old<-old[!old$Tree %in% "Func",]
old<-old[,!colnames(old) %in% c("Assemblage","Metric","Tree")]
colnames(old)[3]<-"Phylo.Relatedness"

#test species
ctor.old<-old[old$Species %in% "Coeligena.torquata",]
ctor.new<-new[new$Species %in% "Coeligena.torquata",]

ctor.old$A<-"Old"
ctor.new$A<-"New"

d<-rbind_all(list(ctor.old,ctor.new))

ggplot(d,aes(x=Locality,y=Phylo.Relatedness,col=A,shape=factor(P_A))) + geom_point(size=5) + geom_line(aes(group=Locality)) + theme_bw() + facet_wrap(~Hyp)

ggplot(d,aes(x=Phylo.Relatedness,fill=A)) + geom_histogram(position="dodge") + facet_wrap(~Hyp)

ggplot(d,aes(x=Locality,y=Suitability,fill=A)) + geom_bar(stat="identity",position="dodge") + theme_bw() + facet_wrap(~Hyp)

table(d$Phylo.Relatedness,d$A)

#something still about the fit?


#split data into types of assemblages
#drop Env + Dispersal for observed, its just a data subset
sdat<-split(new,list(new$P_R,new$Hyp),drop=TRUE)

#for now grab 10000 random samples and remove inf values
sdat<-lapply(sdat,function(x){
  x<-x[is.finite(x$Phylo.Relatedness),]
  return(x)
})


plotlist<-lapply(sdat,function(x){
  p<-ggplot(x,aes(x=Phylo.Relatedness,y=P_A)) + geom_smooth(method="glm",family="binomial",formula=y~poly(x,2)) + theme_bw() + labs("x=Distance to closest related species in an assemblage",y="Probability of presence") + ggtitle("")
  return(p)
})

#title plots
for (x in 1:length(plotlist)){
  plotlist[[x]]<-plotlist[[x]]+ggtitle(names(plotlist)[x])
}

do.call(grid.arrange,plotlist)


#get some comparable assemblages
pold<-old[old$P_R=="P_Apred" & old$Hyp %in% "Env",]
pnew<-new[new$P_R=="P_Apred" & new$Hyp %in% "Env",]

dim(pold)
dim(pnew)

head(pold)
head(pnew)

#get one assemblage
newl<-pnew[pnew$Locality %in% 10172,]
oldl<-pold[pold$Locality %in% 10172,]

oldl$Species[!oldl$Species %in% newl$Species]

head(oldl)
head(newl)

oldl$T<-"Old"
newl$T<-"New"

oldl
newl

a<-oldl[oldl$P_A==1,]$Species 
b<-newl[newl$P_A==1,]$Species
a[!a%in% b]
