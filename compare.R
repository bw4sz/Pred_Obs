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

ggplot(d,aes(x=Phylo.Relatedness,fill=A)) + geom_histogram(position="dodge")

ggplot(d,aes(x=Locality,y=Suitability,fill=A)) + geom_bar(stat="identity",position="dodge") + theme_bw() + facet_wrap(~Hyp)

table(d$Phylo.Relatedness,d$A)
