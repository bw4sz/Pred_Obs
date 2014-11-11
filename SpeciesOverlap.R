##Predicted versus Realized Communities
#Begun August 16th 2012 - Ben Weinstein. Stony Brook University

#Load required libraries
require(raster)
require(reshape2)
require(ggplot2)
require(picante)
require(dismo)
require(ape)
require(doSNOW)
require(gdistance)
require(stringr)
require(foreach)
require(boot)
require(mapproj)
require(maptools)
require(rasterVis)
require(lmtest)

#Set dropbox path
droppath<-"C:/Users/Ben/Dropbox/"

#Set git path
gitpath<-"C:/Users/Ben/Documents/Pred_Obs/"

#Load image if desired
load(paste(droppath,"Thesis/Pred_Realized/Assemblages/Threshold0.2/Results/Run.RData",sep=""))

#################################################
#Source Function Script
source(paste(gitpath,"SpeciesOverlapSourceFunctions.R",sep=""))
#################################################

##############################################################################
#Niche Models Need to Be Run, see SDM.R script, or source the below function, careful this will take several hours
cell_size=0.1
source(paste(gitpath,"SDMupdated.R",sep=""))
#SDM_SP(cell_size)

#A SiteXspp table needs to be made at the same cell resolution
#This function will automatically look up the above niche models and create a observed matrix of 
#presence absence, based on the range maps, at the desired extent and resolution
#source("C:/Users/Ben/Documents/Pred_Obs/Shapefile_siteXspp.R")

##############################################################################

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

#Bring in Phylogenetic Data
trx<-read.nexus(paste(gitpath,"InputData\\ColombiaPhylogenyUM.tre",sep=""))
spnames<-read.table(paste(gitpath,"InputData\\SpNameTree.txt",sep="") , sep = "\t", header = TRUE)

#Replace tip.label with Spnames#
#replace the tiplabels with periods, which is the biomod default
trx$tip.label<-gsub("_",".",as.character(spnames$SpName))
ctrx<-cophenetic(trx)

####bring in Clade data
clades<-read.csv(paste(gitpath,"InputData//CladeList.txt",sep=""),header=FALSE)[,-1]
colnames(clades)<-c("Clade","Genus","Species","double","English")
clades<-clades[,1:5]

#Change the syntax on clades so it matches output, replace . with space
clades$double.<-sub(" ",".",clades$double)

#Bring in spatial data for the sites
#Extract env information for each site
Sites.sp<-SpatialPointsDataFrame(Sites[,c("LongDecDeg","LatDecDeg")],Sites)

################################################
#If we are using assemblage lists as importants start here
#################################################

siteXspp<-t(read.csv(paste(gitpath,"InputData//siteXspp.csv",sep=""),row.names=1))
#If you are using assemblages, set this flag, it will help analysis below
#method<-"Assemblage"

#Bring in niche models from the script SDM.R, get the folder from cell size arguemnt
setwd("D:\\Niche_Models")
setwd(paste(getwd(),cell_size,sep="/"))

#Biomod Consensus ensemble niche models
niche<-list.files(getwd(),pattern="ensemble.gri",full.name=T,recursive=T)

#Just from the current env predictions.
niche<-niche[grep("current",niche,value=FALSE)]

#Name the list of suitability models
names(niche)<-lapply(niche,function(path){
  split.1<-strsplit(path,"/")[[1]][4]
})

###Bring in trait data
morph <- read.csv(paste(gitpath,"InputData//MorphologyShort.csv",sep=""),na.strings="9999")

#just get males
morph.male<-morph[morph$Sex=="Macho",c("SpID","ExpC","Peso","AlCdo")]
morph.complete<-morph.male[complete.cases(morph.male),]

#aggregate for species
agg.morph<-aggregate(morph.complete,list(morph.complete$SpID),mean)
mon<-agg.morph[,-2]
colnames(mon)<-c("Species","Bill","Mass","WingChord")
rownames(mon)<-gsub(" ",".",mon[,1])
mon<-mon[,-1]

#principal component traits and get euclidean distance matrix
means <- apply(mon, 2, mean)

Bill <- mon$Bill - means["Bill"]/sd(mon$Bill)
Mass <- mon$Mass - means["Mass"]/sd(mon$Mass)
WingChord <- (mon$WingChord - means["WingChord"])/sd(mon$WingChord)

z.scores <- data.frame(Bill, Mass, WingChord)
rownames(z.scores) <- rownames(mon)

trait_pc <- as.matrix(dist(z.scores, method = "euclidean"))

#################################################
#Q1 Do we want to subset the models for number of presences in the assemblage lists?
#################################################
#Only retrieve niche models that have a reasonable number of presences in the assemblage lists
#What does the threshold look like?
#plot(sapply(1:100,function(x){
 # sum(apply(siteXspp,1,sum) > x)
#}))

#Not sure if i want to threshold it, if we want to make predicted assemblages, it isn't relevant if we have assemblage information, as long as suitability models performed well
curr_sp<-names(which(apply(siteXspp,1,sum) > 0))
#################################################

#Get the niche models with more than 0 presence points in the assemblage lists

niche_torun<-niche[names(niche) %in% curr_sp]

#Create Spatial Points object of the rownmaes
site.raster<-raster(niche[[1]])
res(site.raster)<-cell_size

#Create PA matrix for the raster
raster.localities<-rasterize(y=site.raster,SpatialPoints(Sites.sp),fun=function(x,...) length(x))

#Which cell is each site in?
cellSites<-extract(raster.localities,Sites.sp,cell=T)

#Split cell
head(cellSites<-data.frame(Sites.sp,cellSites))
splitCellSites<-split(cellSites,factor(cellSites$cells))

#How many duplicate communities are there, ie. number of assemblages per cell
table(sapply(splitCellSites,nrow))

siteXspp.raster<-sapply(splitCellSites,function(x){
  if(nrow(x)== 1) out<-siteXspp[,colnames(siteXspp) %in% x$Community]
  if(nrow(x)>1 ) {
    (out<-apply(siteXspp[,colnames(siteXspp) %in% x$Community],1,sum) > 0)*1
  }
  return(out)
})

#Get the xy lat long of the cells with sites in it
cellSitesXY<-xyFromCell(raster.localities,as.numeric(colnames(siteXspp.raster)),spatial=TRUE)
rownames(cellSitesXY@coords)<-colnames(siteXspp.raster)

##Important to get a grasp on data overlap between the different sources
congruence<-melt(list(Trait=rownames(trait_pc),Phylo=rownames(ctrx),Assemblage=rownames(siteXspp.raster),Suitability=names(niche),Clades=clades$double.))
write.csv(table(congruence),paste(droppath,"Thesis/Pred_Realized/Assemblages//DataOverlap.csv",sep=""))

print(colSums(table(congruence)))

#################################################
#Cost Path Analysis
#################################################
# #Import Friction layer, this can be changed later if we want a more fine grained 
# elevr<-raster(paste(droppath,"Thesis/Pred_Realized/etopo2\\w001001.adf",sep=""))
# 
# #cut it generally by the extent of the points, less waste
# elev.c<-crop(elevr,extent(raster(niche[[1]]))*1.3)
# 
# #only want above sealevel....
# elev.c[elev.c < 0]<-NA
# 
# #For now aggregate elev raster to reduce complexity
# elev.ca<-aggregate(elev.c,2)
# 
# #Find shortest cost path, for all sites, to dataheavy to bring into raster
# #This is find for assemblage mode, but for range map mode, skip down to after euclidean
# cl<-makeCluster(8,"SOCK")
# registerDoSNOW(cl)
# costPath.list<-foreach(x = 1:length(cellSitesXY)) %dopar% {
#   require(raster)
#   require(gdistance)
#   #pick the original site. 
#   orig<-cellSitesXY[x,]
#   
#   #What elevation is the origin
#   elev_origin<-extract(elev.ca,orig)[[1]]
#   if(is.na(elev_origin)) elev_origin<-0
#   
#   #Get the difference between the origin elevation and every cell in the raster
#   elev_diff<-abs(elev_origin-elev.ca)
#   
#   #create a the transition matrix where permeability is high where elev difference is low
#   trans_diff<-transition(elev_diff,function(x) 1/mean(x),8)
#   
#   #Diagonal Cell Correction, diagonals are father away than adjacent cells. 
#   slope <- geoCorrection(trans_diff)
#   
#   #Remember this cost surface is only valid for this site as the origen, ie. we want to create just this column in the cost distance matrix
#   #Cost Distance
#   cdist<-costDistance(slope,orig,cellSitesXY)
#   #labelling is key here.
#   
#   return(list(cdist))}
# stopCluster(cl)
# 
# #Format the cost path matrix
# m.costlist<-melt(costPath.list)
# CostPathMatrix<-log(cast(m.costlist,X2~L1)[,-1])
# rownames(CostPathMatrix)<-colnames(siteXspp.raster)
# colnames(CostPathMatrix)<-colnames(siteXspp.raster)

###############################################################
#Data Retrieval Complete, Begin Analysis
###############################################################

#sink overnight runs if needed
sink(paste(droppath,"Thesis/Pred_Realized/Overnight.txt",sep=""))

#######################################################
#Begin Analysis on Species Presence and Relatedness
#######################################################

#Create list of assemblages
sp.lists<-apply(siteXspp.raster,2,function(x){
  out<-names(x[which(x==1)])
})

######################################################
#Dispersal Limits
######################################################

#Get Costpath distribution for each species
costThresh<-apply(siteXspp.raster,1,function(y){
  sites<-names(y[y==1])
  #get pairwise combination of sites
  siteCombo<-expand.grid(sites,sites)
  outH<-apply(siteCombo,1,function(x){
    CostPathMatrix[x["Var1"],x["Var2"]]
  })
  #remove infinite values
  outH<-outH[is.finite(outH)]
  #return largest distance for the species, this will be the cost path threshold
  return(max(outH))})

hist(costThresh)

#sink overnight runs if needed
sink(paste(droppath,"Thesis/Pred_Realized/Overnight.txt",sep=""))

out_container<-list()
ord<-seq(0,.25,.05)
for (thresh in ord){
print(thresh)

fold<-paste(droppath,paste("Thesis/Pred_Realized/Assemblages/Threshold",thresh,sep=""),sep="")

#Create output folders
#Read in source function to draw predicted and realized assemblages and match relatedness

dir.create(fold)
dir.create(paste(fold,"Species",sep="/"))

#load models to file, makes it more transferable than keeping them on disk

#Choose the model output, for now i'm getting the mean ensemble model
modlist<-lapply(niche,raster,values=TRUE)

#or read from file
#writeRaster(paste(gitpath,"modelStack",sep=""),x=stack(modlist),bylayer=FALSE)
modlist<-stack(paste(gitpath,"modelStack",sep=""))

##########################################################
#Perform Predicted v Realized Function on all Niche Models
##########################################################

#Please note i've turned off the filter for which species to run based on number of presences in an assemblages, see line 105
cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)
system.time(all_models<-foreach(x=1:nlayers(modlist),.export=c("sp.lists","siteXspp.raster"),.errorhandling="pass",.packages=c("raster","maptools","reshape","dismo","picante","ggplot2","rasterVis","lattice")) %dopar%{  
  pred_realized(mod=modlist[[x]],thresh.suit=thresh,dispersal=FALSE,plots=FALSE)
})
stopCluster(cl)

print("SpeciesRelatednessandPresenceTable_Nodispersal")

#If you want to start here, read in the data from file

#combine all datasets
#This naming function would need to changed on other systems
names(all_models)<-names(modlist)

#Remove models that failed. 
working_models<-all_models[!sapply(sapply(all_models,nrow),is.null)]
failed_models<-all_models[sapply(sapply(all_models,nrow),is.null)]
write.csv(names(failed_models),"FailedSpecies.csv")

#melt and name list
#remove any species that failed. 
all.species.data<-melt(working_models,id.var=names(all_models[[1]]))

##############Repeat for no dispersal conditions##############

#Please note i've turned off the filter for which species to run based on number of presences in an assemblages, see line 105
cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)
system.time(all_modelsDD<-foreach(x=1:nlayers(modlist),.errorhandling="pass",.export=c("sp.lists","siteXspp.raster"),.packages=c("raster","maptools","reshape","dismo","picante","ggplot2","rasterVis","lattice")) %dopar%{  
  pred_realized(mod=modlist[[x]],thresh.suit=thresh,dispersal=TRUE,plots=FALSE)
})
stopCluster(cl)

print("SpeciesRelatednessandPresenceTable_dispersal")

#If you want to start here, read in the data from file

#combine all datasets
#This naming function would need to changed on other systems
names(all_modelsDD)<-names(modlist)

#Remove models that failed. 
working_modelsD<-all_modelsDD[!sapply(sapply(all_modelsDD,nrow),is.null)]
failed_modelsD<-all_modelsDD[sapply(sapply(all_modelsDD,nrow),is.null)]
write.csv(names(failed_modelsD),"FailedSpeciesD.csv")

#melt and name list
#remove any species that failed. 
all.species.dataD<-melt(working_modelsD,id.var=names(all_modelsDD[[1]]))

#Merge together as a list, name and melt

all.species.dataBoth<-list(Env=all.species.data,Env_Dispersal=all.species.dataD)
all.species.data<-melt(all.species.dataBoth,id.vars=colnames(all.species.data))

colnames(all.species.data)[13]<-"Hyp"
####################################################################
#Combine Models for each species and create data structure ready for plotting
#####################################################################

#Set working directory to output folder
dir.create(paste(fold,"Results",sep="/"))
setwd(paste(fold,"Results",sep="/"))

print(head(all.species.data))

#Visualize model outputs
ggplot(all.species.data,aes(x=Suitability,fill=P_A)) + geom_density(alpha=.5) + theme_bw() + labs(fill="") + facet_wrap(~Hyp)
ggsave("AllspeciesBinarySuitability.svg",height=8,width=10) 

#Name metrics column
colnames(all.species.data)[c(5,6)]<-c("Metric","value")
all.species.data<-all.species.data[,!colnames(all.species.data) %in% c("L1","P_A")]

#add in which clade each focal species is
PA_mult2<-merge(all.species.data,clades[,-c(3,4,5)],by.x="Species",by.y="double.")

#contigency table
contin<-table(PA_mult2$P_Apred,PA_mult2$PA_Binary,PA_mult2$Hyp,PA_mult2$Tree,deparse.level=2)

#melt the presence absence matrix for plotting
PA_m2<-melt(PA_mult2,measure.vars=c("P_Apred","PA_Binary"))

#Name the columns
colnames(PA_m2)[c(13,14)]<-c("P_R","P_A")


#distribution of trait distances.

with(PA_m2,hist(PA_m2[Tree=="Func","value"]))

PA_m2<-PA_m2[is.finite(PA_m2$value),]


#remove all na's
PA_m2<-PA_m2[!is.na(PA_m2$value),]
#Final data format

PA_m2<-PA_m2[!(PA_m2$Hyp=="Env_Dispersal" & PA_m2$P_R=="PA_Binary"),]

#Final data format
head(PA_m2)

#take out values greater than 95% quartile
trait_out<-with(PA_m2,quantile(PA_m2[Tree=="Func","value"],.95,na.rm=T))
phylo_out<-with(PA_m2,quantile(PA_m2[Tree=="Phylo","value"],.95,na.rm=T))

PA_m2[PA_m2$Tree=="Func" & PA_m2$value > as.numeric(trait_out),"value"]<-NA 
PA_m2[PA_m2$Tree=="Phylo" & PA_m2$value > as.numeric(phylo_out),"value"]<-NA 

PA_m2<-PA_m2[!is.na(PA_m2$value),]


###################################################################
##########Data Visualization: What is the relationship between presence and relatedness?
###################################################################

#run results and model function
#mods<-LogT(PA_m2)

########Just keep species with a reasonable number of presence points
pflat<-dcast(PA_m2,...~Hyp+Tree+P_R,value.var="P_A")

#species with atleast 10 presences
l<-aggregate(pflat$Env_Phylo_PA_Binary,list(pflat$Species),sum,na.rm=TRUE)

keep<-l[l$x > 15,]

#replot with just common species
commsp<-PA_m2[PA_m2$Species %in% keep$Group.1,]

#refactor species
commsp$Species<-droplevels(as.factor(commsp$Species))

#remove florisuga for the moment, bad data?
commsp<-commsp[!commsp$Species %in% "Florisuga.mellivora",]

commsp<-PA_m2[PA_m2$Species %in% keep$Group.1,]

#refactor species
commsp$Species<-as.factor(commsp$Species)

#Change the word Func to Trait
commsp$Tree<-as.character(commsp$Tree)
commsp[commsp$Tree %in% "Func","Tree"]<-"Trait"

############Model FIT##########

#compute psuedo R^2

#compute psuedo R
sapply(mods,function(x) (x$null.deviance - x$deviance)/x$null.deviance)
sapply(mods_nopoly,function(x) (x$null.deviance - x$deviance)/x$null.deviance)


#
##Looking at the logistic modeling output, not all models should be binomial
p<-ggplot(commsp[!(commsp$Tree %in% "Phylo" & commsp$P_R %in% "PA_Binary") & !(commsp$Tree %in% "Func" & commsp$Hyp %in% "Env" & commsp$P_R %in% "P_Apred"),],aes(x=value,y=as.numeric(P_A),col=P_R,linetype=Hyp)) + geom_smooth(method="glm",family="binomial") 
p<- p + labs(col="Assemblage") + scale_color_discrete(labels=c("Predicted","Observed")) 
p<-p+theme_bw() + xlab("") + ylab("") + scale_y_continuous(limits=c(0,1)) + labs(linetype="Predicted Assemblage",x="Distance to nearest species in assemblage",y="Probability of Occurrence")
p<-p + facet_grid(~Tree,scales="free_x") + geom_smooth(data=commsp[(commsp$Tree %in% "Phylo" & commsp$P_R %in% "PA_Binary") | (commsp$Tree %in% "Func" & commsp$Hyp %in% "Env" & commsp$P_R %in% "P_Apred"),],family="binomial",formula=(y~poly(x,2)),method="glm")
print(p) 
ggsave(plot=p,"ModelFits_linklihoodtest.svg",height=6,width=10,dpi=300)

#With panel bins
ph<-ggplot(commsp[!(commsp$Tree %in% "Phylo" & commsp$P_R %in% "PA_Binary") & !(commsp$Tree %in% "Func" & commsp$Hyp %in% "Env" & commsp$P_R %in% "P_Apred" ) ,],aes(x=value,y=as.numeric(P_A),col=P_R,linetype=Hyp)) 
ph<-ph + geom_smooth(data=commsp[(commsp$Tree %in% "Phylo" & commsp$P_R %in% "PA_Binary") | (commsp$Tree %in% "Func" & commsp$Hyp %in% "Env" & commsp$P_R %in% "P_Apred"),],family="binomial",formula=(y~poly(x,2)),method="glm")
ph<-ph + facet_wrap(~Tree,scales="free_x")
ph<-ph +geom_bin2d(data=commsp[commsp$P_R=="PA_Binary",],aes(x=value,y=as.numeric(P_A),fill= ..density..*100),col="black",linetype=1)+ geom_smooth(method="glm",family="binomial") 
ph <- ph + scale_fill_continuous(low="grey90",high="black")
ph<- ph + labs(col="Assemblage") + scale_color_discrete(labels=c("Predicted","Observed")) 
ph<-ph+theme_bw() + scale_y_continuous(limits=c(0,1)) +  labs(linetype="Predicted Assemblage Type",x="Distance to nearest species in assemblage",y="Probability of Occurrence",fill="Point Density (%)")
print(ph) 
ggsave(plot=ph,"Panel1bins.svg",height=6,width=10,dpi=300)
ggsave(plot=ph,"Panel1bins.jpeg",height=6,width=10,dpi=300)

##########################
#Species Plots
##########################

#only plot species with converging models
#splitD<-split(PA_m2,list(PA_m2$P_R,PA_m2$Tree,PA_m2$Hyp,PA_m2$Species),drop=TRUE)

# mod.table<-lapply(splitD,function(x){
#   mod<-glm(data=x,as.numeric(P_A)~value,family="binomial")
#   modpoly<-glm(data=x,as.numeric(P_A)~poly(value,2),family="binomial") 
#   l<-lrtest(mod,modpoly)
#   ltp<-round(l[["Pr(>Chisq)"]][2],2)
#   pmod<-round(summary(mod)$coefficients[2,4],2)
# 
#   pmodpoly<-round(summary(modpoly)$coefficients[3,4])
#   data.frame(Species=unique(mod$data$Species),Hyp=unique(mod$data$Hyp),Tree=unique(mod$data$Tree),P_R=unique(mod$data$P_R),ratioTest=ltp,pmod,pmodpoly)
# })
# 
# modt<-rbind.fill(mod.table)
# sapply(mod.table,function(x){
#   sum(x$P < .05,na.rm=TRUE)
# })
# 
# #which species should be discarded through lack of model fit
# tocast<-head(melt(modt,measure.vars=c("pmod","pmodpoly")))
# modPR<-dcast(tocast,...~P_R+variable)
# 
# require(plyr)
# modt<-rbind.fill(mod.table)
# 
# ###How many species to include in the model, needs to have significant estimates for both predicted and observed 
# remove<-modt[modt$Tree %in% "Phylo"& (modt$pmod > 0.05 | modt$pmodpoly > 0.05),"Species"]

#replot with just common species

##########Species level data, screen for presences

#Easiest to visualize if pre split into phylo and trait
phylo_sp<-commsp[commsp$Tree %in% "Phylo",]
trait_sp<-commsp[commsp$Tree %in% "Trait",]

#Phylo With panel bins
ph<-ggplot(phylo_sp,aes(x=value,y=as.numeric(P_A),col=P_R,linetype=Hyp)) 
ph<-ph + facet_wrap(~Clade,scales="free")
ph<-ph + geom_smooth(family="binomial",formula=(y~poly(x,2)),method="glm")
ph <- ph+geom_bin2d(data=phylo_sp[phylo_sp$P_R=="PA_Binary",],aes(x=value,y=as.numeric(P_A),fill= ..density..*100),col="black",linetype=1)
ph<-ph+ scale_fill_continuous(low="grey90",high="black")
ph<- ph + labs(col="Assemblage") 
ph<-ph+theme_bw() + scale_y_continuous(limits=c(0,1)) +  labs(linetype="Predicted Assemblage",x="Distance to nearest species in assemblage",y="Probability of Occurrence",fill="Point Density (%)")
print(ph) 

<<<<<<< HEAD
=======
ph<-ggplot(trait_sp,aes(x=value,y=as.numeric(P_A),col=P_R,linetype=Hyp)) 
ph<-ph + facet_wrap(~Clade,scales="free")
ph<-ph + geom_smooth(family="binomial",formula=(y~poly(x,2)),method="glm")
ph <- ph+geom_bin2d(data=trait_sp[trait_sp$P_R=="PA_Binary",],aes(x=value,y=as.numeric(P_A),fill= ..density..*100),col="black",linetype=1)
ph<-ph+ scale_fill_continuous(low="grey90",high="black")
ph<- ph + labs(col="Assemblage") 
ph<-ph+theme_bw() + scale_y_continuous(limits=c(0,1)) +  labs(linetype="Predicted Assemblage",x="Distance to nearest species in assemblage",y="Probability of Occurrence",fill="Point Density (%)")
print(ph) 


###################################################################
##########Data Visualization: What is the relationship between presence and relatedness?
###################################################################

############MODEL Fit

require(lme4)

splitD<-split(PA_m2,list(PA_m2$P_R,PA_m2$Tree,PA_m2$Hyp),drop=TRUE)
lmm<-glmer(data=splitD[[1]],family="binomial",P_A~value +(1|Species))

summary(lmm)

names(mods)<-names(splitD)

R2<-lapply(mods,function(x){
  1-x$deviance/x$null.deviance
})

###########Model Fitting of suitability alone###################
#split data into combinations of Assemblage and Tree
#Comapre glm models with polynomials, loglik ratio test, is p < 0.05 we can reject the null that the simpler model explains the data bette


mods_nopoly<-lapply(splitD,function(x){
  mod_poly2<-glm(data=x,as.numeric(P_A)~value,family="binomial") 
})
names(mods_nopoly)<-names(splitD)

lapply(mods_nopoly,summary)
summary(mods_nopoly)

>>>>>>> 5ca833f9eea769acd2241a40a2a5fddad75d8c6c

#write to file
save.image("Run.Rdata")
out_container[[which(ord %in% thresh)]]<-contin
}


###species richness in predicted and observed assemblages

split.D<-split(all.species.data,list(all.species.data$Tree,all.species.data$Hyp))
x<-split.D[[1]]

#Predicted siteXspp
siteXspp.pred<-dcast(x,Species~Locality,value.var="P_Apred")

#set rownames
rownames(siteXspp.pred)<-siteXspp.pred$Species
siteXspp.pred<-siteXspp.pred[,-1]

richnessPred<-apply(siteXspp.pred,2,sum,na.rm=TRUE)

#Observed siteXspp
siteXspp.obs<-dcast(x,Species~Locality,value.var="PA_Binary")

#set rownames
rownames(siteXspp.obs)<-siteXspp.obs$Species
siteXspp.obs<-siteXspp.obs[,-1]

richnessObs<-apply(siteXspp.obs,2,sum,na.rm=TRUE)

#plot differences
ggplot() + geom_histogram(aes(richnessPred),fill="red",alpha=.6) + geom_histogram(aes(richnessObs),fill="blue",alpha=.6)
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/RichnessCompare.jpeg",sep=""),height=8,width=10,dpi=300)

##################################################
#How does the threshold effect the relatedness and presence relationship?
##################################################

#read from file outputs

ord<-seq(0,.25,.05)
out_container<-list()
for (thresh in ord){
  print(thresh)
  
  fold<-paste(droppath,paste("Thesis/Pred_Realized/Assemblages/Threshold",thresh,sep=""),sep="")

  setwd(paste(fold,"Results",sep="/"))
  
  load("Run.Rdata")
  index<-which(ord %in% thresh)
  print(index)
  out_container[[index]]<-contin
}
out_container
names(out_container)<-ord


###Calculate sensitivity score for each 
melt_out<-melt(out_container)

colnames(melt_out)<-c("Predicted","Observed","Hyp","Tree","value","threshold")

#take the tree metrics out, they are idenitcal
melt_out<-melt_out[melt_out$Tree %in% "Phylo",]

#combine predictged and observed

melt_out$Combo<-paste(melt_out$Predicted,melt_out$Observed,sep=" ")
ggplot(melt_out,aes(x=as.numeric(threshold),y=value,col=Hyp)) + geom_line() + facet_wrap(~Combo,scales="free") + labs(col="Model") + xlab("Suitability Threshold") + ylab("Number of Sites")
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/ModelSensitivity.svg",sep=""),height=10,width=8,dpi=300)
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/ModelSensitivity.jpeg",sep=""),height=8,width=10,dpi=300)


##sensitivity score for each threshold and hypothesis

split.melt<-split(melt_out,list(melt_out$threshold,melt_out$Hyp))

#calculate True positive/ true positive + false negative -> sensitivity
#calcualte true negative / false negatice + true negative -> specificity

#all splits are the same row order. 

spliteval<-rbind.fill(lapply(split.melt,function(x){
  sensitivity<-x[4,"value"]/(x[4,"value"] + x[3,"value"])
  specificity<-(x[1,"value"]/(x[1,"value"] + x[2,"value"]))
  #create dataframe of outputs
  df.out<-data.frame(Threshold=unique(x$threshold),Hyp=unique(x$Hyp),sensitivity,specificity)
  }))

#plot sensitivity and specificity
ggplot(spliteval,aes(y=sensitivity,x=specificity,col=Threshold,shape=Hyp)) + geom_line(col="black",alpha=.7) + geom_point(size=4) + ylab("sensitivity") 
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/ModelMetrics.svg",sep=""),height=10,width=8,dpi=300)
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/ModelMetrics.jpeg",sep=""),height=8,width=10,dpi=300)


### i want to minimize the number of predicted absence and presences, and maximize the number of predicted 1,1

spliteval2<-rbind.fill(lapply(split.melt,function(x){
  truepresence<-x[4,"value"]/sum(x$value)*100
  falsepresence<-x[3,"value"]/sum(x$value)*100
  trueabsence<-x[1,"value"]/sum(x$value)*100
  falseabsence<-x[2,"value"]/sum(x$value)*100
  #create dataframe of outputs
  df.out<-data.frame(Threshold=unique(x$threshold),Hyp=unique(x$Hyp),truepresence,falsepresence,trueabsence,falseabsence)
}))


#plot sensitivity and specificity

mval<-melt(spliteval2)

ggplot(mval,aes(y=value,col=variable,x=Threshold,shape=Hyp)) + geom_line(col="black",alpha=.7) + geom_point(size=4) + facet_wrap(~Hyp) + ylab("True Positive")
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/ModelMetrics2.svg",sep=""),height=10,width=8,dpi=300)
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/ModelMetrics2.jpeg",sep=""),height=8,width=10,dpi=300)
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/ModelMetrics2.pdf",sep=""),useDingbats=FALSE,height=8,width=10,dpi=300)

ggplot(spliteval2,aes(x=trueabsence,y=truepresence,col=Threshold,shape=Hyp)) + geom_line(col="black",alpha=.7) + geom_point(size=4) + facet_wrap(~Hyp) + ylab("True Positive") + 
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/ModelMetrics3.jpeg",sep=""),height=8,width=10,dpi=300)

#plot distance to 1,1 which is the goal of the model fit.

spliteval$Dist_1<-as.matrix(dist(rbind(cbind(1,1),cbind(spliteval$sensitivity,spliteval$specificity))))[1,][-1]
ggplot(spliteval,aes(x=Threshold,y=Dist_1,shape=Hyp)) + geom_line(col="black",alpha=.7) + geom_point(size=4) + facet_wrap(~Hyp)


#Get the thresholds for each species
fil<-list.files(paste(droppath,"Thesis/Pred_Realized/Assemblages",sep=""),pattern="SuitThreshold.csv",recursive=TRUE,full.names=TRUE)
out_container<-lapply(fil,read.csv,row.names=1)
head(m.Thresh<-melt(out_container))
m.Thresh$variable<-as.numeric(str_match(fil,pattern="Assemblages/\\w+(.\\w+)")[,2])
m.Thresh$L1<-str_match(fil,pattern="Species/(\\w+.\\w+)")[,2]

#the NA's are 0's, they don't have decimals
m.Thresh$variable[is.na(m.Thresh$variable)]<-0

ggplot(m.Thresh,aes(x=as.factor(variable),y=value)) + geom_boxplot() + theme_bw() + labs(y="Suitability",x="Percentage of Localities") 
ggsave("SpeciesThresholds.svg",height=8,width=10)

ggsave("SpeciesThresholds.jpeg",height=8,width=10)

###threshold and results figure 

ord<-seq(0,.25,.05)
out_graph<-list()
for (thresh in ord){
  print(thresh)
  
  fold<-paste(droppath,paste("Thesis/Pred_Realized/Assemblages/Threshold",thresh,sep=""),sep="")
  
  setwd(paste(fold,"Results",sep="/"))
  
  load("Run.Rdata")
  
  #Final data format
  head(PA_m2)
  
  PA_m2<-PA_m2[is.finite(PA_m2$value),]
  
  #remove extreme species from morphological dataset
  extremesp<-which(PA_m2$Species %in% c("Ensifera.ensifera","Patagona.gigas") & PA_m2$Tree =="Func")
  
  #remove extreme species
  PA_m2<-PA_m2[-extremesp,]
  
  PA_m2<-PA_m2[!(PA_m2$Hyp=="Env_Dispersal" & PA_m2$P_R=="PA_Binary"),]
  index<-which(ord %in% thresh)
  print(index)
  out_graph[[index]]<-PA_m2
}

names(out_graph)<-ord

ot<-melt(out_graph,id.vars=colnames(PA_m2))

ph<-ggplot(ot,aes(x=value,y=as.numeric(P_A),col=P_R,linetype=Hyp)) + geom_smooth(method="glm",family="binomial",formula=(y~(poly(x,2)))) 
ph<- ph + labs(col="Assemblage") + scale_color_discrete(labels=c("Predicted","Observed")) 
ph<-ph+theme_bw() + xlab("") + ylab("") + scale_y_continuous(limits=c(0,1))
ph<-ph + facet_grid(L1~Tree,scales="free_x")
print(ph)
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/Thresh_Graph.pdf",sep=""),useDingbats=FALSE,height=10,width=7,dpi=300)
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/Thresh_Graph.svg",sep=""),height=10,width=7,dpi=300)
ggsave(paste(droppath,"Thesis/Pred_Realized/Assemblages/Thresh_Graph.jpeg",sep=""),height=10,width=7,dpi=300)

save.image(paste(droppath,"Thesis/Pred_Realized/Pred_Obs.RData",sep=""))

#Create output table of species, number of samples, model score, variable importance, and sensitivity threshold
m.Thresh<-m.Thresh[m.Thresh$variable %in% .2,]
colnames(m.Thresh)<-c("Threshold","Suitability Score","Species")

#match format of species column
m.Thresh$Species<-gsub("\\.", " ",m.Thresh$Species)

#combine with model_eval
mall<-dcast(model_eval,Species~variable+Metric)
merge1<-merge(m.Thresh,mall,by="Species")

#number of localities
nloc<-melt(rec)
colnames(nloc)<-c("Species","Number_Localities")

species_tab<-merge(merge1,nloc,by="Species")

write.csv(species_tab,"species_tab.csv")



sink()