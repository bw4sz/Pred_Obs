#Load required libraries
library(reshape2)
require(ggplot2)
library(picante)
library(dismo)
library(ape)
library(doSNOW)
library(gdistance)
library(ggbiplot)
library(foreach)
library(boot)
library(maptools)
library(rasterVis)
library(knitr)
library(vegan)
library(gridExtra)
library(R2jags)
library(stringr)
library(dplyr)
library(scales)
library(ggbiplot)

opts_chunk$set(warning=FALSE,message=FALSE)
opts_chunk$set(cache=FALSE, fig.path='figure/',fig.width=10,echo=TRUE)

#Set git path
gitpath<-"C:/Users/Ben/Documents/Pred_Obs/"
droppath<-"C:/Users/Ben/Dropbox/"

setwd(gitpath)

#env=c('bio_1','bio_5','bio_12')
env<-paste("bio_",seq(1,19,1),sep="")
#Load image if desired
#load(paste(droppath,"Thesis\\Pred_Realized\\Assemblages\\Threshold0.05\\Results\\","PredictedRealized_Revis.RData",sep=""))

#Bring in Phylogenetic Data
trx<-read.tree(paste(gitpath,"InputData\\hum294.tre",sep=""))

#format tips
new<-str_extract(trx$tip.label,"(\\w+).(\\w+)")
#get duplicates
trx<-drop.tip(trx,trx$tip.label[duplicated(new)])

#name tips.
trx$tip.label<-str_extract(trx$tip.label,"(\\w+).(\\w+)")

ctrx<-cophenetic(trx)

ctrx<-ctrx/max(ctrx)

#Bring in the assemblages, cleaned from JP
Sites<-read.csv(paste(gitpath,"InputData//Sites.csv",sep=""),row.names=1)

####bring in Clade data
clades<-read.csv(paste(gitpath,"InputData//CladeList.txt",sep=""),header=FALSE)[,-1]
colnames(clades)<-c("Clade","Genus","Species","double","English")
clades<-clades[,1:5]

#Change the syntax on clades so it matches output, replace . with space
clades$double.<-sub(" ",".",clades$double)

#Bring in spatial data for the sites
#Extract env information for each site
Sites.sp<-SpatialPointsDataFrame(Sites[,c("LongDecDeg","LatDecDeg")],Sites)

#site by species lists
siteXspp<-t(read.csv(paste(gitpath,"InputData//siteXspp.csv",sep=""),row.names=1))

###Source Function Script

source(paste(gitpath,"SpeciesOverlapSourceFunctions.R",sep=""))
source(paste(gitpath,"SDMupdated.R",sep=""))

#Niche Models
cell_size=0.1
inLocalities<-read.csv("InputData/MASTER_POINTLOCALITYarcmap_review.csv")
envfolder<-"C:\\Users\\Ben\\Dropbox\\Thesis\\Pred_Realized\\EnvLayers\\Cropped\\"
savefolder<-"C:/Users/Ben/Dropbox/Thesis/Pred_Realized/NicheModels/vartest"

#If you have already installed, let's start here.
#Call the packages we are going to need in this tutorial
library(biomod2)
library(maptools)
library(raster)
library(rgdal)
library(stringr)
library(foreach)

#set a working directory, where do we want to save files
#save locally for now
setwd(savefolder)

dir.create(paste(getwd(),cell_size,sep="/"))
setwd(paste(getwd(),cell_size,sep="/"))
dir.create("logs")
#To perform the biomod, you must have three pieces of data
#1) Presence Absence Matrix
#2) Input localities in a 2 column matrix
#3) Environmental Variables - unclear whether this should be masked or not

#################################

#1) presence absence data matrix for the desired species
#For the graham lab, or users of the dropbox, just change the user below to your name, or the whatever the dropbox path is:"

#Lets go get the presence data on hummingbird distributions
PA<-inLocalities

#Just take the columns you want. 
PAdat<-PA[,colnames (PA) %in% c("RECORD_ID","SPECIES","COUNTRY","LOCALITY","LATDECDEG","LONGDECDEG","Decision","SpatialCheck","MapDecision")]

#There is one errant record.
PAdat<-PAdat[!PAdat$LONGDECDEG==-6,]

myExpl<-c()
#Import environmental data from worldclim, three variables

if(env=='all'){
  e<-list.files(envfolder,".gri")
  for (x in 1:length(e)){
    myExpl[[x]]<- paste(envfolder,e[[x]],sep="")
  }
} else {
  for(x in 1:length(env)){
    myExpl[[x]]<-paste(envfolder,env[[x]],".gri",sep="")
  }
}


myExpl<-stack(myExpl)

#Just get the clean localities
loc_clean<-PAdat[PAdat$SpatialCheck=="Y" & PAdat$MapDecision %in% levels(PAdat$MapDecision)[!levels(PAdat$MapDecision) %in% "REJECT"],]

#To be safe, take out the bogota localities, likely just bad musueum records
loc_clean<-loc_clean[!loc_clean$LOCALITY %in% "BOGOTA",]

extPoint<-SpatialPoints(loc_clean[,c("LONGDECDEG","LATDECDEG")])
#exte<-extent(c(-81.13411,-68.92061,-5.532386,11.94902))

#Crop by this layer, 
myExpl<-crop(myExpl,extent(extPoint)*1.2)
res(myExpl)

#Set Cell size
####################################
fact<-cell_size/res(myExpl)
####################################

#Set cell size
myExpl<-aggregate(myExpl,fact)

exte<-extent(c(-81.13411,-68.92061,-5.532386,11.94902))

#Cut by the extent

#Crop by this layer, 
myExpl.crop<-stack(crop(myExpl,exte))

#correlation
cort<-as.matrix(myExpl.crop)
cort<-cort[!is.na(cort[,1]),]

m<-melt(cor(cort[,c("bio_1","bio_4","bio_12","bio_15")]))
ggplot(data=m) + geom_tile(aes(x=X1,y=X2,fill=value<0.4))
ggbiplot(princomp(cort),scale=1)

#create a list of all env to project into
projEnv<-list(myExpl.crop)
names(projEnv)<-c("current")

#######################
#Which species to run
#######################

#How many records per species?
rec<-table(loc_clean$SPECIES)

#let's grab some species that have been checked and have more than 20 input localities 
spec<-names(rec[which(rec >= 20)])

###################################################
#remove any species that have already been run
###################################################

#name the list with the correct species names from file
#get all the niche model data
niche<-list.files(getwd(),pattern="ensemble.gri",full.name=T,recursive=T)

completed<-sapply(strsplit(str_extract(niche,pattern="(\\w+.\\w+)/proj_current"),"/"),function(x){
  x[1]
})

#Only run a species if it is NOT run in all models
spec<-spec[!spec %in% gsub("\\."," ",completed)]
paste("Species to be modeled",spec,sep=": ")

