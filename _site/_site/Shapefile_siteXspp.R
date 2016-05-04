#Goal to take underlying distribution maps and create a siteXspp table and spatial information
#This is a source file called from the SpeciesOverlap Script

require(rgdal)
require(raster)
require(maptools)
require(stringr)
require(reshape2)

#Read in range maps
fil<-list.files(pattern="_pl.shp$","C:/Users/Ben/Dropbox/Lab paper 1 Predicted vs observed assemblages/Trochilidae/",full.names=TRUE)
filnam<-list.files(pattern="_pl.shp$","C:/Users/Ben/Dropbox/Lab paper 1 Predicted vs observed assemblages/Trochilidae/")

#Bring in niche models from the script SDM.R
setwd("D:\\Niche_Models")
setwd(paste(getwd(),cell_size,sep="/"))

#Biomod Consensus ensemble niche models
niche<-list.files(getwd(),pattern="TotalConsensus_EMbyROC.gri",full.name=T,recursive=T)

#Just from the current env predictions.
niche<-niche[grep("current",niche,value=FALSE)]

#Create blank raster, set extent in the future
r<-raster(niche[[1]])

#set extent if desired.

#Turn each shapefile into a raster map
species_rasterize<-list()

#Loop through all files
for (x in 1:length(fil)){
#load in file
  try(sp_shp<-readShapePoly(fil[[x]],delete_null_obj=TRUE))

#turn to raster
sp_ras<-rasterize(sp_shp,r)

#Which cells are 1's
presence_cells<-Which(sp_ras,cell=TRUE)

#check extent, return NA if no presences
if(length(presence_cells)==0) {
  species_rasterize[[x]]<-(NA); next}

#get xy of presence cells
species_rasterize[[x]]<-presence_cells
}

names(species_rasterize)<-strsplit(filnam,"_pl.shp")

#remove species with just NA (length==1)
species_full<-species_rasterize[!sapply(species_rasterize,length)==1]

#yields  how many species?
length(species_full)

#melt list into a dataframe
m.species<-melt(species_full)

siteXspp.rangemaps<-as.data.frame.array(table(m.species$L1,m.species$value))

####bring in Clade data
clades<-read.csv("C:\\Users\\Ben\\Dropbox//Shared Ben and Catherine//DimDivEntire//Files for Analysis//CladeList.txt",header=FALSE)[,-1]
colnames(clades)<-c("Clade","Genus","Species","double","English")

#take the first four letters of the species and genus to match the datatable
clades$abr<-tolower(paste(str_match(clades$Genus,"^.{0,4}"),str_match(clades$Species,"^.{0,4}"),sep="_"))

#fix rownames
rownames(siteXspp.rangemaps)<-sapply(rownames(siteXspp.rangemaps),function(x){
  if(sum(clades$abr %in% x)==0) return(x)
  new_name<-as.character(clades[clades$abr %in% x,"double"])
  gsub(" ",".",new_name)
  })

colnames(siteXspp.rangemaps)<-as.numeric(colnames(siteXspp.rangemaps))
#the spatial information for those cells which have species
sp.cells<-xyFromCell(r,as.numeric(colnames(siteXspp.rangemaps)))
sp.cells<-data.frame(sp.cells,cellnumber=as.numeric(colnames(siteXspp.rangemaps)))

#write to file
write.csv(siteXspp.rangemaps,"C:/Users/Ben/Dropbox/Thesis/Pred_Realized/Range_Maps//siteXspp_rangemaps.csv")
write.csv(sp.cells,"C:/Users/Ben/Dropbox/Thesis/Pred_Realized/Range_Maps//CellSites_Rangemaps.csv")


