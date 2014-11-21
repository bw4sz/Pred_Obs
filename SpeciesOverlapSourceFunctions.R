#Apply predicted v realized function to all species
#Define function
pred_realized<-function(mod,thresh.suit,dispersal,PA_phylo,plots=TRUE,loc_clean=loc_clean){
  
  setwd(paste(fold,"Species",sep="/"))
  
  names(mod)<-strsplit(names(mod[[1]]),"_")[[1]][1]
  
  #Create a directory to store outputs
  dir.create(names(mod))
  
  setwd(names(mod))
  
  #Get the number of records per site
  records_site<-apply(siteXspp.raster,1,sum,na.rm=TRUE)
  
  #Get row in siteXspp.raster for the model species
  rowS<-siteXspp.raster[rownames(siteXspp.raster) %in% names(mod),]
  
  #get siteXspp matrix just for the cells
  #for each species get the list of sites its present and absent
  PA<-list(names(rowS[which(rowS==1)]),names(rowS[which(rowS==0)]))
  names(PA)<-c("Present","Absent")
  
  #If species has not presence in our communities
  if(sum(sapply(PA,is.null)*1)==2){
    PA$Absent<-as.numeric(colnames(rowS))
  }
  
  #Sample raster with each the presence and absence points
  Com<-data.frame(cellSitesXY,colnames(siteXspp.raster))
  colnames(Com)<-c("Long","Lat","Cell")
  
  ######################################Step 1, what is the predicted probability at present v absent sites?
  ###Build a function that gets the presence and absence localities and samples raster at that location
  
  PA_niche<-lapply(PA,function(f){ 
    if(length(f)==0){return(data.frame(Locality=NA,Suitability=NA))}
    
    #Turn into a spatial point
    Com.f<-Com[Com$Cell %in% factor(f),]
    rownames(Com.f)<-Com.f$Cell
    species.points<-SpatialPoints(Com.f[,-3])
    
    #Get the values of the niche model
    niche.values<-extract(mod,species.points)
    out<-data.frame(f,niche.values)
    colnames(out)<-c("Locality","Suitability")
    return(out)})
  
  #melt and cast niche value data so that we have every combinition of species/localitiy/presence
  giant<-melt(PA_niche,id.var="Locality")
  species.data<-dcast(giant,...~variable)
  colnames(species.data)<-c("Locality","P_A","Suitability")
  
  #add in species name
  species.data$Species<-names(mod)
  
  #This isn't worth much, just showing that the niche models work well.
  if(plots==TRUE){
  giant.p<-ggplot(data=species.data,aes(x=P_A,y=Suitability)) + geom_boxplot() + xlab("") + theme_bw()+ggtitle(names(mod))

  filname<-paste(names(mod),"PA_suitability.jpeg",sep="")
  ggsave(plot=giant.p,filname,height=8,width=4)
  }
   
  #For the threshold column, what is the distribution of presence sites
  #match formating for species names
  sp.loc<-loc_clean[loc_clean$SPECIES %in%  gsub("\\."," ",as.character(unique(species.data$Species))),]
  sp.loc<-SpatialPointsDataFrame(sp.loc[,c("LONGDECDEG","LATDECDEG")],sp.loc)
  
  #draw suitability for occurence points
  site_suit<-extract(mod,sp.loc)
  suit_cut<-quantile(na.omit(site_suit),thresh.suit)
  
  #plot to file
  if(plots==TRUE){
  qplot(site_suit) + geom_vline(aes(xintercept=suit_cut),linetype="dashed",col="Red")
  ggsave("SuitThresholds.jpeg",height=8,width=11,dpi=300)
  }
  
  #Predicted Presence absence Column
  species.data$P_Apred<-(species.data$Suitability > suit_cut)*1
  
  #Turn Presence absence to 0-1 and do logistic regressions?
  species.data$PA_Binary<-(species.data$P_A=="Present")*1
  
  #Map distributions for presence and absence
  if(plots==TRUE){
    
  ggplot(species.data,aes(x=Suitability,fill=P_A)) + geom_density(alpha=.5) + geom_vline(xintercept=suit_cut,linetype="dashed") + theme_bw() + labs(fill="")
  ggsave("SpeciesDensity.svg",height=8,width=10)
  }
  
  ######################################################
  #####Cost Path Filter for predicted presence##########
  ######################################################
  if(dispersal){
    
    #get the cost path threshold
    cthres<-costThresh[names(mod)]
    
    #get the distance from every obs site to every pred site
    pred_p<-levels(species.data[species.data$P_Apred==1,"Locality"])
    obs_p<-levels(species.data[species.data$PA_Binary==1,"Locality"])
    
    #Min distance from any obs site to a pred site
    minSite<-sapply(pred_p,function(y){
      combs<-expand.grid(y,obs_p)
      #remove to cell to itself
      combs<-combs[!combs$Var2 %in% combs$Var1,]  
      outP<-apply(combs,1,function(x){
        CostPathMatrix[x["Var1"],x["Var2"]]
      })
      if(length(outP) < 2){return(min(outP[is.finite(outP)]))}
      return(min(outP[is.finite(outP)])
      )
    })
    
    #Which predicted sites are outside the threshold
    toRemove<-minSite[minSite > cthres]
    
    #Turn any predicted presence and absences that are outside the threshold
    species.data[species.data$Locality %in% names(toRemove),"P_Apred"]<-0
    }
  
  #################################################
  #Mapping of Predictor v Realized Assemblages
  #################################################
  
  #   #Get the XY coordinates of all the cells that have assemblage information
  xy<-data.frame(cellSitesXY@coords)
  colnames(xy)<-c("LongDecDeg","LatDecDeg")
  
  #Merge with species information
  species.data<-merge(species.data,xy,by.x="Locality",by.y="row.names")
  spat.species<-SpatialPointsDataFrame(species.data[,c("LongDecDeg","LatDecDeg")],species.data)
  
  if(plots==TRUE){
    
  writeSpatialShape("SpatialSpecies",x=spat.species)
  
  #Write threshold to file
  write.csv(suit_cut,"SuitThreshold.csv")
  }
  
  #write model to file
  if(plots==TRUE){
    
  writeRaster(mod,"Suitability.tif",overwrite=TRUE)
  write.csv(species.data,"SpeciesData.csv")}
  
  return(species.data)
}

co_occur<-function(siteXsppm){
  
  #reshape data
  colnames(siteXsppm)<-c("Species","Locality","P_A")
  
  f<-split(siteXsppm,siteXsppm$Locality)
  
  #remove the localities with only 1 species
  f<-f[sapply(f,function(x) sum(x$P_A)) > 1]
  
  #Get closest related species in each site
  closest<-lapply(f,function(x){
    #species present
    pres<-x$Species[x$P_A ==1]
    dis<-ctrx[rownames(ctrx) %in% trx$tip.label,colnames(ctrx) %in% pres]
    apply(dis,1,function(y){
      as.matrix(min(y[!y==0]))
    })
  }) 
  
  #melt each list, a bit ugly
  tomerge<-lapply(closest,function(x){
    dm<-melt(as.matrix(x))[,-2]
    colnames(dm)<-c("Species","Phylo.Relatedness")
    return(dm)})
  
  tomerge<-melt(tomerge,id.var=c("Species"))
  
  tomerge<-dcast(tomerge,...~variable)
  colnames(tomerge)<-c("Species","Locality","Phylo.Relatedness")
  
  #need to be characters not factors
  siteXsppm$Locality<-as.character(siteXsppm$Locality)
  tomerge$Locality<-as.character(tomerge$Locality)
  
  PA_phylo<-merge(siteXsppm,tomerge,by=c("Locality","Species"))
  PA_phylo$P_A<-as.numeric(PA_phylo$P_A)
  
  #remove Inf lengths, no co-occurrence
  PA_phylo[!is.finite(PA_phylo$Phylo.Relatedness),"Phylo.Relatedness"]<-NA
  
  return(PA_phylo)
}