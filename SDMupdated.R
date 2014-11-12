#########Biomod2 - Ben Weinstein, Stony Brook University 10/11/2012

#Install packages - only needs to be done the first run
#install.packages("biomod2", repos = "http://R-Forge.R-project.org", dependencies = TRUE)
#install.packages("gam")


#Wrap this into a function to be called from another script

SDM_SP<-function(cell_size,inLocalities,envfolder,savefolder){
  
  #If you have already installed, let's start here.
  #Call the packages we are going to need in this tutorial
  library(biomod2)
  library(maptools)
  library(ggplot2)
  library(reshape)
  library(raster)
  library(rgdal)
  library(stringr)
  library(foreach)
  library(GGally)
  library(plyr)
  
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
  
  #Import environmental data from worldclim, three variables
  myExpl <- c(paste(envfolder,"bio_1.bil",sep="/"),
              paste(envfolder,"bio_12.bil",sep="/"),
              paste(envfolder,"bio_15.bil",sep="/"))
  
  
  myExpl<-stack(myExpl)
  
  #Just get the clean localities
  loc_clean<-PAdat[PAdat$SpatialCheck=="Y" & PAdat$MapDecision %in% levels(PAdat$MapDecision)[!levels(PAdat$MapDecision) %in% "REJECT"],]
  
  #To be safe, take out the bogota localities, likely just bad musueum records
  loc_clean<-loc_clean[!loc_clean$LOCALITY %in% "BOGOTA",]
  
  extPoint<-SpatialPoints(loc_clean[,c("LONGDECDEG","LATDECDEG")])
  #exte<-extent(c(-81.13411,-68.92061,-5.532386,11.94902))
  
  #Crop by this layer, 
  myExpl<-crop(myExpl,extPoint)
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
  
  cl<-snow::makeCluster(10,"SOCK")
  doSNOW::registerDoSNOW(cl)
  system.time(niche_loop<-foreach::foreach(x=1:10,.packages=c("reshape2","biomod2","plyr"),.errorhandling="pass") %dopar% {
    sink(paste("logs/",paste(spec[[x]],".txt",sep=""),sep=""))
    
    #remove sites that have no valid records
    #For the moment, only get the clean records from Decision, or the cleaned map localities. 
    #print(paste("Start Time is",Sys.time()))
    ###############Step 1) Get presence records for species
    ###############Step 1) Get presence records for species
    PA_species<-loc_clean[loc_clean$SPECIES %in% spec[x],]
    
    #How many records do we have?
    records<-nrow(PA_species)
    
    #Have the species names as a vector that will helpful later
    spname<-levels(factor((PA_species$SPECIES)))
    
    #get unique localities 
    pts<-aggregate(PA_species,list(PA_species$LOCALITY,PA_species$LATDECDEG,PA_species$LONGDECDEG),FUN=mean)
    pts.sp<-SpatialPointsDataFrame(pts[,c("LONGDECDEG","LATDECDEG")],pts)
    
    #########################The above gets you presence only records
    #To get presence/psuedoabscence  
    #Get the presence absence matrix
    loc_matrix<-table(loc_clean$LOCALITY,loc_clean$SPECIES)
    
    #Select the species you'd like
    sp_matrix<-as.data.frame(melt(loc_matrix[,spec[x]]))
    unique.loc<-unique(loc_clean[,c("LOCALITY","LONGDECDEG","LATDECDEG")])
    p_a<-merge(sp_matrix,unique.loc,by.x="row.names",by.y="LOCALITY")
    p_a<-p_a[!duplicated(p_a),]
    
    #name the columns
    colnames(p_a)<-c("Locality","Response","LONGDECDEG","LATDECDEG")
    p_a[p_a$Response > 1,"Response"]<-1
    p_a<-p_a[!p_a$Locality=="",]
    p_a<-aggregate(p_a,list(p_a$Locality),FUN=mean)
    
    #the 0 are not true absences, they are NA's psuedoabsences
    p_a[p_a$Response == 0,"Response"]<-NA
    
    #we want 2000 psuedoabsences, randomly pick 2000 NA rows. 
    if(length(which(is.na(p_a$Response))) < 2000){
      psuedos<-which(is.na(p_a$Response))
    }
    
    if(length(which(is.na(p_a$Response))) > 2000){
      psuedos<-sample(which(is.na(p_a$Response)),2000)
    }
              
              pres_pts<-which(!is.na(p_a$Response))
              
              #Format the data
              p_a<-p_a[c(pres_pts,psuedos),]
              
              
              myBiomodData <- BIOMOD_FormatingData(resp.var =  p_a[,"Response"],
                                                   expl.var = stack(myExpl.crop),
                                                   resp.xy =  p_a[,c("LONGDECDEG","LATDECDEG")],
                                                   resp.name = gsub(" ","_",spec[x]))
              
              print(myBiomodData)
              
              
              #Define modeling options
              myBiomodOption <- BIOMOD_ModelingOptions(    
                MAXENT = list( path_to_maxent.jar = savefolder,
                               maximumiterations = 200,
                               visible = TRUE,
                               linear = TRUE,
                               quadratic = TRUE,
                               product = TRUE,
                               threshold = TRUE,
                               hinge = TRUE,
                               lq2lqptthreshold = 80,
                               l2lqthreshold = 10,
                               hingethreshold = 15,
                               beta_threshold = -1,
                               beta_categorical = -1,
                               beta_lqp = -1,
                               beta_hinge = -1,
                               defaultprevalence = 0.5)
              )
              
              
              #Give current project a name, so we can go get the files later
              projnam<-'current'
              myBiomodModelOut<-BIOMOD_Modeling( myBiomodData, 
                                                 models = c("GBM","GLM","MAXENT"), 
                                                 models.options = myBiomodOption, 
                                                 NbRunEval=1, 
                                                 DataSplit=80, 
                                                 Yweights=NULL, 
                                                 VarImport=3, 
                                                 #models.eval.meth = c('ROC'),
                                                 models.eval.meth = c('ROC',"TSS","KAPPA"),
                                                 SaveObj = TRUE )
              
              # get all models evaluation                                     
              myBiomodModelEval <- get_evaluations(myBiomodModelOut)
              
              # print the dimnames of this object
              dimnames(myBiomodModelEval)
              
              # let's print the ROC and TSS scores of all selected models, get the mean value for all the combined runs.
              stat<-myBiomodModelEval["ROC","Testing.data",,"Full",]
              full_stat<-myBiomodModelEval[,"Testing.data",,"Full",]
              
              #need to write this to file
              filename<-paste(paste(getwd(),gsub(" ",".",spec[x]),sep="/"),"ModelRocEval.csv",sep="/")
              write.csv(cbind(spec[x],stat),filename)
              
              filename<-paste(paste(getwd(),gsub(" ",".",spec[x]),sep="/"),"ModelEval.csv",sep="/")
              write.csv(cbind(spec[x],full_stat),filename)
              
              #Let's look at variable importance
              m.var<-melt(get_variables_importance(myBiomodModelOut)[,,"Full",])
              c.var<-cast(m.var,X1~X2)
              
              #Write variable importance to file
              filename<-paste(paste(getwd(),gsub(" ",".",spec[x]),sep="/"),"VarImportance.csv",sep="/")
              write.csv(cbind(c.var,spec[x]),filename)
              
              # projection over the globe under current conditions  
              
              #Ensemble model
              myBiomodEM <- BIOMOD_EnsembleModeling( 
                modeling.output = myBiomodModelOut,
                chosen.models = 'all',
                em.by='all',
                eval.metric = c('ROC'),
                eval.metric.quality.threshold = c(0.75),
                prob.mean = T,
                prob.cv = T,
                prob.ci = T,
                prob.ci.alpha = 0.05,
                prob.median = T,
                committee.averaging = T,
                prob.mean.weight = T,
                prob.mean.weight.decay = 'proportional' )
              
              #################
              #Project SDM into env projections
              #All the gcm and current worldclim layer (first) are put together in a list
              #Loop through this list, only run if it has not been run before
              #name it corrently
              
              bio_project<-function(GCM,nam){
                paste("Running Env", nam)
                #if(!gsub(" ",".",spec[x]) %in% completed_models[[nam]]){
                myBiomodProjection<- BIOMOD_Projection(
                  modeling.output = myBiomodModelOut,
                  new.env = projEnv[[1]],
                  proj.name = "current",
                  selected.models = 'all',
                  binary.meth = 'ROC',
                  compress = 'xz',
                  clamping.mask = T)
                
                #Ensemble projection
                EnsBas<-BIOMOD_EnsembleForecasting(projection.output = myBiomodProjection, EM.output = myBiomodEM)
                
              }
              
              mapply(bio_project,projEnv,names(projEnv))
              ##################################
              
              #end file output
              #print(paste("End Time is",system.time()))
              sink()
              return(stat)
  })
stopCluster(cl)

print("ModelsComplete")

############################Get model evaluation stats
###########################

#Get the ROC model evaluation from file
model_eval<-list.files(full.name=TRUE,recursive=T,pattern="ModelRocEval.csv")
model_eval<-rbind.fill(lapply(model_eval,read.csv))
colnames(model_eval)<-c("Model","Species","Stat")
model_eval<-melt(model_eval,id.var=c("Model","Species","Stat"))

#Plot
ggplot(model_eval, aes(x=Species,y=Model,fill=Stat)) + geom_tile() + scale_fill_gradient("ROC",limits=c(0,1),low="blue",high="red",na.value="white") + theme(axis.text.x=element_text(angle=-90))
ggsave("ModelROCEvaluations.jpeg",height=5,width=20,dpi=300)

model_thresh<-sapply(seq(.5,.95,.05),function(x){
  table(model_eval$Stat > x,model_eval$Model)["TRUE",]
})

colnames(model_thresh)<-seq(.5,.95,.05)
model_thresh<-melt(model_thresh)

names(model_thresh)<-c("Model","ROC_Threshold","Number_of_Species")
ggplot(model_thresh,aes(x=ROC_Threshold,y=Number_of_Species,col=Model)) + geom_line(size=1.5) + geom_point() + geom_text(aes(label=Number_of_Species),vjust=4,size=5) + theme_bw()
ggsave("ModelThresholding.pdf",dpi=300,height=8,width=8)

#Get the all model evaluations from file
model_eval<-list.files(full.name=TRUE,recursive=T,pattern="ModelEval.csv")
model_eval<-rbind.fill(lapply(model_eval,read.csv))
colnames(model_eval)<-c("Metric","Species","GBM","GLM","MAXENT")
model_eval<-melt(model_eval,id.var=c("Metric","Species"))

#Plot with facets
p<-ggplot(model_eval, aes(x=Species,y=variable,fill=value)) + geom_tile() + facet_wrap(~Metric) + theme(axis.text.x=element_text(angle=-90))
p + scale_fill_gradient("Score",limits=c(0,1),low="blue",high="red",na.value="white")+ theme_bw() + theme(axis.text.x = element_blank()) 
ggsave("ModelEvaluations.jpeg",height=6,width=11,dpi=300)

#How correlated are model scores
cast_mods<-cast(model_eval, Species~...)

#GBM correlations
ggpairs(cast_mods[,c(2,5,8)],lower=list(continuous="smooth",method="lm"))

svg("GBMmetriccor.svg",height=10,width=10)
ggpairs(cast_mods[,c(2,5,8)],lower=list(continuous="smooth",method="lm"))
dev.off()

#GLM correcation
ggpairs(na.omit(cast_mods[,c(3,6,9)]),na.rm=TRUE,lower=list(continuous="smooth",method="lm"))

svg("GLMmetriccor.svg",height=10,width=10)
ggpairs(na.omit(cast_mods[,c(3,6,9)]),na.rm=TRUE,lower=list(continuous="smooth",method="lm"))
dev.off()

#Maxent Correlation
ggpairs(cast_mods[,c(4,7,10)],lower=list(continuous="smooth",method="lm"))

svg("Maxentmetriccor.svg",height=10,width=10)
ggpairs(cast_mods[,c(4,7,10)],lower=list(continuous="smooth",method="lm"))
dev.off()

ggplot(model_eval,aes(x=variable,y=value)) + facet_wrap(~Metric) + geom_boxplot()+ geom_smooth(method="lm") + theme_bw() + labs(y="Score",x="Model")


#Get the variable importance from file
varI<-list.files(full.name=TRUE,recursive=T,pattern="VarImportance.csv")
varI<-rbind.fill(lapply(varI,read.csv))
varI<-varI[,-1]

#Melt variable for plotting
mvar<-melt(varI)
colnames(mvar)<-c("Bioclim","Species","Model","value")

#Plot variable importance across all models
ggplot(mvar, aes(x=Species,y=Bioclim,fill=value)) + geom_tile() + scale_fill_gradient(limits=c(0,1),low="blue",high="red",na.value="white") + theme(axis.text.x=element_text(angle=-90)) + facet_grid(Model ~ .) + ggtitle("Variable Importance")

ggsave("VariableImportance.jpeg",height=5,width=15,dpi=300)
}

print("SDM Function Defined")