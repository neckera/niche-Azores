#############
# load climatic and spec data #
source('directory/1. Extract climate data in occurrences_20112020.R', chdir = TRUE)

#############

#######################
# prepare climate data #
# visually explore the correlation of the bioclim variables
library(corrplot)
library(RColorBrewer)
M <-cor(climate[,1:19])
corrplot(M, type="upper", order="hclust", col=brewer.pal(n=8, name="RdYlBu"))

library(vegan)
# standardize the variables selected
env<-decostand(climate[,c("bio1","bio2","bio12","bio15")],method="standard")
# create a column in the climatic df containing the coordinates as text
env$temp<-paste(climate[,c("Lon_New")], climate[,c("Lat_New")])
# create a column in spec with the same structure as the previous to paste the climate data on it
spec$temp<-paste(spec[,c("Lon_New")], spec[,c("Lat_New")])
# merge the climatic data in the spec df. This way, we remove the points with NA in climate 
spec<-merge(spec,env,by="temp")
rm(env)
#######################

###############################
# List to save results #
volumes<-list()
rare.res.25<-list()
spec.coor<-list()
library(hypervolume)

# create function rare, which calculates the hypervolume of species according to n and size parameters
# n defines the number of times the hypervolume has to be calculated
# size defines the number of occurrences used to calculate the hypervolume
rare<-function(x,n,size){
  sub.volume<-list()
  if(nrow(x)>=size){
    for(i in 1:n){
      sub.env <- x[sample(1:nrow(x), size),]	
      sub.volume[[paste("run",i)]]<-hypervolume_svm(sub.env, verbose = FALSE)
    }} else {sub.volume<-NA}
  sub.volume
}

# extract species with at least 25 occurrences
temp<-names(table(spec$names_repl))[table(spec$names_repl)>24]

N<-rep(NA, length(temp))
names(N)<-temp

for(i in temp){
  # extract the climatic and presence data of an species
    env <- spec[as.character(spec$names_repl)==i,c("bio1","bio2","bio12","bio15")]
  N[names(N)==i]<-nrow(env)
    # extract the median value for the climatic data
  spec.coor[[i]]<-apply(env,2,median)
  # calculate the hypervolume using svm method and save it. Verbose = F stops the output of the function
  volumes[[i]] <- hypervolume_svm(env, verbose = FALSE)
  
  # sample 100 times for each species #
  rare.res.25[[i]]<-rare(env,100,25)
  rm(env)
  
  # show the progress: which sp was the last calculated and the percentage
  cat(paste("species = ", i,", ", round(which(temp==i)/length(temp)*100,1), " %"), "\n")
}
###############################


###############################
# create a dataframe to save the results of the hypervolume size calculation
dati<-data.frame(name=names(unlist(lapply(volumes,function(x){x@Volume}))),
                 volume=NA,dispersal=NA,N=NA)

# unlist hypervolume data. Dispersal is the dispersal type from Schaefer et al. 2011, 
# N is the number of points included, 
# volume is the mean of the 100 repetitions
for(i in dati$name){
  
  temp2<-as.matrix(spec[as.character(spec$names_repl)==i,][1,c("Cat.Schaef")])
  dati[dati$name==i,"dispersal"] <- temp2[1]
  dati[dati$name==i,"N"]<-N[names(N)==i]
  dati[dati$name==i,"volume"]<- mean(unlist(lapply(rare.res.25[[i]],function(x){x@Volume})))
  dati[dati$name==i,"volume.sd"]<- sd(unlist(lapply(rare.res.25[[i]],function(x){x@Volume})))
}
# remove temporal objects
rm(temp2,i,temp)
###############################
# save resulting object as an .RData
save.image("NicheSizeResults.RData")