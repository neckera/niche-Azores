# load libraries #
library(tidyverse)
library(ggplot2)
library(hypervolume)

# load results from script 2 #
load("NicheSizeResults.RData")

# create a dataframe containing all existing combinations of the climatic varibales
bio<-spec[!duplicated(spec$temp),c("temp","bio1", "bio2", "bio12", "bio15")] # already transformed bio coordinates

# Create a list to save niche.fill results
niche.fill.25 <- list()

# temp is a vector with the name of the species
temp<-names(volumes)
# n is the number of repetitions (100 in this case)
n<-100

# nested loops to calculate niche fill #
for (i in temp){
  for (j in 1:n){
    included <- hypervolume_inclusion_test(hv = rare.res.25[[i]][j][[1]], 
                                           # points are all island coordinates for our 4D (not geographical but climatic)
                                           points = bio[,2:5], 
                                           # fraction of random points sampled from the hypervolume for the stochastic inclusion test
                                           reduction.factor = 1, 
                                           # type of calculation. Accurate is KDE (better but slower)
                                           fast.or.accurate ="fast", 
                                           # use only if selected "fast" before
                                           fast.method.distance.factor = 1,
                                           # use only if selected "accurate" before
                                           # accurate.method.threshold = quantile(volumes[[i]]@ValueAtRandomPoints,0.5), 
                                           # print process if T, silent if F
                                           verbose = FALSE)
    niche.fill.25[[i]][[j]] <- data.frame(name=i,points.in.niche=sum(included),share.included=sum(bio$temp[included] %in% spec$temp[spec$names_repl==i])/sum(included))
  }
  
  # calculate hypervolume of the island
  volumes.island <- hypervolume_svm(bio[,2:5], verbose = FALSE)  
  
  # save results in list
  niche.fill[[i]]<-do.call(rbind.data.frame, niche.fill.25[[i]])
  
  # show the progress: which sp was the last calculated and the percentage
  cat(paste("species = ", i,", ", round(which(temp==i)/length(temp)*100,1), " %"), "\n")
}

# remove non needed objects
rm(bio, climate, coord, islands, M, newproj, niche.fill.25, rare.res.25, spec.coor, volumes, volumes.island, 
   b, directory, i, included, j, n, N, temp, rare, temp2)

# save results #
save.image("NicheFillResults.RData")