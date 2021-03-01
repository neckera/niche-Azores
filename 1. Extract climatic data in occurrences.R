rm(list=ls())
library("raster")

# create a list to save all the data
islands <- list()

# create an object to save the directory
directory <- ""
# save new CRS projection for the data
newproj  <-  crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# loop to load all the data
for (i in dir(directory)){
  # load each of the bioclim variables (19 per island)
  for (b in 1:19) {eval(parse(text=paste("bio",b,"  <-  raster('", directory,i,"/Normal_Observado/bio",b,".tif')",sep=""))) 	}
  
  # Create an empty stack
  clim_cur  <-  stack() 
  # ???
  for (b in 1:19) {eval(parse(text=paste("clim_cur  <-  addLayer(clim_cur,bio",b,")",sep="")))}
  # save the names of ecah climatic raster layer
  names(clim_cur)  <-  paste("bio",1:19,sep="")
  # reproject the raster
  clim_cur1  <-  projectRaster(clim_cur, crs=newproj)
  # save the reprojected raster in the islands list
  islands[[i]] <- clim_cur1
  # remove all the objects created
  rm(clim_cur,clim_cur1)
  rm(list=paste("bio",1:19,sep=""))
  # plot in which island is the loop
  cat(paste("island = ", i), "\n")
}

# clean the environment of non-needed things
rm(i,b,directory)

# read the occurrence data
spec <- read.csv("")

# extract the coordinates of the cells where there is AT LEAST one occurrence
coord<-spec[!duplicated(paste(spec[,c("Lon_New")], spec[,c("Lat_New")])),c("Lon_New","Lat_New")]

# create a list to save the climatic data per each occurrence
climate.list <- list()

# create a loop to extract the climatic data per each occurrence
for(i in names(islands)){
  # extract the climatic data per each cell with at least an occurrence
  a<-raster::extract(islands[[i]], coord)	
  # save a as a data.frame instead of matrix
  a <- as.data.frame(a)
  # create a column with the name of the island
  a$island <- i	
  # merge the coordinates of the data with the climatic data
  a <- cbind(a,coord)
  # save only the climatic data of the coordinates (those which are no "NA")
  climate.list[[i]] <- a[!is.na(a$bio1),]
  # remove all the objects created
  rm(a)
}


# create a DF containg the climatic data + name of the island + coordinates of the climatic data
climate  <-  do.call("rbind", climate.list)
# remove the previously created list
rm(climate.list)