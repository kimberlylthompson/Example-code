#########################################################################
###########            Prediction summaries              ################
###########      Mean, Median and CV for each model's    ################
###########                Raster Stack                  ################
###########     (not yet with greenhouse correction)     ################
###########   for models that did not include latitude   ################
###########   as a predictor, but do have mean density   ################
#########################################################################

# Clear workspace
rm(list = ls() ) 
gc() #releases memory

library(raster)
library(sp)
library(lubridate)

#########################################################################
# Read in the raster stacks for each model

setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Uncorrected Stacks")

# Load BRT prediction rasters
ext_density <- stack("Uncorrected External Model.tif")

unc_H0 <- stack("Uncorrected H0 Model.tif")

unc_H3 <- stack ("Uncorrected H3 Model.tif") 

unc_H5 <- stack ("Uncorrected H5 Model.tif") 


########################################################################
# Load in the states shapefile and make sure everything matches up
#Import shapefile of these states
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/GIS")
states<-rgdal::readOGR("Great_Lakes_States.shp" )

# Convert projection so that it does match:
states <- spTransform( states, proj4string( unc_H0 ) )


# Make functions to add states to each layer in the stack

# Just states
fun <- function() {
  plot(states, add=TRUE)
} 


#############################################################################
########### Taking the mean, median, and CV for each model  #################
#############################################################################

# List the uncorrected raster stacks
path <- "L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Uncorrected Stacks/"
predictions <- list.files(path, pattern = "Uncorrected")

# Write the new functions for mean, median, and cv to make sure that NA values are ignored
mean.na <- function (x) {
  mean(x, na.rm=TRUE)
}

median.na <- function (x) {
  median(x, na.rm=TRUE)
}

cv.na <- function (x) {
  cv(x, na.rm=TRUE)
}

start.time <- Sys.time() # Takes about 3.4 hours

for (i in 1:length(predictions)) {
  #Set the working directory
  setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Uncorrected Stacks")
  
  #Read in the stack
  modelstack <- stack(predictions[[i]])
  
  #Take the mean, median, and CV of the stack
  mean_stack <- overlay(modelstack, fun=mean.na)
  median_stack <- overlay(modelstack, fun=median.na)
  cv_stack <- overlay(modelstack, fun=cv.na)
  
  #Write the rasters to a folder
  setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Uncorrected Summaries")
  
  writeRaster(mean_stack, paste("Mean_", unlist(predictions[[i]]), sep=""), overwrite=TRUE)
  writeRaster(median_stack, paste("Median_", unlist(predictions[[i]]), sep=""), overwrite=TRUE)
  writeRaster(cv_stack, paste("CV_", unlist(predictions[[i]]), sep=""), overwrite=TRUE)
  
  print(i)
}

end.time <- Sys.time()
end.time - start.time
