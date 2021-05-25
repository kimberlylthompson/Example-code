#########################################################################
###########            Corrections to mean, median,      ################
###########                  and CV summaries            ################
###########         for each model to account for        ################
###########           artefact of the greenhouse         ################
#########################################################################

# Clear workspace
rm(list = ls() ) 
gc() #releases memory

library(raster)
library(sp)
library(lubridate)


#########################################################################
#############            Greenhouse Correction           ################
#########################################################################

# To determine how to do the greenhouse correction compare the plots of 
# external without density and H0 to see whether the greenhouse causes
# overprediction or underpredicion - result is that it actually varies by cover
# type.

# Goal of correction is to get the house 0 prediction to be equivalent to that 
# of the external (because this takes away the effect of the greenhouse structure)

# formula: Correction = House 0 - external
# therefore: House 0 - Correction = external
# extending to other houses: House 3 (or House 5) - Correction = predictions that remove
# effect of the greenhouse

#############################
#######   Mean     ##########
#############################

# Read in the raster layers
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Uncorrected Summaries")

# Load BRT prediction rasters - start with mean
ext_mean <- raster::raster("Mean_Uncorrected External Model.tif")
H0_mean <- raster::raster("Mean_Uncorrected H0 Model.tif")
H3_mean <- raster::raster("Mean_Uncorrected H3 Model.tif")
H5_mean <- raster::raster("Mean_Uncorrected H5 Model.tif")


# Obtaining the greenhouse correction: H0-external
ghcorrection <- H0_mean - ext_mean


# Apply the correction formula to each of the greenhouse mean rasters
H0_corrected <- H0_mean - ghcorrection
H3_corrected <- H3_mean - ghcorrection
H5_corrected <- H5_mean - ghcorrection

# Write the corrected rasters
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Corrected Summaries")

writeRaster(H0_corrected, "H0_mean_corrected.tif", overwrite = TRUE)
### NoTE - H0_corrected is now equivalent to external (current subnivium)

writeRaster(H3_corrected, "H3_mean_corrected.tif", overwrite = TRUE)
writeRaster(H5_corrected, "H5_mean_corrected.tif", overwrite = TRUE)

# What this does is make H0_corrected equivalent to ext
# To check this, I can take the difference between them ad it should be 
# a raster with only 0 values
# check <- H0_corrected - ext_mean
# plot(check, addfun=fun, main="Correction Check")
# 
# # it works!

#############################
#######   Median   ##########
#############################

# Read in the raster layers
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Uncorrected Summaries")

# Load BRT prediction rasters
ext_median <- raster::raster("Median_Uncorrected External Model.tif")
H0_median <- raster::raster("Median_Uncorrected H0 Model.tif")
H3_median <- raster::raster("Median_Uncorrected H3 Model.tif")
H5_median <- raster::raster("Median_Uncorrected H5 Model.tif")


# Obtaining the greenhouse correction: H0-external no density
ghcorrection.med <- H0_median - ext_median

# Now this means that if I'm trying to achieve equivalence between H0 and external
# and starting with H0, that the correction formula would be
# House 0 (or H3 or H5) - ghcorrection = external

# Apply the correction formula to each of the greenhouse mean rasters
H0med_corrected <- H0_median - ghcorrection.med
H3med_corrected <- H3_median - ghcorrection.med
H5med_corrected <- H5_median - ghcorrection.med

# Write the corrected rasters
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Corrected Summaries")

writeRaster(H0med_corrected, "H0_median_corrected.tif", overwrite = TRUE)
writeRaster(H3med_corrected, "H3_median_corrected.tif", overwrite = TRUE)
writeRaster(H5med_corrected, "H5_median_corrected.tif", overwrite = TRUE)


#############################
#######   CV       ##########
#############################

# Read in the raster layers
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Uncorrected Summaries")

# Load BRT prediction rasters
ext_cv <- raster::raster("CV_Uncorrected External Model.tif")
H0_cv <- raster::raster("CV_Uncorrected H0 Model.tif")
H3_cv <- raster::raster("CV_Uncorrected H3 Model.tif")
H5_cv <- raster::raster("CV_Uncorrected H5 Model.tif")


# Obtaining the greenhouse correction: H0-external no density
ghcorrection.cv <- H0_cv - ext_cv

# Now this means that if I'm trying to achieve equivalence between H0 and external
# and starting with H0, that the correction formula would be
# House 0 (or H3 or H5) - ghcorrection = external

# Apply the correction formula to each of the greenhouse mean rasters
H0cv_corrected <- H0_cv - ghcorrection.cv
H3cv_corrected <- H3_cv - ghcorrection.cv
H5cv_corrected <- H5_cv - ghcorrection.cv

# Write the corrected rasters
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Corrected Summaries")

writeRaster(H0cv_corrected, "H0_cv_corrected.tif", overwrite = TRUE)
writeRaster(H3cv_corrected, "H3_cv_corrected.tif", overwrite = TRUE)
writeRaster(H5cv_corrected, "H5_cv_corrected.tif", overwrite = TRUE)



########################################################################
# This takes care of correcting the summaries, but I still need to make
# corrections to the individual rasters for each day for each model.


######### Corrections for each individual day for each model
######### Each day's correction formula will be 
######### ghcorrection <- H0 (for that day) - ext (for that day)
######### H0 - ghcorrection = corrected day raster

# Clear workspace
rm(list = ls() ) 
gc() #releases memory

# Make a list of the daily raster stacks
path <- "L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/"
external.list <- list.files(paste(path, "External", sep=""), pattern = "X")
H0.list <- list.files(paste(path, "House 0", sep=""), pattern = "X")
H3.list <- list.files(paste(path, "House 3", sep=""), pattern = "X")
H5.list <- list.files(paste(path, "House 5", sep=""), pattern = "X")


######################################################################################
####                               HOUSE 0                                        ####
######################################################################################

start_time <- Sys.time() 

for (i in 1:length(external.list)) {
  
  # Read in external day's layer
  setwd(paste(path, "External", sep=""))
  external <- raster::raster(external.list[[i]])
  
  # Read in H0 day's layer
  setwd(paste(path, "House 0", sep=""))
  H0 <- raster::raster(H0.list[[i]])
  
  # Find the greenhouse correction
  correction <- H0 - external
  
  # Adjust the H0 raster by the correction
  corrected.external <- H0 - correction
  
  # Write the resulting corrected raster
  setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Corrected Dailies/House 0")
  
  # To get the file name, split the string of the predictor list
  partialname <- substr(unlist(H0.list[[i]]), 1, 11)
  
  writeRaster(corrected.external, paste(partialname, "_H0cor.tif", sep = ""), overwrite=TRUE)
  
  print(i)
  
}


rm(correction)
rm(external)
rm(H0)

######################################################################################
####                               HOUSE 3                                        ####
######################################################################################

for (i in 1:length(external.list)) {
  
  # Read in external day's layer
  setwd(paste(path, "External", sep=""))
  external <- raster::raster(external.list[[i]])
  
  # Read in H0 day's layer
  setwd(paste(path, "House 0", sep=""))
  H0 <- raster::raster(H0.list[[i]])
  
  # Read in H3 day's layer
  setwd(paste(path, "House 3", sep=""))
  H3 <- raster::raster(H3.list[[i]])
  
  # Find the greenhouse correction
  correction <- H0 - external
  
  # Adjust the H3 raster by the correction
  corrected.H3 <- H3 - correction
  
  # Write the resulting corrected raster
  setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Corrected Dailies/House 3")
  
  # To get the file name, split the string of the predictor list
  partialname <- substr(unlist(H0.list[[i]]), 1, 11)
  
  writeRaster(corrected.H3, paste(partialname, "_H3cor.tif", sep = ""), overwrite=TRUE)
  
  print(i)
  
}



######################################################################################
####                               HOUSE 5                                        ####
######################################################################################



for (i in 1:length(external.list)) {
  
  # Read in external day's layer
  setwd(paste(path, "External", sep=""))
  external <- raster::raster(external.list[[i]])
  
  # Read in H0 day's layer
  setwd(paste(path, "House 0", sep=""))
  H0 <- raster::raster(H0.list[[i]])
  
  # Read in H5 day's layer
  setwd(paste(path, "House 5", sep=""))
  H5 <- raster::raster(H5.list[[i]])
  
  # Find the greenhouse correction
  correction <- H0 - external
  
  # Adjust the H5 raster by the correction
  corrected.H5 <- H5 - correction
  
  # Write the resulting corrected raster
  setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Corrected Dailies/House 5")
  
  # To get the file name, split the string of the predictor list
  partialname <- substr(unlist(H0.list[[i]]), 1, 11)
  
  writeRaster(corrected.H5, paste(partialname, "_H5cor.tif", sep = ""), overwrite=TRUE)
  
  print(i)
  
}

end_time <- Sys.time()
end_time - start_time
