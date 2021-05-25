#########################################################################
###########   Comparing the distributions of             ################
###########              subnivium extent                ################
###########           in each warming scenario           ################
#########################################################################

# Clear workspace
rm(list = ls() ) 
gc() #releases memory

library(raster)
library(sp)
library(lubridate)
library(tidyr)
library(car)
library(ggplot2)
library(FSA)
library(rcompanion)
library(PMCMRplus)
library(gridExtra)
library(grid)


# Read in the median extent summary for each treatment
# (current [corrected H0], House 3, and House 5)

setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Spatial Prediction/No Lat With Density/Corrected Summaries")

median0 <- raster :: raster("H0_median_corrected.tif")

median3 <- raster :: raster("H3_median_corrected.tif")

median5 <- raster :: raster("H5_median_corrected.tif")


###########################################################################
#      make adjustments to each raster to get them ready for calculations #
###########################################################################

# Make a stack of the mean and median rasters
rasterlist <- stack(median0, median3, median5)

# Because of the greenhouse correction, some cells in 3 and 5 will have values less than 0 and greater than 1
# Reclassify these so that all values are between 0 and 1
# Only do this for the mean and median rasters!
rasterlist[rasterlist > 1] <- 1
rasterlist[rasterlist < 0] <- 0

# Remove the single raster layers for a cleaner workspace
rm(median0, median3, median5)

# For plotting transform each raster into degrees so that the axes will look
# better and be clearer to understand
wgs84 = sp:: CRS("+init=epsg:4326")
rasterlist <- projectRaster(from=rasterlist, crs =  wgs84)

# Load in the states shapefile and make sure everything matches up
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/GIS")
states<-rgdal::readOGR("Great_Lakes_States.shp" )

# Convert projection so that it does match:
states <- spTransform( states, proj4string( rasterlist ) )

# Read in the latitude raster so that I can use it to crop the predictive map
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Latitude")
latitude.deg <- raster :: raster("Latitude_degrees.tif") 

# Make sure it matches the spatial projection of the degrees predictive rasters
latitude.deg <- projectRaster(from=latitude.deg, crs =  wgs84)

# Load in great lakes to add to plot
# Load in the states shapefile and make sure everything matches up
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/GIS")
lakes<-rgdal::readOGR("Great_Lakes.shp" )

# Convert projection so that it does match:
lakes <- spTransform( lakes, proj4string( rasterlist ) )

# crop the states
newstates <- raster::crop(states, extent(-97.2291, -75.65702, 41.5, 49.37173))


# Crop the predictive raster to the states extent and then mask it so only states show
rasterlist <- crop(rasterlist, extent(newstates))
rasterlist <- raster :: mask(rasterlist, newstates)


#####################################################
#######      Testing the distributions         ######
#####################################################

### None of the data from the 3 scenarios are normally distributed. Therefore, it is not possible to use either a t-test
### or a 1-way anova. The nonparametric alternative to the one-way anova is the Kruskal-Wallis rank rum test.
### https://rcompanion.org/rcompanion/d_06.html - Information of KW tests in R and post hoc testing
### However, upon further reading at https://statistics.laerd.com/spss-tutorials/kruskal-wallis-h-test-using-spss-statistics.php
### one assumption of the Kruskal wallis test is that there is independence of observations, so no relationship between the
### observations in each group with no cell being in more than one scenario - This is violated because of spatial autocorrelation
### as well as multiple measurements of the same cells.
### Therefore, I need to switch to the Friedman test, which is the non-parametric alternative to the one-way ANOVA
### with repeated measures. - It's also a rank test so the null should remain the same.
### Friedman testin R: https://rcompanion.org/handbook/F_10.html

### The Friedman test determines if there are differences among groups for two-way data structured in a specific way,
### namely in an unreplicated complete block design.  In this design, one variable serves as the treatment or group variable,
### and another variable serves as the blocking variable.  It is the differences among treatments or groups that we are 
### interested in.  We aren't necessarily interested in differences among blocks, but we want our statistics to take into
### account differences in the blocks.  In the unreplicated complete block design, each block has one and only one
### observation of each treatment. 

###### In this explanation - the blocks are equivalent to the cells. Differences (or similarities) between them need to be 
###### taken into account due to spatial autocorrelation.


## If the the distributions of values for each group have similar shape and spread:
##     Null hypothesis:  The medians of values for each group are equal.
##     Alternative hypothesis (two-sided): The medians of values for each group are not equal.

# Convert each median extent raster to a dataframe
H0 <- as.data.frame(rasterlist[[1]], xy=TRUE)
H3 <- as.data.frame(rasterlist[[2]], xy=TRUE)
H5 <- as.data.frame(rasterlist[[3]], xy=TRUE)

# Add an ID column for each which will represent the unique lat long combo
H0$ID <- seq(from = 1, to=length(H0$x), by=1)
H3$ID <- seq(from = 1, to=length(H3$x), by=1)
H5$ID <- seq(from = 1, to=length(H5$x), by=1)

# Merge dataframes
combined <- merge(H0, H3, by.x = c("x", "y", "ID"), by.y = c("x", "y", "ID"), all=TRUE )
combined <- merge(combined, H5, by.x = c("x", "y", "ID"), by.y = c("x", "y", "ID"), all=TRUE )

# Make sure that ID column is a factor
combined$ID <- factor(combined$ID)

# Reshape data into long form 
combined2 <- gather(combined, Treatment, Probability, H0_median_corrected:H5_median_corrected, factor_key=TRUE)

# Make sure that ID column is a factor
combined2$ID <- factor(combined2$ID)

# Removed the rows in long form data frame where there are NA values
combined2 <- combined2[complete.cases(combined2),]

# Remove the rows in wide form data frame where there are NA values
combined <- combined[complete.cases(combined),]


########### Friedman Test ###################

# Medians and descriptive statistics
Summarize(Probability ~ Treatment, data = combined2)

# Convert combined2 to a matrix
combined2.matrix <- as.matrix(combined2)

friedman.test(Probability ~ Treatment | ID, data = combined2.matrix)

# post hoc Conover test
# Use wide form data frame for the conover test

# Take out the xy columns
combined <- combined[,-c(1:2)]

# Convert to a matrix
conover.matrix <- data.matrix(combined)

# Take out ID column of the matrix
conover.matrix <- conover.matrix[,-1]

# Conover test (same result when using default holm padjust method)
PT <- posthoc.friedman.conover.test(y = conover.matrix,
                                    p.adjust.method = "fdr")


#####################################################
#######      creating histogram figures        ######
#####################################################

extent <- ggplot() +
  geom_histogram(data=combined, aes(x=H5_median_corrected, fill = 'H5'), color = "grey3", binwidth = 0.05) +
  geom_histogram(data=combined, aes(x=H3_median_corrected, fill = 'H3'), color = "grey3", binwidth = 0.05) +
  geom_histogram(data=combined, aes(x=H0_median_corrected, fill = 'H0'), color = "grey3", binwidth = 0.05) +
  geom_vline(xintercept=0.2316998, linetype="dashed", color="indianred3", size = 2) + #H5 median value
  geom_vline(xintercept=0.766747, linetype="dashed", color="maroon1", size = 2) + #H3 median value
  geom_vline(xintercept=0.8138023, linetype="dashed", color="dodgerblue3", size = 2) + #H0 median value
  theme_bw() +
  theme(axis.text.x = element_text(size=22, face="bold")) +
  theme(axis.text.y = element_text(size=22, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) +
  scale_x_continuous(name="\nPredicted Probability") +
  scale_y_continuous(name="Count\n", limits=c(0, 278000),
                     breaks=c(0, 50000, 100000, 150000, 200000, 250000),
                     labels=c("0", "50,000", "100,000", "150,000", "200,000", "250,000")) +
  scale_fill_manual(name='', values = c('H5' = "indianred3",
                                        'H3' = "maroon1",
                                        'H0' = "dodgerblue3"),
                    labels = c("Current conditions", "+3°C", "+5°C")) +
  # theme(aspect.ratio=1) + #aspect ratio expressed as y/x
  # theme(legend.position = c(0.85, 0.60)) +
  # theme(legend.position = "none") +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = (list(size=3))))


setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/Change in Extent")
ggsave("Extent Histogram.jpg", plot=extent, device = "jpeg",
       width=9, height=5, dpi=600)


#### Density plot which looks better than the histogram

# Determining which color scheme looks best - I like this one the best
extent2 <- ggplot() +
  geom_density(data=combined2, aes(x=Probability, fill = Treatment), alpha = 0.3, position = "identity", color = "grey3") +
  geom_vline(xintercept=0.2316998, linetype="dashed", color="#FC4E07", size = 2, alpha = 0.8) + #H5 median value
  geom_vline(xintercept=0.766747, linetype="dashed", color="#E7B800", size = 2, alpha = 0.8) + #H3 median value
  geom_vline(xintercept=0.8138023, linetype="dashed", color="#00AFBB", size = 2, alpha = 0.8) + #H0 median value
  # theme_bw() +
  theme(axis.text.x = element_text(size=16, color="black")) +
  theme(axis.text.y = element_text(size=16, color="black")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="black")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="black")) +
  scale_x_continuous(name="\nPredicted Probability",
                     breaks = c(0, 0.25, 0.50, 0.75, 1),
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_continuous(name="Density\n", limits=c(0, 7),
                     breaks=c(0, 2, 4, 6),
                     labels=c("0", "2", "4", "6")) +
  scale_fill_manual(name='', values = c('H5_median_corrected' = "#FC4E07",
                                        'H3_median_corrected' = "#E7B800",
                                        'H0_median_corrected' = "#00AFBB"),
                    labels = c(" Current conditions", " +3°C", " +5°C")) +
  # theme(aspect.ratio=1) + #aspect ratio expressed as y/x
  # theme(legend.position = c(0.85, 0.60)) +
  # theme(legend.position = "none") +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = (list(size=3))))


setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/Change in Extent")
ggsave("Extent Density plot.jpg", plot=extent2, device = "jpeg",
       width=9, height=3.75, dpi=600)



##########################################################################
############### For graphing with the same height and width  #############
###############           as the extent figure               #############
##########################################################################

### Need to clear the memory and run the code in 
### L:\LabMemberFolders\KimberlyThompson\Ch 2 Climate\Analysis\R Code\Spatial Predictions\Comparing distributions of duration scenarios

### Then clear everything except p0 object
rm(combined, combined2, combined2.matrix, conover.matrix, duration2, H0, H3, H5, lakes, latitude.deg, newstates, PT,
   rasterlist, states, wgs84)

### Then re-run lines 25 - 141 of this script and 214 - 239

# p.ex = extent2

p.ex <- ggplotGrob(extent2)

#Make extent width and height match p0 
p.ex$widths<-p0$widths
p.ex$heights<-p0$heights

#Save the extent plot
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/Change in Extent")
ggsave("Extent Density plot_correct hw.jpg", plot=grid.draw(p.ex), device = "jpeg",
       width=9, height=3.75, dpi=600)



#######################################################################################
############ Creating Individual Distribution Plots to use as insets in the maps ######
#######################################################################################

extent.0 <- ggplot() +
  geom_density(data=combined, aes(x=H0_median_corrected, fill = 'H0'), alpha = 0.7, position = "identity", color = "grey3") +
  # geom_vline(xintercept=0.2316998, linetype="dashed", color="#FC4E07", size = 2, alpha = 0.8) + #H5 median value
  # geom_vline(xintercept=0.766747, linetype="dashed", color="#E7B800", size = 2, alpha = 0.8) + #H3 median value
  geom_vline(xintercept=0.8138023, linetype="dashed", color="purple4", size = 2, alpha = 0.8) + #H0 median value
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_fill_manual(name='', values = c('H0' = "purple4"),
                    labels = '') +
  theme(legend.position = "none")

setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/Change in Extent")
ggsave("Extent Density plot H0 ONLY.jpg", plot=extent.0, device = "jpeg",
       width=5, height=3.75, dpi=600)

extent.3 <- ggplot() +
  geom_density(data=combined, aes(x=H3_median_corrected, fill = 'H3'), alpha = 0.7, position = "identity", color = "grey3") +
  # geom_vline(xintercept=0.2316998, linetype="dashed", color="#FC4E07", size = 2, alpha = 0.8) + #H5 median value
  geom_vline(xintercept=0.766747, linetype="dashed", color="purple4", size = 2, alpha = 0.8) + #H3 median value
  # geom_vline(xintercept=0.8138023, linetype="dashed", color="purple4", size = 2, alpha = 0.8) + #H0 median value
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_fill_manual(name='', values = c('H3' = "purple4"),
                    labels = '') +
  theme(legend.position = "none")

setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/Change in Extent")
ggsave("Extent Density plot H3 ONLY.jpg", plot=extent.3, device = "jpeg",
       width=5, height=3.75, dpi=600)

extent.5 <- ggplot() +
  geom_density(data=combined, aes(x=H5_median_corrected, fill = 'H5'), alpha = 0.7, position = "identity", color = "grey3") +
  geom_vline(xintercept=0.2316998, linetype="dashed", color="purple4", size = 2, alpha = 0.8) + #H5 median value
  # geom_vline(xintercept=0.766747, linetype="dashed", color="purple4", size = 2, alpha = 0.8) + #H3 median value
  # geom_vline(xintercept=0.8138023, linetype="dashed", color="purple4", size = 2, alpha = 0.8) + #H0 median value
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_fill_manual(name='', values = c('H5' = "purple4"),
                    labels = '') +
  theme(legend.position = "none")

setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/Change in Extent")
ggsave("Extent Density plot H5 ONLY.jpg", plot=extent.5, device = "jpeg",
       width=5, height=3.75, dpi=600)



###############################################################################################
###################### Making individual histograms for supplemental figures ##################
###############################################################################################

extent.hist0 <- ggplot() +
  # geom_histogram(data=combined, aes(x=H5_median_corrected, fill = 'H5'), color = "grey3", binwidth = 0.05) +
  # geom_histogram(data=combined, aes(x=H3_median_corrected, fill = 'H3'), color = "grey3", binwidth = 0.05) +
  geom_histogram(data=combined, aes(x=H0_median_corrected, fill = 'H0'), color = "grey3", binwidth = 0.05) +
  # geom_vline(xintercept=0.2316998, linetype="dashed", color="indianred3", size = 2) + #H5 median value
  # geom_vline(xintercept=0.766747, linetype="dashed", color="maroon1", size = 2) + #H3 median value
  geom_vline(xintercept=0.8138023, linetype="dashed", color="dodgerblue3", size = 2) + #H0 median value
  theme_bw() +
  theme(axis.text.x = element_text(size=22, face="bold")) +
  theme(axis.text.y = element_text(size=22, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) +
  scale_x_continuous(name="\nPredicted Probability", limits=c(-0.05, 1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_continuous(name="Count\n", limits=c(0, 278000),
                     breaks=c(0, 50000, 100000, 150000, 200000, 250000),
                     labels=c("0", "50,000", "100,000", "150,000", "200,000", "250,000")) +
  scale_fill_manual(name='', values = c('H0' = "dodgerblue3"),
                    labels = c('')) +
  # theme(aspect.ratio=1) + #aspect ratio expressed as y/x
  # theme(legend.position = c(0.85, 0.60)) +
  theme(legend.position = "none")


setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/Change in Extent")
ggsave("Extent Histogram H0.jpg", plot=extent.hist0, device = "jpeg",
       width=7, height=5, dpi=600)


extent.hist3 <- ggplot() +
  # geom_histogram(data=combined, aes(x=H5_median_corrected, fill = 'H5'), color = "grey3", binwidth = 0.05) +
  geom_histogram(data=combined, aes(x=H3_median_corrected, fill = 'H3'), color = "grey3", binwidth = 0.05) +
  # geom_histogram(data=combined, aes(x=H0_median_corrected, fill = 'H0'), color = "grey3", binwidth = 0.05) +
  # geom_vline(xintercept=0.2316998, linetype="dashed", color="indianred3", size = 2) + #H5 median value
  geom_vline(xintercept=0.766747, linetype="dashed", color="dodgerblue3", size = 2) + #H3 median value
  # geom_vline(xintercept=0.8138023, linetype="dashed", color="dodgerblue3", size = 2) + #H0 median value
  theme_bw() +
  theme(axis.text.x = element_text(size=22, face="bold")) +
  theme(axis.text.y = element_text(size=22, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) +
  scale_x_continuous(name="\nPredicted Probability", limits=c(-0.05, 1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_continuous(name="Count\n", limits=c(0, 278000),
                     breaks=c(0, 50000, 100000, 150000, 200000, 250000),
                     labels=c("0", "50,000", "100,000", "150,000", "200,000", "250,000")) +
  scale_fill_manual(name='', values = c('H3' = "dodgerblue3"),
                    labels = c('')) +
  # theme(aspect.ratio=1) + #aspect ratio expressed as y/x
  # theme(legend.position = c(0.85, 0.60)) +
  theme(legend.position = "none")


setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/Change in Extent")
ggsave("Extent Histogram H3.jpg", plot=extent.hist3, device = "jpeg",
       width=7, height=5, dpi=600)


extent.hist5 <- ggplot() +
  geom_histogram(data=combined, aes(x=H5_median_corrected, fill = 'H5'), color = "grey3", binwidth = 0.05) +
  # geom_histogram(data=combined, aes(x=H3_median_corrected, fill = 'H3'), color = "grey3", binwidth = 0.05) +
  # geom_histogram(data=combined, aes(x=H0_median_corrected, fill = 'H0'), color = "grey3", binwidth = 0.05) +
  geom_vline(xintercept=0.2316998, linetype="dashed", color="dodgerblue3", size = 2) + #H5 median value
  # geom_vline(xintercept=0.766747, linetype="dashed", color="dodgerblue3", size = 2) + #H3 median value
  # geom_vline(xintercept=0.8138023, linetype="dashed", color="dodgerblue3", size = 2) + #H0 median value
  theme_bw() +
  theme(axis.text.x = element_text(size=22, face="bold")) +
  theme(axis.text.y = element_text(size=22, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) +
  scale_x_continuous(name="\nPredicted Probability", limits=c(-0.05, 1.05),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.00),
                     labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_continuous(name="Count\n", limits=c(0, 278000),
                     breaks=c(0, 50000, 100000, 150000, 200000, 250000),
                     labels=c("0", "50,000", "100,000", "150,000", "200,000", "250,000")) +
  scale_fill_manual(name='', values = c('H5' = "dodgerblue3"),
                    labels = c('')) +
  # theme(aspect.ratio=1) + #aspect ratio expressed as y/x
  # theme(legend.position = c(0.85, 0.60)) +
  theme(legend.position = "none")


setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/Change in Extent")
ggsave("Extent Histogram H5.jpg", plot=extent.hist5, device = "jpeg",
       width=7, height=5, dpi=600)





































### Kruskal Wallis and Nemenyi - based on above notes KW does not seem to the the test to use bc the assumption
### of independence between raster cells is violated
# https://rcompanion.org/rcompanion/d_06.html
# Medians and descriptive statistics
library(FSA)
Summarize(Probability ~ Treatment, data = combined2)

# Stacked histograms of values across groups
library(lattice)
histogram(~ Probability | Treatment, data=combined2, layout = c(1, 3))

# Kruskal-Wallis test
kruskal.test(Probability ~ Treatment, data = combined2)
# There is a significant difference in the distributions of values among groups based on the significant Kruskal-Wallis test.
# Only in cases where the distributions in each group are similar can a significan KW test be interpreted as a difference in
# medians. H0 and H3 have bimodal distributions where H5 is slightly bimodal but more zero-inflated.

# Nemenyi test for mulitple comparisons
# Zar, J.H. 2010. Biostatistical Analysis, 5th ed.  Pearson Prentice Hall: Upper Saddle River, NJ.
library(DescTools)
posthoc <- NemenyiTest(x = combined2$Probability,
                       g = combined2$Treatment,
                       dist = "tukey")

posthoc

# Incidentally, I got the same results when I performed Pairwise Mann-Whitney U Tests






