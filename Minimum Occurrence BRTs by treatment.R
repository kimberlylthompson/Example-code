####################################################################################
###### Analysis of winter 2016-17, December - March subnivium presence/absence ######
#####              data using Boosted Regression Trees                        #######
####  FOR DATA WHICH TOOK THE MINIMUM OCCURRENCE OF PRESENCE INSTANCES FROM   #######
####                        ALL TREATMENT PROBES                              #######
#####################################################################################

# Minimum Occurrence: If one probe shows subnivium presence, then that day for that treatment is assigned a 1.

# UpdateD to rerun model diagnostics to make sure I am using the correct settings for the models
# in light of the additional changes we decided to make to the parameters
# 1st that we are not including latitude
# 2nd that we are going to include a site & house-specific density that comes from
# our own measures of depth, and SNODAS's measure of external density. 


########## clean workspace and load required packages ####################
###########################################################################
#####clean workspace to improve efficiency: ###
rm(list = ls() ) 
gc() #releases memory

####### load relevant packages ###
library(lubridate)
library( ggplot2 ) #fancy plots
library( dplyr ) #data manipulation
library( dismo ) #Species distribution analysis package. Package used in Elith's 2017 tutorial
library( gbm) #Boosted regression tree package that is used in Elith's 2008 tutorial
########## end of package loading ###########



##############################################################
########    data loading            ##########################
##############################################################

setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Presence Absence Min Occur Full Data")

presence<-as.data.frame(read.csv(file="02-04-19 Master_min_presence.csv"),
                        colClasses=c("integer", "character", "integer",
                                     rep("character", 2), rep("numeric", 7), rep("character", 2),
                                     rep("numeric", 8))) 
#Fix the dates
presence$Date<-gsub("/", "-", presence$Date)
presence$Date<-parse_date_time(presence$Date, c("%Y-%m-%d", "%m-%d-%Y"), tz = "US/Central")
presence$Date<-as.POSIXct(presence$Date, tz = "US/Central", format = "%Y-%m-%d")


#### Model each treatment separately ####

pres0 <- presence[presence$Treat == "0",]
pres3 <- presence[presence$Treat == "3",]
pres5 <- presence[presence$Treat == "5",]
presext <- presence[presence$Treat == "ext",]

#Add ID column to each dataframe to aid in subsetting
pres0$ID <- seq(from=1, to=1089, by=1)
pres3$ID <- seq(from=1, to=1089, by=1)
pres5$ID <- seq(from=1, to=1089, by=1)
presext$ID <- seq(from=1, to=1089, by=1)


######################################################################################
####                              Model Explanation                               ####
######################################################################################

# Within each treatment dataset, two sets of presence/absence data are necessary,
# one for model training (building) and one for model testing (evaluation).
# Starting with approximately 70% of data for model building and 30% for model testing.

## First, identify the optimal number of trees. ##

## As per the article by Elith and Leathwick ##
## found here, https://cran.r-project.org/web/packages/dismo/vignettes/brt.pdf, look at the ##
## total number of observations you have and the total number of presences, and decide whether ##
## you think there is enough data to model interactions with reasonable complexity. They say that if yes, ##
## then a learning rate of 0.01 could be a reasonable starting point.                       ##

# data = training dataset, gbm.x = indices of names of predictor columns in data
# gbm.y = index of response variable
# family = use Bernoulli for presence absence (binomial)
# tree complexity = set the complexity of individual trees (models interactions --> 2 = two-way; 3 = three-way)
# learning rate = sets the weight applied to individual trees (lower the learning rate, the better performance in model building)
# bag fraction = controls the stochasticity in the model, bag.frac = 1 results in deterministic model, optimal is between 0.5-0.75

# Keep in mind these models are stochastic and so slightly different each time you run them unless you make them deterministic
# by using a bag fraction of 1. 

# Goal is the joint optimization of the number of trees, the learning rate, and tree complexity
# (and minimal predictive deviance)

# IMPORTANT: Splitting the data into training and testing data is necessary in order to more broadly explore
# the optimal combination of number of trees, learning rate, and tree complexity
# Systematically alter tc, lr, and the bag fraction and compare the results.
# Since the external treatment has the most complexity associated with it, do this systematic approach on
# this treatment, and then apply that to the other treatments where things like wind are minimized.





######################################################################################
####  FINDING THE IDEAL SETTINGS FOR THE MODEL USING THE EXTERNAL PROBE DATASET  ####
######################################################################################

#####  Split data into training and testing datasets with approximately 70/30 split

# Total number of sites sampled in presence/absence data:
M <- max( presext$ID )

# Define number of data points to use for testing:
TP <- 326 # ~30% of data points
# Select which rows to use for testing for our presence/absence data:
t.parows <- sample( x = 1:M, size = TP, replace = FALSE )

presext.train <- presext[ -t.parows, ]

# View
head( presext.train ); dim( presext.train )

# Create testing dataset:
presext.test <- presext[ t.parows, ]

#### Systematic testing of different numbers of trees and learning rates

# Create a list of learning rates to test
lr <- c(0.05, 0.01, 0.005, 0.001, 0.0005)
tree.list <- seq(100, 20000, by=100)


################ TC = 1 ###########################

#Create an empty dataframe in which to store the predicted deviances
deviance.tc1<-data.frame(lr05=numeric(200), lr01=numeric(200), lr005=numeric(200), lr001=numeric(200),
                        lr0005=numeric(200), nt=seq(100, 20000, by=100))

# Model each different learning rate and store them in a data frame
for (i in 1:length(lr)) {
  
  mod <- gbm.fixed(data=presext.train, gbm.x = c(6, 7, 10:13), gbm.y = 1,
                   family = "bernoulli", tree.complexity = 1, learning.rate = lr[i],
                   n.trees = 20000)
  
  # Create matrix of predictions, each column = predictions from the mode - for example, the predictions
  # in column 5 are for tree.list[5]=500 trees
  pred <- gbm :: predict.gbm(mod, presext.test, n.trees = tree.list, type="response")
  
  # Calculate deviance of all of these results and store them in a dataframe:
  for (j in 1:length(deviance.tc1[,i])) {
    deviance.tc1[,i][j] <- calc.deviance(presext.test$Presence, pred[,j], calc.mean = T)
  }
}

min(deviance.tc1, na.rm=TRUE)
which(deviance.tc1$lr05 == min(deviance.tc1, na.rm=TRUE))
deviance.tc1[6,]
#Minimum is 0.790  at 600 trees at lr of 0.05

# Graph the results of the learning rate comparison for tc=1
tc1 <- ggplot() +
  geom_line(data=deviance.tc1, aes(x=nt, y=lr05, linetype="0.05"), size=1) +
  geom_line(data=deviance.tc1, aes(x=nt, y=lr01, linetype="0.01"), size=1) +
  geom_line(data=deviance.tc1, aes(x=nt, y=lr005, linetype="0.005"), size=1) +
  geom_line(data=deviance.tc1, aes(x=nt, y=lr001, linetype="0.001"), size=1) +
  geom_line(data=deviance.tc1, aes(x=nt, y=lr0005, linetype="0.0005"), size=1) +
  geom_hline(yintercept=min(deviance.tc1[,c(1:5)]), linetype="solid", color="red") +
  geom_vline(xintercept=600, linetype="solid", color="green") +
  #geom_errorbar(width=0.2, size=2) + geom_point(size=8) +
  theme_bw() +
  theme(axis.text.x = element_text(size=20, face="bold")) +
  theme(axis.text.y = element_text(size=20, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) +
  scale_x_continuous(name="\nNumber of Trees", limits=c(0, 20000),
                     breaks=c(0, 5000, 10000, 15000, 20000)) +
  scale_y_continuous(name="Predictive Deviance\n") +
  scale_linetype_manual(name='', values = c('0.05' = "twodash", 
                                            '0.01' = "solid",
                                            '0.005' = "longdash",
                                            '0.001' = "dotted",
                                            '0.0005' = "dotdash"),
                        labels = c("0.0005", "0.001", "0.005", "0.01", "0.05")) +
  theme(aspect.ratio=1) + #aspect ratio expressed as y/x
  # theme(legend.position = c(0.85, 0.60)) +
  theme(legend.position = "none") +
  theme(legend.text = element_text(size=14)) +
  guides(colour = guide_legend(override.aes = (list(size=3)))) +
  annotate("text", x = 300, y = 1.18, label = "(a)",
           size=7) +
  annotate("text", x = 3500, y = 1.18, label = "tc = 1",
          fontface="italic", size=7) +
  annotate("text", x = 14500, y = 1.18, label = "min=0.79, lr=0.05",
           fontface="italic", size=7) +
  annotate("text", x = 17500, y = 1.14, label = "nt=600",
           fontface="italic", size=7)
  

setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/BRT/02-04-19 Minimum Occurrence")
ggsave("Ext - Deviance tc1 nt20000.jpg", plot=tc1, device = "jpeg",
       width=7, height=5, dpi=300)



################ TC = 2 ###########################

#Create an empty dataframe in which to store the predicted deviances
deviance.tc2<-data.frame(lr05=numeric(200), lr01=numeric(200), lr005=numeric(200), lr001=numeric(200),
                         lr0005=numeric(200), nt=seq(100, 20000, by=100))

# Model each different learning rate and store them in a data frame
for (i in 1:length(lr)) {
  
  mod <- gbm.fixed(data=presext.train, gbm.x = c(6, 7, 10:13), gbm.y = 1,
                   family = "bernoulli", tree.complexity = 2, learning.rate = lr[i],
                   n.trees = 20000)
  
  # Create matrix of predictions, each column = predictions from the mode - for example, the predictions
  # in column 5 are for tree.list[5]=500 trees
  pred <- gbm :: predict.gbm(mod, presext.test, n.trees = tree.list, type="response")
  
  # Calculate deviance of all of these results and store them in a dataframe:
  for (j in 1:length(deviance.tc2[,i])) {
    deviance.tc2[,i][j] <- calc.deviance(presext.test$Presence, pred[,j], calc.mean = T)
  }
}

min(deviance.tc2, na.rm=TRUE)
which(deviance.tc2$lr01 == min(deviance.tc2, na.rm=TRUE))
deviance.tc2[31,]
#Minimum is 0.70 at 3100 trees at lr of 0.01

# Graph the results of the learning rate comparison for tc=1
tc2 <- ggplot() +
  geom_line(data=deviance.tc2, aes(x=nt, y=lr05, linetype="0.05"), size=1) +
  geom_line(data=deviance.tc2, aes(x=nt, y=lr01, linetype="0.01"), size=1) +
  geom_line(data=deviance.tc2, aes(x=nt, y=lr005, linetype="0.005"), size=1) +
  geom_line(data=deviance.tc2, aes(x=nt, y=lr001, linetype="0.001"), size=1) +
  geom_line(data=deviance.tc2, aes(x=nt, y=lr0005, linetype="0.0005"), size=1) +
  geom_hline(yintercept=min(deviance.tc2[,c(1:5)], na.rm=TRUE), linetype="solid", color="red") +
  geom_vline(xintercept=3100, linetype="solid", color="green") +
  #geom_errorbar(width=0.2, size=2) + geom_point(size=8) +
  theme_bw() +
  theme(axis.text.x = element_text(size=20, face="bold")) +
  theme(axis.text.y = element_text(size=20, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) +
  scale_x_continuous(name="\nNumber of Trees", limits=c(0, 20000),
                     breaks=c(0, 5000, 10000, 15000, 20000)) +
  scale_y_continuous(name="Predictive Deviance\n", breaks=c(0.9, 1.2, 1.5)) +
  scale_linetype_manual(name='', values = c('0.05' = "twodash", 
                                            '0.01' = "solid",
                                            '0.005' = "longdash",
                                            '0.001' = "dotted",
                                            '0.0005' = "dotdash"),
                        labels = c("0.0005", "0.001", "0.005", "0.01", "0.05")) +
  theme(aspect.ratio=1) + #aspect ratio expressed as y/x
  theme(legend.position = c(0.85, 0.50)) +
  # theme(legend.position = "none") +
  theme(legend.text = element_text(size=14)) +
  guides(colour = guide_legend(override.aes = (list(size=3)))) +
  annotate("text", x = 300, y = 1.75, label = "(b)",
           size=7) +
  annotate("text", x = 3500, y = 1.75, label = "tc = 2",
           fontface="italic", size=7) +
  annotate("text", x = 14500, y = 1.75, label = "min=0.70, lr=0.01",
           fontface="italic", size=7) +
  annotate("text", x = 17300, y = 1.63, label = "nt=3100",
           fontface="italic", size=7)


setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/BRT/02-04-19 Minimum Occurrence")
ggsave("Ext - Deviance tc2 nt20000.jpg", plot=tc2, device = "jpeg",
       width=7, height=5, dpi=300)



################ TC = 3 ###########################

#Create an empty dataframe in which to store the predicted deviances
deviance.tc3<-data.frame(lr05=numeric(200), lr01=numeric(200), lr005=numeric(200), lr001=numeric(200),
                         lr0005=numeric(200), nt=seq(100, 20000, by=100))

# Model each different learning rate and store them in a data frame
for (i in 1:length(lr)) {
  
  mod <- gbm.fixed(data=presext.train, gbm.x = c(6, 7, 10:13), gbm.y = 1,
                   family = "bernoulli", tree.complexity = 3, learning.rate = lr[i],
                   n.trees = 20000)
  
  # Create matrix of predictions, each column = predictions from the mode - for example, the predictions
  # in column 5 are for tree.list[5]=500 trees
  pred <- gbm :: predict.gbm(mod, presext.test, n.trees = tree.list, type="response")
  
  # Calculate deviance of all of these results and store them in a dataframe:
  for (j in 1:length(deviance.tc3[,i])) {
    deviance.tc3[,i][j] <- calc.deviance(presext.test$Presence, pred[,j], calc.mean = T)
  }
}

min(deviance.tc3, na.rm=TRUE)
which(deviance.tc3$lr01 == min(deviance.tc3, na.rm=TRUE))
deviance.tc3[23,]
#Minimum is 0.65 at 2300 trees at lr of 0.01

# Graph the results of the learning rate comparison for tc=3
tc3 <- ggplot() +
  geom_line(data=deviance.tc3, aes(x=nt, y=lr05, linetype="0.05"), size=1) +
  geom_line(data=deviance.tc3, aes(x=nt, y=lr01, linetype="0.01"), size=1) +
  geom_line(data=deviance.tc3, aes(x=nt, y=lr005, linetype="0.005"), size=1) +
  geom_line(data=deviance.tc3, aes(x=nt, y=lr001, linetype="0.001"), size=1) +
  geom_line(data=deviance.tc3, aes(x=nt, y=lr0005, linetype="0.0005"), size=1) +
  geom_hline(yintercept=min(deviance.tc3[,c(1:5)], na.rm=TRUE), linetype="solid", color="red") +
  geom_vline(xintercept=2300, linetype="solid", color="green") +
  #geom_errorbar(width=0.2, size=2) + geom_point(size=8) +
  theme_bw() +
  theme(axis.text.x = element_text(size=20, face="bold")) +
  theme(axis.text.y = element_text(size=20, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) +
  scale_x_continuous(name="\nNumber of Trees", limits=c(0, 20000),
                     breaks=c(0, 5000, 10000, 15000, 20000)) +
  scale_y_continuous(name="Predictive Deviance\n") +
  scale_linetype_manual(name='', values = c('0.05' = "twodash", 
                                            '0.01' = "solid",
                                            '0.005' = "longdash",
                                            '0.001' = "dotted",
                                            '0.0005' = "dotdash"),
                        labels = c("0.0005", "0.001", "0.005", "0.01", "0.05")) +
  theme(aspect.ratio=1) + #aspect ratio expressed as y/x
  # theme(legend.position = c(0.85, 0.60)) +
  theme(legend.position = "none") +
  theme(legend.text = element_text(size=14)) +
  guides(colour = guide_legend(override.aes = (list(size=3)))) +
  annotate("text", x = 300, y = 2, label = "(c)",
           size=7) +
  annotate("text", x = 3500, y = 2, label = "tc = 3",
           fontface="italic", size=7) +
  annotate("text", x = 14500, y = 2, label = "min=0.65, lr=0.01",
           fontface="italic", size=7) +
  annotate("text", x = 17300, y = 1.85, label = "nt=2300",
           fontface="italic", size=7)


setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/BRT/02-04-19 Minimum Occurrence")
ggsave("Ext - Deviance tc3 nt20000.jpg", plot=tc3, device = "jpeg",
       width=7, height=5, dpi=300)





################ TC = 4 ###########################

#Create an empty dataframe in which to store the predicted deviances
deviance.tc4<-data.frame(lr05=numeric(200), lr01=numeric(200), lr005=numeric(200), lr001=numeric(200),
                         lr0005=numeric(200), nt=seq(100, 20000, by=100))

# Model each different learning rate and store them in a data frame
for (i in 1:length(lr)) {
  
  mod <- gbm.fixed(data=presext.train, gbm.x = c(6, 7, 10:13), gbm.y = 1,
                   family = "bernoulli", tree.complexity = 4, learning.rate = lr[i],
                   n.trees = 20000)
  
  # Create matrix of predictions, each column = predictions from the mode - for example, the predictions
  # in column 5 are for tree.list[5]=500 trees
  pred <- gbm :: predict.gbm(mod, presext.test, n.trees = tree.list, type="response")
  
  # Calculate deviance of all of these results and store them in a dataframe:
  for (j in 1:length(deviance.tc4[,i])) {
    deviance.tc4[,i][j] <- calc.deviance(presext.test$Presence, pred[,j], calc.mean = T)
  }
}

min(deviance.tc4, na.rm=TRUE)
which(deviance.tc4$lr01 == min(deviance.tc4, na.rm=TRUE))
deviance.tc4[19,]
#Minimum is 0.62 at 1900 trees at lr of 0.01

# Graph the results of the learning rate comparison for tc=1
tc4 <- ggplot() +
  geom_line(data=deviance.tc4, aes(x=nt, y=lr05, linetype="0.05"), size=1) +
  geom_line(data=deviance.tc4, aes(x=nt, y=lr01, linetype="0.01"), size=1) +
  geom_line(data=deviance.tc4, aes(x=nt, y=lr005, linetype="0.005"), size=1) +
  geom_line(data=deviance.tc4, aes(x=nt, y=lr001, linetype="0.001"), size=1) +
  geom_line(data=deviance.tc4, aes(x=nt, y=lr0005, linetype="0.0005"), size=1) +
  geom_hline(yintercept=min(deviance.tc4[,c(1:5)], na.rm=TRUE), linetype="solid", color="red") +
  geom_vline(xintercept=1900, linetype="solid", color="green") +
  #geom_errorbar(width=0.2, size=2) + geom_point(size=8) +
  theme_bw() +
  theme(axis.text.x = element_text(size=20, face="bold")) +
  theme(axis.text.y = element_text(size=20, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) +
  scale_x_continuous(name="\nNumber of Trees", limits=c(0, 20000),
                     breaks=c(0, 5000, 10000, 15000, 20000)) +
  scale_y_continuous(name="Predictive Deviance\n") +
  scale_linetype_manual(name='', values = c('0.05' = "twodash", 
                                            '0.01' = "solid",
                                            '0.005' = "longdash",
                                            '0.001' = "dotted",
                                            '0.0005' = "dotdash"),
                        labels = c("0.0005", "0.001", "0.005", "0.01", "0.05")) +
  theme(aspect.ratio=1) + #aspect ratio expressed as y/x
  # theme(legend.position = c(0.85, 0.60)) +
  theme(legend.position = "none") +
  theme(legend.text = element_text(size=14)) +
  guides(colour = guide_legend(override.aes = (list(size=3)))) +
  annotate("text", x = 300, y = 1.95, label = "(d)",
           size=7) +
  annotate("text", x = 3500, y = 1.95, label = "tc = 4",
           fontface="italic", size=7) +
  annotate("text", x = 14500, y = 1.95, label = "min=0.62, lr=0.01",
           fontface="italic", size=7) +
  annotate("text", x = 17300, y = 1.80, label = "nt=1900",
           fontface="italic", size=7)


setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/BRT/02-04-19 Minimum Occurrence")
ggsave("Ext - Deviance tc4 nt20000.jpg", plot=tc4, device = "jpeg",
       width=7, height=5, dpi=300)


######## Based on the above systematic alteration of the learning rates and tree ######
######## complexity, the optimal combination of settings is tc=3, lr=0.005, and  ######
######## nt= 5200 for the EXTERNAL PROBES. TC4 also has a lower deviance and     ######
######## it does have at least 1000 trees, but thinking about a tree complexity #####
######## meaning essentially a 4-way interaction seems excessive.


######################### BAG FRACTION = STOCHASTICITY ###############################

# Now vary the bag fraction 
bf<-c(0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75)

#for the first potential tc and lr combination (tc=3 and lr=0.005)
#Create an empty dataframe in which to store the predicted deviances
deviance.bf<-data.frame(bf25=numeric(200), bf30=numeric(200), bf35=numeric(200), bf40=numeric(200),
                        bf45=numeric(200), bf50=numeric(200), bf55=numeric(200), bf60=numeric(200),
                        bf65=numeric(200), bf70=numeric(200), bf75=numeric(200), nt=seq(100, 20000, by=100))

# Model each different learning rate and store them in a data frame
for (i in 1:length(bf)) {
  
  mod <- gbm.fixed(data=presext.train, gbm.x = c(6, 7, 10:13), gbm.y = 1,
                   family = "bernoulli", tree.complexity = 3, learning.rate = 0.005,
                   n.trees = 20000, bag.fraction = bf[i])
  
  # Create matrix of predictions, each column = predictions from the mode
  pred <- gbm :: predict.gbm(mod, presext.test, n.trees = tree.list, type="response")
  
  # Calculate deviance of all of these results and store them in a dataframe:
  for (j in 1:length(deviance.bf[,i])) {
    deviance.bf[,i][j] <- calc.deviance(presext.test$Presence, pred[,j], calc.mean = T)
  }
}

min(deviance.bf, na.rm=TRUE)
which(deviance.bf[,5] == min(deviance.bf, na.rm=TRUE))
deviance.bf[54,]

#Minimum is 0.654 at 4600 trees at bf of 0.45 - 02/05/19
# Rerunning plots for manuscript - minimum is 0.65 at 5400 trees at bf of 0.45

# Graph the results of the bag fraction comparison
bf <- ggplot() +
  geom_line(data=deviance.bf, aes(x=nt, y=bf25, color="bf25"), size=1) +
  geom_line(data=deviance.bf, aes(x=nt, y=bf30, color="bf30"), size=1) +
  geom_line(data=deviance.bf, aes(x=nt, y=bf35, color="bf35"), size=1) +
  geom_line(data=deviance.bf, aes(x=nt, y=bf40, color="bf40"), size=1) +
  geom_line(data=deviance.bf, aes(x=nt, y=bf45, color="bf45"), size=1) +
  geom_line(data=deviance.bf, aes(x=nt, y=bf50, color="bf50"), size=1) +
  geom_line(data=deviance.bf, aes(x=nt, y=bf55, color="bf55"), size=1) +
  geom_line(data=deviance.bf, aes(x=nt, y=bf60, color="bf60"), size=1) +
  geom_line(data=deviance.bf, aes(x=nt, y=bf65, color="bf65"), size=1) +
  geom_line(data=deviance.bf, aes(x=nt, y=bf70, color="bf70"), size=1) +
  geom_line(data=deviance.bf, aes(x=nt, y=bf75, color="bf75"), size=1) +
  geom_hline(yintercept=min(deviance.bf[,c(1:11)], na.rm=TRUE), linetype="solid", color="red") +
  geom_vline(xintercept=5400, linetype="solid", color="green") +
  #geom_errorbar(width=0.2, size=2) + geom_point(size=8) +
  theme_bw() +
  theme(axis.text.x = element_text(size=20, face="bold")) +
  theme(axis.text.y = element_text(size=20, face="bold")) +
  theme(axis.title.x = element_text(size=22, face="bold", color="gray30")) +
  theme(axis.title.y = element_text(size=22, face="bold", color="gray30")) +
  scale_x_continuous(name="\nNumber of Trees", limits=c(0, 20000),
                     breaks=c(0, 5000, 10000, 15000, 20000)) +
  scale_y_continuous(name="Predictive Deviance\n") +
  scale_colour_manual(name='', values = c('bf25' = "red", 
                                          'bf30' = "orange",
                                          'bf35' = "yellow",
                                          'bf40' = "green",
                                          'bf45' = "blue",
                                          'bf50' = "purple",
                                          'bf55' = "darkgreen",
                                          'bf60' = "turquoise2",
                                          'bf65' = "brown3",
                                          'bf70' = "darkgoldenrod3",
                                          'bf75' = "hotpink"), 
                      labels = c("0.25", "0.30", "0.35", "0.40", "0.45", "0.50",
                                 "0.55", "0.60", "0.65", "0.70", "0.75")) +
  theme(aspect.ratio=1) + #aspect ratio expressed as y/x
  # theme(legend.position = c(0.85, 0.60)) +
  # theme(legend.position = "none") +
  theme(legend.text = element_text(size=14)) +
  guides(colour = guide_legend(override.aes = (list(size=3)))) +
  # annotate("text", x = 2200, y = 1, label = "bf = 0.45",
  #          fontface="italic", size=7) +
  # annotate("text", x = 12500, y = 1, label = "min=0.65, lr=0.01",
  #          fontface="italic", size=7) +
  # annotate("text", x = 12500, y = 0.9, label = "nt=5400, tc=3",
  #          fontface="italic", size=7)


setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/BRT/02-04-19 Minimum Occurrence")
ggsave("Ext - Deviance BF nt20000.jpg", plot=bf, device = "jpeg",
       width=7, height=5, dpi=300)



# Based on the above, ideal settings are as follows (because deviance is slightly minimized with a tree complexity of 3):
# learning rate = 0.005
# tree complexity = 3
# bag fraction = 0.45



######################################################################################
####                              External Model                                 ####
######################################################################################

#Rerun test and training separation so that new random sets (different from those used
#above in the optimal settings identification) can be generated

#####  Split data into training and testing datasets with approximately 70/30 split

# Total number of sites sampled in presence/absence data:
M <- max( presext$ID )

# Define number of data points to use for testing:
TP <- 326 # ~30% of data points
# Select which rows to use for testing for our presence/absence data:
t.parows <- sample( x = 1:M, size = TP, replace = FALSE )

presext.train <- presext[ -t.parows, ]

# View
head( presext.train ); dim( presext.train )

# Create testing dataset:
presext.test <- presext[ t.parows, ]

########## Model

#BRT model with optimal settings identified from first part of code

modext <- gbm.step(data=presext.train, gbm.x = c(6, 7, 10, 11, 13, 22), gbm.y = 1,
                 family = "bernoulli", tree.complexity = 3, learning.rate = 0.005,
                 max.trees = 7500, bag.fraction = 0.45)

# Be careful of partial dependence plots (plotting the effect of each single variable)
# if predictors are correlated or if there are interactions

###### Identifying interactions ######
find.int <- dismo::gbm.interactions(modext)
find.int$interactions
find.int$rank.list

gbm.perspec(modext, 1, 6)

gbm.perspec(modext, 5, 3)

#Variable importance
modext$contributions

#partial dependence plots
gbm.plot(modext)

### Predicting to testing data - MODEL evaluation ###

preds <- gbm::predict.gbm(modext, presext.test, n.trees=modext$gbm.call$best.trees, type="response")
calc.deviance(obs=presext.test$Presence, pred=preds, calc.mean=TRUE)

d <- cbind(presext.test$Presence, preds)
pres <- d[d[,1]==1, 2]
abs <- d[d[,1]==0, 2]
e <- evaluate(p=pres, a=abs)
e

plot(e, 'ROC')

# Save this graph in case I need to present them in supplemental materials
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/BRT/AUC Plots/Minimum Occurrence")
jpeg('External Model AUC Curve.jpg', width=7, height=5, units= 'in', res=300)
plot(e, 'ROC')
dev.off()


modext$cv.statistics

# In order to make prettier graphs of any interactions and/or partial depency plots, make a new dataframe that has the 
# training data, and the fitted values. Should then be able to use them with ggplot.
# Add columns for fitted, fitted.vars, and residuals

ext <- cbind(presext.train, modext$fitted, modext$fitted.vars, modext$residuals, modext$fit)
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Modeling Results/Fitted values for graphing/Minimum Occurrence")
write.csv(ext, "External Model Fitted Values.csv", row.names = FALSE)






######################################################################################
####                              House 0 Model                                 ####
######################################################################################

#####  Split data into training and testing datasets with approximately 70/30 split

# Total number of sites sampled in presence/absence data:
M <- max( pres0$ID )

# Define number of data points to use for testing:
TP <- 326 # ~30% of data points
# Select which rows to use for testing for our presence/absence data:
t.parows <- sample( x = 1:M, size = TP, replace = FALSE )

pres0.train <- pres0[ -t.parows, ]

# View
head( pres0.train ); dim( pres0.train )

# Create testing dataset:
pres0.test <- pres0[ t.parows, ]

########## Model

#BRT model with optimal settings identified from first part of code

mod0 <- gbm.step(data=pres0.train, gbm.x = c(6, 7, 10, 11, 13, 22), gbm.y = 1,
                   family = "bernoulli", tree.complexity = 3, learning.rate = 0.005,
                   max.trees = 7500, bag.fraction = 0.45)

# Be careful of partial dependence plots (plotting the effect of each single variable)
# if predictors are correlated or if there are interactions

###### Identifying interactions ######
find.int <- dismo::gbm.interactions(mod0)
find.int$interactions
find.int$rank.list

gbm.perspec(mod0, 5, 3)

#Variable importance
mod0$contributions

#partial dependence plots
gbm.plot(mod0)

### Predicting to testing data - MODEL evaluation ###

preds <- gbm::predict.gbm(mod0, pres0.test, n.trees=mod0$gbm.call$best.trees, type="response")
calc.deviance(obs=pres0.test$Presence, pred=preds, calc.mean=TRUE)

d <- cbind(pres0.test$Presence, preds)
pres <- d[d[,1]==1, 2]
abs <- d[d[,1]==0, 2]
e <- evaluate(p=pres, a=abs)
e

plot(e, 'ROC')

# Save this graph in case I need to present them in supplemental materials
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/BRT/AUC Plots/Minimum Occurrence")
jpeg('House 0 Model AUC Curve.jpg', width=7, height=5, units= 'in', res=300)
plot(e, 'ROC')
dev.off()


mod0$cv.statistics

# In order to make prettier graphs of any interactions and/or partial depency plots, make a new dataframe that has the 
# training data, and the fitted values. Should then be able to use them with ggplot.
# Add columns for fitted, fitted.vars, and residuals

H0 <- cbind(pres0.train, mod0$fitted, mod0$fitted.vars, mod0$residuals, mod0$fit)
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Modeling Results/Fitted values for graphing/Minimum Occurrence")
write.csv(H0, "House 0 Model Fitted Values.csv", row.names = FALSE)


######################################################################################
####                              House 3 Model                                 ####
######################################################################################


#####  Split data into training and testing datasets with approximately 70/30 split

# Total number of sites sampled in presence/absence data:
M <- max( pres3$ID )

# Define number of data points to use for testing:
TP <- 326 # ~30% of data points
# Select which rows to use for testing for our presence/absence data:
t.parows <- sample( x = 1:M, size = TP, replace = FALSE )

pres3.train <- pres3[ -t.parows, ]

# View
head( pres3.train ); dim( pres3.train )

# Create testing dataset:
pres3.test <- pres3[ t.parows, ]

########## Model

#BRT model with optimal settings identified from first part of code

mod3 <- gbm.step(data=pres3.train, gbm.x = c(6, 7, 10, 11, 13, 22), gbm.y = 1,
                 family = "bernoulli", tree.complexity = 3, learning.rate = 0.005,
                 max.trees = 7500, bag.fraction = 0.45)

# Be careful of partial dependence plots (plotting the effect of each single variable)
# if predictors are correlated or if there are interactions

###### Identifying interactions ######
find.int <- dismo::gbm.interactions(mod3)
find.int$interactions
find.int$rank.list

gbm.perspec(mod3, 5, 3)
gbm.perspec(mod3, 6, 3)

#Variable importance
mod3$contributions

#partial dependence plots
gbm.plot(mod3)

### Predicting to testing data - MODEL evaluation ###

preds <- gbm::predict.gbm(mod3, pres3.test, n.trees=mod3$gbm.call$best.trees, type="response")
calc.deviance(obs=pres3.test$Presence, pred=preds, calc.mean=TRUE)

d <- cbind(pres3.test$Presence, preds)
pres <- d[d[,1]==1, 2]
abs <- d[d[,1]==0, 2]
e <- evaluate(p=pres, a=abs)
e

plot(e, 'ROC')

# Save this graph in case I need to present them in supplemental materials
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/BRT/AUC Plots/Minimum Occurrence")
jpeg('House 3 Model AUC Curve.jpg', width=7, height=5, units= 'in', res=300)
plot(e, 'ROC')
dev.off()


mod3$cv.statistics

# In order to make prettier graphs of any interactions and/or partial depency plots, make a new dataframe that has the 
# training data, and the fitted values. Should then be able to use them with ggplot.
# Add columns for fitted, fitted.vars, and residuals

H3 <- cbind(pres3.train, mod3$fitted, mod3$fitted.vars, mod3$residuals, mod3$fit)
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Modeling Results/Fitted values for graphing/Minimum Occurrence")
write.csv(H3, "House 3 Model Fitted Values.csv", row.names = FALSE)



######################################################################################
####                              House 5 Model                                 ####
######################################################################################

#####  Split data into training and testing datasets with approximately 70/30 split

# Total number of sites sampled in presence/absence data:
M <- max( pres5$ID )

# Define number of data points to use for testing:
TP <- 326 # ~30% of data points
# Select which rows to use for testing for our presence/absence data:
t.parows <- sample( x = 1:M, size = TP, replace = FALSE )

pres5.train <- pres5[ -t.parows, ]

# View
head( pres5.train ); dim( pres5.train )

# Create testing dataset:
pres5.test <- pres5[ t.parows, ]

########## Model

#BRT model with optimal settings identified from first part of code

mod5 <- gbm.step(data=pres5.train, gbm.x = c(6, 7, 10, 11, 13, 22), gbm.y = 1,
                 family = "bernoulli", tree.complexity = 3, learning.rate = 0.005,
                 max.trees = 7500, bag.fraction = 0.45)

# Be careful of partial dependence plots (plotting the effect of each single variable)
# if predictors are correlated or if there are interactions

###### Identifying interactions ######
find.int <- dismo::gbm.interactions(mod5)
find.int$interactions
find.int$rank.list

gbm.perspec(mod5, 3, 6)
gbm.perspec(mod5, 6, 3)


#Variable importance
mod5$contributions

#partial dependence plots
gbm.plot(mod5)

### Predicting to testing data - MODEL evaluation ###

preds <- gbm::predict.gbm(mod5, pres5.test, n.trees=mod5$gbm.call$best.trees, type="response")
calc.deviance(obs=pres5.test$Presence, pred=preds, calc.mean=TRUE)

d <- cbind(pres5.test$Presence, preds)
pres <- d[d[,1]==1, 2]
abs <- d[d[,1]==0, 2]
e <- evaluate(p=pres, a=abs)
e

plot(e, 'ROC')

# Save this graph in case I need to present them in supplemental materials
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Graphs/BRT/AUC Plots/Minimum Occurrence")
jpeg('House 5 Model AUC Curve.jpg', width=7, height=5, units= 'in', res=300)
plot(e, 'ROC')
dev.off()


mod5$cv.statistics

# In order to make prettier graphs of any interactions and/or partial depency plots, make a new dataframe that has the 
# training data, and the fitted values. Should then be able to use them with ggplot.
# Add columns for fitted, fitted.vars, and residuals

H5 <- cbind(pres5.train, mod5$fitted, mod5$fitted.vars, mod5$residuals, mod5$fit)
setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Modeling Results/Fitted values for graphing/Minimum Occurrence")
write.csv(H5, "House 5 Model Fitted Values.csv", row.names = FALSE)


########################################################################################
##### SAVE MODEL OBJECTS SO THAT THEY CAN BE RELOADED FOR PREDICTIONS LATER ############
########################################################################################

setwd("L:/LabMemberFolders/KimberlyThompson/Ch 2 Climate/Analysis/Modeling Results/Model Objects/Minimum Occurrence")
saveRDS(modext, "External_BRT_Model.rds")
saveRDS(mod0, "House0_BRT_Model.rds")
saveRDS(mod3, "House3_BRT_Model.rds")
saveRDS(mod5, "House5_BRT_Model.rds")


