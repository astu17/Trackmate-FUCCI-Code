# Select proper folder ; this should be a folder within a major experiment-defining folder because # this is where the ensuing code will drop all the graphs
setwd("/Users/YOU/INSERT/PATH/HERE")

# You'll need the following theme to make nice black background single-track line graphs. They look pretty solid!
# Use ggplot2 themes to mess around with the look however you'd like :) 

library(ggplot2)
library(gridExtra)
theme_black = function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
    )
}

#Now, load your libraries
library(tidyr)
library(tidyverse)
library(tidyselect)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(dplyr)
library(dbplyr)
library(lubridate)
library(ggpubr)
library(sqldf)
library(cowplot)
library(kernlab)
library(e1071)
library(data.table)
library(ggbiplot)
library(ggmap)
library(broom)
library(lubridate)
library(stringr)

# read in all files
# change input file as needed

# Get all treatment conditions into individual data tables. Example dataframes are MRC5FDMSO1 # and MRC5FATMi1

dfMRC5FDMSO1 <- fread("/YOUR/PATH/HERE/TrackMate DMSO File in a simple txt form.txt", header = T)
dfMRC5FATMi1 <- fread("/YOUR/PATH/HERE/TrackMate ATMi File in a simple txt form.txt", header = T)

#If there are some tracks that you know will mess with the data (e.g., high autofluoresence) you can exlude them with the function below.
#Remove the # from the code to make the line functional
#dfMRC5FDMSO1 <- subset(dfMRC5FDMSO1,!(TRACK_ID %in% c('BAD-TRACK-ID-HERE')))

# Put all input dataframes into list; in this case, there is only one dataframe
df_list <- list(dfMRC5FDMSO1,dfMRC5FATMi1)

len_list <- length(df_list)

# change column names for all dfs
for (i in 1:len_list){
  colnames(df_list[[i]])[18] <- 'mCherry'
  colnames(df_list[[i]])[12] <- 'mVenus'
}

#Re-set to working directory where you want the graphs to end up
setwd("/YOUR/DESIRED/PATH/HERE")

# define condition. "[SPACE][TREATMENT][TechnicalReplicateNumber]"
condition = " Low O2 MRC5 DMSO Proliferating Tracks "

# Create a data frame with the file of choice (specify which you'd like to anaylyze)
tdf <- dfMRC5FDMSO1

# Re-name the appropriate columns after the fluoresence intensity their values represent.
# Make sure these columns align correctly. Various Trackmate updates may mess with the spacing a little bit!
# So be careful! This is saying that column 18 is the mCherry values and column 12 is the mVenus values 
colnames(tdf)[18] = 'mCherry'
colnames(tdf)[12] = 'mVenus'

#get tracks of working dataframe
tracks <- as.list(unique(tdf$TRACK_ID))

## Set up the requirements to graph the raw values

# Create a pared-down data table for graphing the raw values; only includes the raw values, Track ID, Frame, mCherry values and mVenus values
# The gather function is saying "ignore the FRAME and TRACK_ID columns", create a Key entitled fluroesence, underwhich the remaining
# column titles will be listed. The Y Value (that is, the mCherry and mVenus values) will be "Fluoresence_Intensity" on the graph.

raw <- tdf %>%
  select(TRACK_ID, FRAME, mCherry, mVenus) %>%
  gather(key = "Fluoresence", value = "Fluoresence_Intensity", -FRAME, -TRACK_ID)

## Set up for Z score and fluoroscore

#combine all values across conditions for each channel 

mcherry_all <- numeric()
mvenus_all <- numeric()
for (i in df_list){
  #put current iteration mCherry values into a lil list
  mcherry_all <- c(mcherry_all, i$mCherry)
  #put current iteration mVenus values into a lil list
  mvenus_all <- c(mvenus_all, i$mVenus)
}


# check distributions
hist(mcherry_all)
hist(mvenus_all)

# sort by value
mcherry_all <- sort(mcherry_all)
mvenus_all <- sort(mvenus_all)

# Set a lower bound percentage; you will take be scoring intensities based off the amount of Standard Deviations above this value.
# Specifically, the Standard Deviations are the those of the values in this bottom percentage
PercentCutoff <- 5

# find number of measurements equal to 5% (or whatever value is set as "PercentCutoff") of total values.
smallestcohort <- (length(mcherry_all)/100)*PercentCutoff

# find maximum measurement for the bottom 5 percent 
mcherry_min <- mcherry_all[1:smallestcohort]
mvenus_min <- mvenus_all[1:smallestcohort]

#convert to scores as defined by number of sd from the maximum value of the bottom percentage cohort (defined in the previous 2 steps), and add those scores as new columns
# If you would rather alter the distrubition of values (e.g. making it easier to see on graphs), you can use the #-out code
#tdf <- tdf %>% mutate(mcherry_fluorscore = (30 + ((mCherry-max(mcherry_min))/sd(mcherry_min))),
#                             mvenus_fluorscore = (30 + (mVenus-max(mvenus_min))/sd(mvenus_min)))

tdf <- tdf %>% mutate(mcherry_fluorscore = (((mCherry-max(mcherry_min))/sd(mcherry_min))),
                      mvenus_fluorscore = ((mVenus-max(mvenus_min))/sd(mvenus_min)))

#smooth out the fluorscores: smoothing needs to occur on a TRACK_ID basis
# shoutout https://philchodrow.github.io/cos_2018/4_advanced_topics/notes.html for the help
# Warnings don't matter too much in the final output, as long as it’s not an error. There
# are also functions that can help find the "optimal" smoothing

# import relevant libraries
library(tidyr)
library(purrr)

# subset tdf to get two new dataframes (one for mCherry, one for mVenus) 
# these contain only the track ID and the fluoroscore (or z score, or wutever)

mvenus_to_nest <- subset(tdf, select = c(TRACK_ID, mvenus_fluorscore, FRAME))

mcherry_to_nest <- subset(tdf, select = c(TRACK_ID, mcherry_fluorscore, FRAME))



# function for LOESS
my_loess_mvenus <- function(data, span){
  loess(mvenus_fluorscore ~ FRAME,  
        data = data, 
        na.action = na.exclude, 
        span = span)
}


my_loess_mcherry <- function(data, span){
  loess(mcherry_fluorscore ~ FRAME,  
        data = data, 
        na.action = na.exclude, 
        span = span)
}


# fit model/smooth. Define the amount of smoothing aaaas youuuu wiiiish
vspanbig <- 0.32
vspanmed <- 0.19
vspansmall <- 0.158


smooth_mvenus1 <- mvenus_to_nest %>% 
  nest(-TRACK_ID) %>% 
  mutate(model = map(data, my_loess_mvenus, span = vspanbig), preds = map(model, augment)) %>% 
  unnest(preds)

smooth_mvenus2 <- mvenus_to_nest %>% 
  nest(-TRACK_ID) %>% 
  mutate(model = map(data, my_loess_mvenus, span = vspanmed), preds = map(model, augment)) %>% 
  unnest(preds)

smooth_mvenus3 <- mvenus_to_nest %>% 
  nest(-TRACK_ID) %>% 
  mutate(model = map(data, my_loess_mvenus, span = vspansmall), preds = map(model, augment)) %>% 
  unnest(preds)


tdf <- cbind(tdf, smooth_mvenus1$.fitted,smooth_mvenus2$.fitted,smooth_mvenus3$.fitted)
colnames(tdf)[46] = 'mVenus_smooth1'
colnames(tdf)[47] = 'mVenus_smooth2'
colnames(tdf)[48] = 'mVenus_smooth3'

mspanbig <- 0.32
mspanmed <- 0.19
mspansmall <- 0.16

smooth_mcherry1 <- mcherry_to_nest %>% 
  nest(-TRACK_ID) %>% 
  mutate(model = map(data, my_loess_mcherry, span = mspanbig), preds = map(model, augment)) %>% 
  unnest(preds)

smooth_mcherry2 <- mcherry_to_nest %>% 
  nest(-TRACK_ID) %>% 
  mutate(model = map(data, my_loess_mcherry, span = mspanmed), preds = map(model, augment)) %>% 
  unnest(preds)

smooth_mcherry3 <- mcherry_to_nest %>% 
  nest(-TRACK_ID) %>% 
  mutate(model = map(data, my_loess_mcherry, span = mspansmall), preds = map(model, augment)) %>% 
  unnest(preds)


tdf <- cbind(tdf, smooth_mcherry1$.fitted,smooth_mcherry2$.fitted,smooth_mcherry3$.fitted)
colnames(tdf)[49] = 'mCherry_smooth1'
colnames(tdf)[50] = 'mCherry_smooth2'
colnames(tdf)[51] = 'mCherry_smooth3'



# For creating z scores, first find the mean and standard deviation of all mcherry and mvenus measurements across conditions.
mean_mCherry <- mean(mcherry_all)
sd_mCherry <- sd(mcherry_all)
mean_mVenus <- mean(mvenus_all)
sd_mVenus <- sd(mvenus_all)

# Compute the z scores for each mCherry and mVenus value, and add a new column for those z score values
tdf <- tdf %>% mutate (mcherry_zscore = (mCherry-mean_mCherry)/sd_mCherry,
                       mvenus_zscore = (mVenus-mean_mVenus)/sd_mVenus)

# The last track is a manually-created background track. The track ID will therefore be the largest.
# So, we want to take the highest Track ID, and put the largest values for mCherry and mVenus under
# mVenusMax and mCherryMax monikers, with specifiers for raw score, Z score, and Fluorscore.

MaxTrackID <- max(tdf$TRACK_ID)
BackgroundTrack <- subset(tdf, TRACK_ID == MaxTrackID)

#HBRc = High Background Correction Factor for mCherry. Default should be 1; if the lense is 5X, you may need to mess with this
HBRc <- 1.2
#HBRv = High Background Correction Factor for mVenus. Default should be 1; if the lense is 5X, you may need to mess with this
HBRv <- 1.1

MaxmCherryRaw <- max(BackgroundTrack$mCherry)*HBRc
MaxmVenusRaw <- max(BackgroundTrack$mVenus)*HBRv

MaxmCherryZ <- max(BackgroundTrack$mcherry_zscore)*HBRc
MaxmVenusZ <- max(BackgroundTrack$mvenus_zscore)*HBRv

MaxmCherryFluor <- max(BackgroundTrack$mcherry_fluorscore)*HBRc
MaxmVenusFluor <- max(BackgroundTrack$mvenus_fluorscore)*HBRv

# create a subset frame that has TRACK_ID, FRAME, the mcherry z scores and the mvenus scores
zscoreplotting <- tdf %>%
  select(TRACK_ID, FRAME, mcherry_zscore, mvenus_zscore) %>%
  gather(key = "Fluoresence", value = "ZScore", -FRAME, -TRACK_ID)

# create a subset frame that has TRACK_ID, FRAME, the mcherry low scores and the mvenus low scores
fluorscoreplotting <- tdf %>%
  select(TRACK_ID, FRAME, mcherry_fluorscore, mvenus_fluorscore) %>%
  gather(key = "Fluoresence", value = "FluoresenceScore", -FRAME, -TRACK_ID)

# create a subset frame that has TRACK_ID, FRAME, and the smoothed fluorscores

smoothfluorscoreplotting1 <- tdf %>%
  select(TRACK_ID, FRAME, mCherry_smooth1, mVenus_smooth1) %>%
  gather(key = "Fluoresence", value = "Smoothed_FluoresenceScore", -FRAME, -TRACK_ID)

smoothfluorscoreplotting2 <- tdf %>%
  select(TRACK_ID, FRAME, mCherry_smooth2, mVenus_smooth2) %>%
  gather(key = "Fluoresence", value = "Smoothed_FluoresenceScore", -FRAME, -TRACK_ID)

smoothfluorscoreplotting3 <- tdf %>%
  select(TRACK_ID, FRAME, mCherry_smooth3, mVenus_smooth3) %>%
  gather(key = "Fluoresence", value = "Smoothed_FluoresenceScore", -FRAME, -TRACK_ID)


# Now let's get graphing! Woo!


# for each track number
for (i in tracks){
  # open jpeg to save the ensuing plot
  filename = paste0(condition, " Track ", i, ".jpeg")
  jpeg(filename, width = 3200, height = 1700)
  
  # subset tdf to only include track i 
  rawsubplot <- subset(raw, TRACK_ID == i, select = c("FRAME", "Fluoresence", "Fluoresence_Intensity"))
  title = paste0(condition, " Track ", i)
  
  # make raw plot for track i
  rawplot <- ggplot(rawsubplot, aes(x = FRAME, y = Fluoresence_Intensity)) + 
    geom_line(aes(color = Fluoresence, linetype = Fluoresence), size = 1.8) + 
    scale_color_manual(values = c("#D7431D", "#3CCD9A")) +
    geom_hline(yintercept = MaxmCherryRaw, linetype = "longdash", colour = "#D7431D") + 
    geom_hline(yintercept = MaxmVenusRaw, linetype = "longdash", colour = "#3CCD9A") + 
    theme(legend.text = element_text(size= 20)) +
    theme(legend.title = element_text(size = 25)) +
    theme(axis.title = element_text(size = 20)) +
    ##next two lines are for black background and no panel grid marks
    theme_black(base_size = 24) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 25)) + ylim(-10,100) + xlim(0, 190)
  
  # subset Zscoreplotting to only include track i
  zscoresubplot <- subset(zscoreplotting, TRACK_ID == i, select = c("FRAME", "Fluoresence", "ZScore"))
  title = paste0(condition, " Track ", i, " Zscore ")
  
  zplot <- ggplot(zscoresubplot, aes(x = FRAME, y = ZScore)) + 
    geom_line(aes(color = Fluoresence, linetype = Fluoresence), size = 1.8) + 
    theme_black() +
    scale_color_manual(values = c("#D7431D", "#3CCD9A")) +
    geom_hline(yintercept = MaxmCherryZ, linetype = "longdash", colour = "#D7431D") + 
    geom_hline(yintercept = MaxmVenusZ, linetype = "longdash", colour = "#3CCD9A") + 
    theme(legend.text = element_text(size= 20)) +
    theme(legend.title = element_text(size = 25)) +
    theme(axis.title = element_text(size = 20)) +
    ##next two lines are for black background and no panel grid marks
    theme_black(base_size = 24) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 25)) + ylim(-10, 100) + xlim(0, 190)
  
  # subset fluorscoreplotting to only include track i
  fluorscoresubplot <- subset(fluorscoreplotting,TRACK_ID == i, select = c("FRAME", "Fluoresence", "FluoresenceScore"))
  title = paste0(condition, " Track ", i, " FluorScore ")
  
  scoreplot <- ggplot(fluorscoresubplot, aes(x = FRAME, y = FluoresenceScore)) + 
    geom_line(aes(color = Fluoresence, linetype = Fluoresence), size = 1.8) + 
    scale_color_manual(values = c("#D7431D", "#3CCD9A")) +
    geom_hline(yintercept = MaxmCherryFluor, linetype = "longdash", colour = "#D7431D") + 
    geom_hline(yintercept = MaxmVenusFluor, linetype = "longdash", colour = "#3CCD9A") + 
    theme(legend.text = element_text(size= 20)) +
    theme(legend.title = element_text(size = 25)) +
    theme(axis.title = element_text(size = 25)) +
    ##next two lines are for black background and no panel grid marks
    theme_black(base_size = 24) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 25)) + ylim(-10, 100) + xlim(0, 190)
  
  #finalplot <- ggarrange(rawplot, zplot, scoreplot,
  #                      labels = c("A", "B", "C"),
  #                    ncol = 3, nrow = 1)
  
  # subset fluorscoreplottingsmooth to only include track i
  
  smooth_fluorscoresubplot1 <- subset(smoothfluorscoreplotting1,TRACK_ID == i, select = c("FRAME", "Fluoresence", "Smoothed_FluoresenceScore"))
  title = paste0(condition, " Track ", i, " Large Smoothed FluorScore (1)")
  
  min_x_value <- min(smooth_fluorscoresubplot1$FRAME)
  max_x_value <- max(smooth_fluorscoresubplot1$FRAME)
  min_y_value <- min(smooth_fluorscoresubplot1$Smoothed_FluoresenceScore)
  max_y_value <- max(smooth_fluorscoresubplot1$Smoothed_FluoresenceScore)
  
  smoothplot1 <- ggplot(smooth_fluorscoresubplot1, aes(x = FRAME, y = Smoothed_FluoresenceScore)) + 
    geom_line(aes(color = Fluoresence, linetype = Fluoresence), size = 1.8) + 
    scale_color_manual(values = c("#D7431D", "#3CCD9A")) +
    geom_hline(yintercept = MaxmCherryFluor, linetype = "longdash", colour = "#D7431D") + 
    geom_hline(yintercept = MaxmVenusFluor, linetype = "longdash", colour = "#3CCD9A") + 
    theme(legend.text = element_text(size= 20)) +
    theme(legend.title = element_text(size = 25)) +
    theme(axis.title = element_text(size = 20)) +
    theme(axis.ticks.length.x = unit(0.5, "cm")) +
    theme(axis.text.x = element_text(size = 0.1, angle = 45)) +
    ##next two lines are for black background and no panel grid marks
    theme_black(base_size = 24) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 25)) + ylim(-10, 100) + xlim(0, 190) +
    scale_x_continuous(breaks = seq(min_x_value, max_x_value, by = 10), minor_breaks = 1) +
    scale_y_continuous(breaks = seq(min_y_value, max_y_value, by = 2), minor_breaks = 1)
  
  smooth_fluorscoresubplot2 <- subset(smoothfluorscoreplotting2,TRACK_ID == i, select = c("FRAME", "Fluoresence", "Smoothed_FluoresenceScore"))
  title = paste0(condition, " Track ", i, " Medium Smoothed FluorScore (2)")
  
  min_x_value <- min(smooth_fluorscoresubplot2$FRAME)
  max_x_value <- max(smooth_fluorscoresubplot2$FRAME)
  min_y_value <- min(smooth_fluorscoresubplot2$Smoothed_FluoresenceScore)
  max_y_value <- max(smooth_fluorscoresubplot2$Smoothed_FluoresenceScore)
  
  smoothplot2 <- ggplot(smooth_fluorscoresubplot2, aes(x = FRAME, y = Smoothed_FluoresenceScore)) + 
    geom_line(aes(color = Fluoresence, linetype = Fluoresence), size = 1.8) + 
    scale_color_manual(values = c("#D7431D", "#3CCD9A")) +
    geom_hline(yintercept = MaxmCherryFluor, linetype = "longdash", colour = "#D7431D") + 
    geom_hline(yintercept = MaxmVenusFluor, linetype = "longdash", colour = "#3CCD9A") + 
    theme(legend.text = element_text(size= 20)) +
    theme(legend.title = element_text(size = 25)) +
    theme(axis.title = element_text(size = 20)) +
    theme(axis.ticks.length.x = unit(0.25, "cm")) +
    theme(axis.text.x = element_text(size = 0.1, angle = 45)) +
    ##next two lines are for black background and no panel grid marks
    theme_black(base_size = 24) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 25)) + ylim(-10, 100) + xlim(0, 190) +
    scale_x_continuous(breaks = seq(min_x_value, max_x_value, by = 10), minor_breaks = 1) +
    scale_y_continuous(breaks = seq(min_y_value, max_y_value, by = 2), minor_breaks = 1)
  
  smooth_fluorscoresubplot3 <- subset(smoothfluorscoreplotting3,TRACK_ID == i, select = c("FRAME", "Fluoresence", "Smoothed_FluoresenceScore"))
  title = paste0(condition, " Track ", i, " Small Smoothed FluorScore (3)")
  
  min_x_value <- min(smooth_fluorscoresubplot3$FRAME)
  max_x_value <- max(smooth_fluorscoresubplot3$FRAME)
  min_y_value <- min(smooth_fluorscoresubplot3$Smoothed_FluoresenceScore)
  max_y_value <- max(smooth_fluorscoresubplot3$Smoothed_FluoresenceScore)
  
  smoothplot3 <- ggplot(smooth_fluorscoresubplot3, aes(x = FRAME, y = Smoothed_FluoresenceScore)) + 
    geom_line(aes(color = Fluoresence, linetype = Fluoresence), size = 1.8) + 
    scale_color_manual(values = c("#D7431D", "#3CCD9A")) +
    geom_hline(yintercept = MaxmCherryFluor, linetype = "longdash", colour = "#D7431D") + 
    geom_hline(yintercept = MaxmVenusFluor, linetype = "longdash", colour = "#3CCD9A") +
    theme(legend.text = element_text(size= 20)) +
    theme(legend.title = element_text(size = 25)) +
    theme(axis.title = element_text(size = 20)) +
    theme(axis.text.x = element_text(size = 0.1, angle = 45)) +
    ##next two lines are for black background and no panel grid marks
    theme_black(base_size = 24) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 25)) + ylim(-10, 100) + xlim(0, 190) +
    scale_x_continuous(breaks = seq(min_x_value, max_x_value, by = 10), minor_breaks = 1) +
    scale_y_continuous(breaks = seq(min_y_value, max_y_value, by = 2), minor_breaks = 1)
  
  
  finalplot <- ggarrange(rawplot, zplot, scoreplot, smoothplot1, smoothplot2, smoothplot3,
                         labels = c("A", "B", "C", "D", "E", "F"),
                         ncol = 3, nrow = 2)
  
  
  print(finalplot)
  dev.off()
}

# open pdf to save the ensuing plot
# This one isn't quite as useful, but
# it will plot a all fluoresence values
# normalized to the max value in each channel
# together on a single graph. Useful for
# getting a quick glimpse of "looks like there's
# activity" or "hm. no activity."

smooth_fluorscoresubplotALL <- subset(smoothfluorscoreplotting2, TRACK_ID != MaxTrackID, select = c("TRACK_ID", "FRAME", "Fluoresence", "Smoothed_FluoresenceScore"))


my.max <- function(x) ifelse( !all(is.nan(x)), max(x, na.rm=T), NA)

max_mCherry <- my.max(smooth_fluorscoresubplotALL$Smoothed_FluoresenceScore[smooth_fluorscoresubplotALL$Fluoresence == "mCherry_smooth2"])
max_mVenus <- my.max(smooth_fluorscoresubplotALL$Smoothed_FluoresenceScore[smooth_fluorscoresubplotALL$Fluoresence == "mVenus_smooth2"])

# Divide mCherry values by maximum mCherry value
smooth_fluorscoresubplotALL$mCherry_scaled <- smooth_fluorscoresubplotALL$Smoothed_FluoresenceScore / max_mCherry

# Divide mVenus values by maximum mVenus value
smooth_fluorscoresubplotALL$mVenus_scaled <- smooth_fluorscoresubplotALL$Smoothed_FluoresenceScore / max_mVenus

# Calculate min and max values
min_x_value <- min(smooth_fluorscoresubplotALL$FRAME)
max_x_value <- max(smooth_fluorscoresubplotALL$FRAME)
NormalizedBackgroundmCherry <- (HBRc*max(BackgroundTrack$mcherry_fluorscore))/max_mCherry
NormalizedBackgroundmVenus <- (HBRv*max(BackgroundTrack$mvenus_fluorscore))/max_mVenus

# For background subtraction
smooth_fluorscoresubplotALL$mCherry_scaled <- smooth_fluorscoresubplotALL$mCherry_scaled - NormalizedBackgroundmCherry
smooth_fluorscoresubplotALL$mVenus_scaled <- smooth_fluorscoresubplotALL$mVenus_scaled - NormalizedBackgroundmVenus


filename <- paste0(condition, " ALL TRACKS", ".pdf")
pdf(filename, width = 50, height = 40)


# Create the ggplot object
smoothplotALL <- ggplot(smooth_fluorscoresubplotALL, aes(x = FRAME, y = Smoothed_FluoresenceScore, color = Fluoresence, linetype = Fluoresence, group = interaction(TRACK_ID, Fluoresence))) + 
  geom_line(aes(y = ifelse(Fluoresence == "mCherry_smooth2", mCherry_scaled, mVenus_scaled)), size = 0.9) +
  scale_color_manual(values = c("#D7431D", "#2AA278")) +
  geom_hline(yintercept = NormalizedBackgroundmCherry, linetype = "longdash", colour = "#B43718") + 
  geom_hline(yintercept = NormalizedBackgroundmVenus, linetype = "longdash", colour = "#1E7166") +
  scale_linetype_manual(values = c("mCherry_smooth2" = "solid", "mVenus_smooth2" = "solid")) +
  theme_classic() +
  theme(legend.text = element_text(size= 30),
        legend.title = element_text(size = 35),
        axis.title = element_text(size = 50),
        axis.ticks.length.x = unit(0.25, "cm"),
        axis.text.x = element_text(size = 40, angle = 45, vjust = -0.0000001),
        axis.text.y = element_text(size = 50),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 1, size = 25)) +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(min_x_value, max_x_value, by = 5), minor_breaks = 1) +
  scale_y_continuous(
    name = "normalized Fluorescence Score"
  )

# Print the plot
print(smoothplotALL)

# Close the PDF device
dev.off()

filename <- paste0(condition, HBRc, HBRv, mspanmed, vspanmed, " ALL TRACKS Background Subtraction", ".pdf")
pdf(filename, width = 50, height = 40)

# Create the ggplot object
smoothplotALL <- ggplot(smooth_fluorscoresubplotALL, aes(x = FRAME, y = Smoothed_FluoresenceScore, color = Fluoresence, linetype = Fluoresence, group = interaction(TRACK_ID, Fluoresence))) + 
  geom_line(aes(y = ifelse(Fluoresence == "mCherry_smooth2", mCherry_scaled, mVenus_scaled)), size = .9) +
  scale_color_manual(values = c("#D7431D", "#2AA278")) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "black", size = 2) + 
  scale_linetype_manual(values = c("mCherry_smooth2" = "solid", "mVenus_smooth2" = "solid")) +
  theme_classic() +
  theme(legend.text = element_text(size = 30),
        legend.title = element_text(size = 35),
        axis.title = element_text(size = 50),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.text.x = element_text(size = 40, angle = 45, vjust = -0.001),
        axis.text.y = element_text(size = 50),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 1, size = 25)) +
  ggtitle(title) +
  scale_x_continuous(breaks = seq(min_x_value, max_x_value, by = 5), minor_breaks = 1) +
  scale_y_continuous(
    name = "normalized Fluorescence Score"
  )

# Print the plot
print(smoothplotALL)

# Close the PDF device
dev.off()

#Plots should pop up in the space you previously designated up at the top!

#############################################################################################

## Below is the first of the cell cycle definer functions, with a whole series of
## mCherry / mVenus smoothed relationships and Trackmate parameters considered

smooth_define_cell_cycle <- function(w,x,y,z,a,b,f) {
  #w = mCherry_smooth
  #x = mVenus_smooth
  #y = MaxmCherryFluor
  #z = MaxmVenusFluor
  #a = SNR_Ch4
  #b = SN4_Ch1
  #f = FRAME
  case_when(
    # if mVenus has an extremely steep positive slope and mCherry slope is nowhere near that steep, and mVenus is just below background or greater and mCherry is lowish: S
    w < (y + 1) & x > (z - 0.25) & (abs(lag(x)/x))/(abs(lag(w)/w)) > 0.5 & (abs(x/lead(x)))/(abs(w/lead(w))) > 0.5 ~ "S",
    # mCherry is above background 
    # but mCherry is decreasing precipitously (for at least 3 frames) and mVenus is stable or increasing (for atleast 2 frames): S
    w > y & w < (y + 4) & abs(lag(lag(w))/lag(w)) > 1.5 & abs(lag(w)/w) > 1.5 & abs(w/lead(w)) > 1.5 & abs(lag(x)/x) > 0.75 & abs(x/lead(x)) > 0.75 ~ "S",
    # mCherry value is above Background, and mVenus is below Background: G1
    w > y & x < z ~ "G1",
    # mCherry value is above Background, and mVenus is within 1 of background: G1
    w > y & x < (z + 1) ~ "G1",
    # mCherry value is above Background, and mVenus is above background BUT 
    # mVenus is very low (within 3 of Background) and 
    # mCherry is high relative to mVenus (mCherry/mVenus is more than 2.5): G1
    w > y & x > z & x < (z + 3) & (w/x) > 2.5 ~ "G1",
    # mCherry and mVenus are above background BUT mVenus is only a little above Background and mCherry is very high above background
    w > y & x > z & x < (z + 3) & w > (y + 7) ~ "G1",
    # mCherry is less than background and mVenus is above background: S
    w < y & x > z ~ "S",
    # mCherry and mVenus are both below background, i.e. a blank cell indicating Cdt1 degradation, and mVenus is increasing for at least 2 frames
    w < y & x < z & abs(lag(x)/x) < 1 & abs(x/lead(x)) < 1 & abs(lag(y)/y) > 0.75 & abs(y/lead(y)) > 0.75 ~"S",
    # if mCherry has a extremely steep positive slope and mVenus slope is nowhere near that steep: G1
    (abs(lag(x)/x))/(abs(lag(w)/w)) < 0.5 & (abs(x/lead(x)))/(abs(w/lead(w))) < 0.5 ~ "G1",
    # mCherry is above background but is decreasing (for at least 2 frames) and mVenus is hovering around Background: S
    w > y & (abs(lag(w)/w) > 1 | abs(lag(lag(w)/w) >1)) & abs(w/lead(w)) > 1 & (x > (z - 1) & x < (z + 1)) ~ "S",
    # mCherry and mVenus are both above background 
    # but mCherry is decreasing (for at least 2 frames) and mVenus is increasing (for atleast 2 frames): S
    w > y & x > z & abs(lag(w)/w) > 1 & abs(w/lead(w)) > 1 & abs(lag(x)/x) < 1 & abs(x/lead(x)) < 1 ~ "S",
    # mCherry is barely above background (within 2) and mVenus is above Background: S
    w > y & w < (y+2) & x > z  ~ "S",
    # mCherry is comfortably above background (more than 2) and mVenus is above background: G2
    w > (y + 2) & x > z  ~ "G2",
    # mCherry and mVenus are both above background but the Signal-to-noise ratio is above 1.2: G1
    w > y & x > z & (a > 1.2 | b < 1 | b > 1.11) ~ "G1",
  )
}


# Fill in those anomalies! I.e. when previous and proceeding frames are the same
# cell cycle stage, but the one in the middle isn't (usually due to debris or
# z plane changes)

fill_in_cell_cycle_gaps <- function(s,w,x,y,z,a,b,f) {
  #s = Cell_Cycle_Smooth values or Cell_Cycle_Smooth_Final values
  #w = mCherry_smooth
  #x = mVenus_smooth
  #y = MaxmCherryFluor
  #z = MaxmVenusFluor
  #a = SNR_Ch4
  #b = SN4_Ch1
  #f = FRAME
  case_when(
    is.na(s) & lag(s) == lead(s) ~ lag(s),
    is.na(s) & is.na(lead(s)) & lag(s) == "G1" ~ "G1",
    is.na(s) & is.na(lead(s)) & lag(s) == "S" ~ "S",
    is.na(s) & is.na(lead(s)) & lag(s) == "G2" ~ "G2",
    lag(s) == "G1" & lead(s) == "G1" & s != "G1" ~ "G1",
    lag(s) == "S" & lead(s) == "S" & s != "S" ~ "S",
    lag(s) == "G2" & lead(s) == "G2" & s != "G2" ~ "G2",
    lag(s) == "G1" & lead(s) == "G1" & is.na(s) ~ "G1",
    s == "G1" ~ "G1",
    s == "S" ~ "S",
    s == "G2" ~ "G2",
    is.na(s) & lag(s) == "S" & lead(s) == "G1" ~ "G2",
    # undefined and previous two frames of mVenus are increasing but next two frames of mVenus are decreasing,
    # mCherry increasing the whole time
    is.na(s) & abs(lag(x)/x) < 1 & abs(lag(lag(x)/x)) < 1 & abs(x/lead(x)) > 1 & abs(lead((x/lead(x)))) & abs(lag(w)/w) < 1 & abs(lag((lag(w)/w))) & abs(w/lead(w)) < 1 & abs(lead((w/lead(w)))) ~ "G2",
    # undefined and mCherry is more than three times more than mVenus
    is.na(s) & (w/x) > 3 ~ "G1",
    # undefined and mCherry and mVenus are below background but both increasing
    is.na(s) & w < y & x < z & abs(lag(w)/w) < 1 & abs(w/lead(w)) < 1 & abs(lag(x)/x) < 1 & abs(x/lead(x)) < 1 ~ "G2",
    # undefined and mCherry and mVenus are below background but mcherry is decreasing and mVenus is increasing
    is.na(s) & w < y & x < z & abs(lag(w)/w) > 1 & abs(w/lead(w)) > 1 & abs(lag(x)/x) < 1 & abs(x/lead(x)) < 1 ~ "S",
    # undefined and mCherry and mVenus are below background but mVenus is brighter than mCherry
    is.na(s) & x > w ~ "S",
    # undefined and mCherry is more than one std dev below background, indicating cdt1 degradation
    is.na(s) & w < (y - 1) ~ "S",
    # undefined and mCherry is bigger than mVenus and mCherry is decreasing and mVenus stays below background for 5 frames
    is.na(s) & x < w & x < z & lead(x) < z & lead(lead(x)) < z & lead(lead(lead(x))) < z & lag(x) < z ~ "G1",
    # all cases from first Cell Cycle Calling Pass, except these ones
    # were not able to be defined, so taking the non-defined ones and passing them
    # through the same parameters, but dividing the definition of the background by 2)
    # mCherry value is above Background, and mVenus is below Background: G1
    is.na(s) & w > (y/2) & x < (z/3) ~ "G1",
    # mCherry value is above Background, and mVenus is within 1 std. dev. of background: G1
    is.na(s) & w > (y/2) & x < ((z/3) + 1) ~ "G1",
    # mCherry value is above Background, and mVenus is above background BUT 
    # mVenus is very low (within 3 std devs of Background) and 
    # mCherry is high relative to mVenus (mCherry/mVenus is more than 3): G1
    is.na(s) & w > (y/2) & x > (z/3) & x < ((z/3) + 6) & ((w/2)/x) > 3 ~ "G1",
    # mCherry is less than background and mVenus is above background: S
    is.na(s) & w < (y/2) & x > (z/3) ~ "S",
    # mCherry is above background but is decreasing (for at least 2 frames) and mVenus is hovering # around Background: S
    is.na(s) & w > (y/2) & (abs(lag(w)/w) > 1 | abs(lag(lag(w)/w) >1)) & abs(w/lead(w)) > 1 & (x > ((z/3) - 1) & x < ((z/3) + 1)) ~ "S",
    # mCherry and mVenus are both above background 
    # but mCherry is decreasing (for at least 2 frames) and mVenus is increasing (for atleast 2 frames): S
    is.na(s) & w > (y/2) & x > (z/3) & abs(lag(w)/w) > 1 & abs(w/lead(w)) > 1 & abs(lag(x)/x) < 1 & abs(x/lead(x)) < 1 ~ "S",
    # mCherry is barely above background (within 2 std devs) and mVenus is above Background: S
    is.na(s) & w > (y/2) & w < ((y/2)+2) & x > z  ~ "S",
    # mCherry is comfortably above background (more than 2 std devs) and mVenus is above background: G2
    is.na(s) & w > ((y/2) + 2) & x > (z/3)  ~ "G2",
    # mCherry and mVenus are both above background but the Signal-to-noise ratio is above 1.2: G1
    is.na(s) & w > (y/2) & x > (z/3) & (a > 1.2 | b < 1 | b > 1.11) ~ "G1",
    # if the value is still undefined after all this, then make it G1
    is.na(s) ~ "G1")
}

# Make some empty columns that'll get filled in with cell cycle score calls

tdf_smooth_predict <- tdf %>% mutate(CellCycle_Smooth_1 = "PlaceHold")
tdf_smooth_predict <- tdf %>% mutate(CellCycle_Smooth_Final_1 = "PlaceHold")
tdf_smooth_predict <- tdf %>% mutate(CellCycle_Smooth_2 = "PlaceHold")
tdf_smooth_predict <- tdf %>% mutate(CellCycle_Smooth_Final_2 = "PlaceHold")
tdf_smooth_predict <- tdf %>% mutate(CellCycle_Smooth_3 = "PlaceHold")
tdf_smooth_predict <- tdf %>% mutate(CellCycle_Smooth_Final_3 = "PlaceHold")
tdf_smooth_predict <- tdf %>% mutate(CellCycle_Smooth_4 = "PlaceHold")
tdf_smooth_predict <- tdf %>% mutate(CellCycle_Smooth_Final_4 = "PlaceHold")

tdf_smooth_predict <- tdf_smooth_predict %>% 
  mutate(CellCycle_Smooth_1 = smooth_define_cell_cycle(mcherry_fluorscore, mvenus_fluorscore, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))

tdf_smooth_predict <- tdf_smooth_predict %>% mutate(CellCycle_Smooth_Final_1 = fill_in_cell_cycle_gaps(CellCycle_Smooth_1,mcherry_fluorscore, mvenus_fluorscore, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))
tdf_smooth_predict <- tdf_smooth_predict %>% mutate(CellCycle_Smooth_Final_1 = fill_in_cell_cycle_gaps(CellCycle_Smooth_Final_1,mcherry_fluorscore, mvenus_fluorscore, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))


tdf_smooth_predict <- tdf_smooth_predict %>% 
  mutate(CellCycle_Smooth_2 = smooth_define_cell_cycle(mCherry_smooth1, mVenus_smooth1, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))

tdf_smooth_predict <- tdf_smooth_predict %>% mutate(CellCycle_Smooth_Final_2 = fill_in_cell_cycle_gaps(CellCycle_Smooth_2,mCherry_smooth1, mVenus_smooth1, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))
tdf_smooth_predict <- tdf_smooth_predict %>% mutate(CellCycle_Smooth_Final_2 = fill_in_cell_cycle_gaps(CellCycle_Smooth_Final_2,mCherry_smooth1, mVenus_smooth1, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))

tdf_smooth_predict <- tdf_smooth_predict %>% 
  mutate(CellCycle_Smooth_3 = smooth_define_cell_cycle(mCherry_smooth2, mVenus_smooth2, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))

tdf_smooth_predict <- tdf_smooth_predict %>% mutate(CellCycle_Smooth_Final_3 = fill_in_cell_cycle_gaps(CellCycle_Smooth_3,mCherry_smooth2, mVenus_smooth2, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))
tdf_smooth_predict <- tdf_smooth_predict %>% mutate(CellCycle_Smooth_Final_3 = fill_in_cell_cycle_gaps(CellCycle_Smooth_Final_3,mCherry_smooth2, mVenus_smooth2, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))

tdf_smooth_predict <- tdf_smooth_predict %>% 
  mutate(CellCycle_Smooth_4 = smooth_define_cell_cycle(mCherry_smooth3, mVenus_smooth3, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))

tdf_smooth_predict <- tdf_smooth_predict %>% mutate(CellCycle_Smooth_Final_4 = fill_in_cell_cycle_gaps(CellCycle_Smooth_4,mCherry_smooth3, mVenus_smooth3, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))
tdf_smooth_predict <- tdf_smooth_predict %>% mutate(CellCycle_Smooth_Final_4 = fill_in_cell_cycle_gaps(CellCycle_Smooth_Final_4,mCherry_smooth3, mVenus_smooth3, MaxmCherryFluor, MaxmVenusFluor, SNR_CH4, SNR_CH1,FRAME))

tdf_smooth_predict <- tdf_smooth_predict %>%
  mutate(CellCycleScore_1 = case_when(
    CellCycle_Smooth_Final_1== "G1" ~ 0,
    CellCycle_Smooth_Final_1 == "S" ~ 1,
    CellCycle_Smooth_Final_1 == "G2" ~ 2
  ))

tdf_smooth_predict <- tdf_smooth_predict %>%
  mutate(CellCycleScore_2 = case_when(
    CellCycle_Smooth_Final_2== "NA" ~ 0,
    CellCycle_Smooth_Final_2== "G1" ~ 0,
    CellCycle_Smooth_Final_2 == "S" ~ 1,
    CellCycle_Smooth_Final_2 == "G2" ~ 2
  ))

tdf_smooth_predict <- tdf_smooth_predict %>%
  mutate(CellCycleScore_3 = case_when(
    CellCycle_Smooth_Final_3== "G1" ~ 0,
    CellCycle_Smooth_Final_3 == "S" ~ 1,
    CellCycle_Smooth_Final_3 == "G2" ~ 2
  ))


tdf_smooth_predict <- tdf_smooth_predict %>%
  mutate(CellCycleScore_4 = case_when(
    CellCycle_Smooth_Final_4== "G1" ~ 0,
    CellCycle_Smooth_Final_4 == "S" ~ 1,
    CellCycle_Smooth_Final_4 == "G2" ~ 2
  ))


#create a small table that has all of the unique TRACK_ID's
uniqueID <- sqldf("select distinct (TRACK_ID) from tdf_smooth_predict") 
# Add a sequential index column to the uniqueID data, so that the first TRACK_ID has a corresponding index value in another column
uniqueID <- tibble::rowid_to_column(uniqueID, "Base_ID")
# multiply all these by 2; this will let them be somewhat spaced out in a graph when graphing by the newly-formed GraphID
uniqueID <- uniqueID %>% mutate(GraphID = (Base_ID) * 2)
# add another column that is the sum of all the cell cycle scores within a given Track_ID
SumOfCellCycleScore <-  aggregate(tdf_smooth_predict$CellCycleScore_2, list(tdf_smooth_predict$TRACK_ID),FUN=sum)
# Merge Cell Cycle Score Sum with uniqueID
names(SumOfCellCycleScore)[names(SumOfCellCycleScore) == "Group.1"] <- "TRACK_ID"
names(SumOfCellCycleScore)[names(SumOfCellCycleScore) == "x"] <- "CellCycleSum"
uniqueID <- merge(uniqueID, SumOfCellCycleScore, by = "TRACK_ID")
# create a new table in which tdfx merges with UniqueID ; join them together based on TRACK_ID. # This will make it so
# every TRACK_ID value has its corresponding GraphID value (and BaseID value) in the adjacent  #column
tdfxz <- merge(uniqueID, tdf_smooth_predict, by = "TRACK_ID")

# If there's lots of background in the beginning of your movie,
# use this code to cut off the first few frames of the analysis
# tdfxz <- tdfxz[tdfxz$FRAME > 10,]

# Now, the tracks will get ordered by CellCycleScore
# so the most active tracks appear up at the top

library(dplyr)

all_frames <- data.frame(FRAME = seq(min(tdfxz$FRAME), max(tdfxz$FRAME)))

tdfxz_filled <- merge(all_frames, tdfxz, by = "FRAME", all.x = TRUE)

# Identify missing values
missing_values <- filter(tdfxz_filled, is.na(CellCycleScore_2))

# Add a count of FRAMES for each GraphID
tdfxz <- tdfxz %>%
  add_count(GraphID, name = "FrameCount")

# Create a column where each track starts at time point 0, to get rid of scattered track effects
tdfxz <- tdfxz %>%
  group_by(TRACK_ID) %>%
  group_modify(~ mutate(.x, Adjusted_Frame = FRAME - min(FRAME, na.rm = TRUE)))

# Order the data by CellCycleSum, FrameCount, and GraphID
tdfxz <- tdfxz %>%
  arrange(CellCycleSum, FrameCount, GraphID) %>%
  mutate(GraphID = factor(GraphID, levels = unique(GraphID)))

# Convert CellCycleScore_1 to a factor with specified levels
tdfxz$CellCycleScore_1 <- factor(tdfxz$CellCycleScore_1, levels = c("0", "1", "2", NA))

# Convert CellCycleScore_2 to a factor with specified levels
tdfxz$CellCycleScore_2 <- factor(tdfxz$CellCycleScore_2, levels = c("0", "1", "2", NA))

# Convert CellCycleScore_3 to a factor with specified levels
tdfxz$CellCycleScore_3 <- factor(tdfxz$CellCycleScore_3, levels = c("0", "1", "2", NA))

# Convert CellCycleScore_4 to a factor with specified levels
tdfxz$CellCycleScore_4 <- factor(tdfxz$CellCycleScore_4, levels = c("0", "1", "2", NA))


#create a combined bunch of colors corresponding to numbers (recognized by their character, not their value)
cols <- c("0" = "#696363", "1" = "#5DD5AB", "2" = "#5DD5AB", "NA" = "#696363")

#By making objects with a specified output as below, we can make nifty subsets for a prettier #legend
G1G <- expression(paste("G"["1"]))
S1S <- expression(paste("S"))
G2G <- expression(paste("G"["2"]))
#MMM <- expression(paste("M"))
Unknown <- expression(paste("Undefined"))

#get all the unique track IDs together
Track_IDs <- uniqueID %>% pull(TRACK_ID)

#make all of them have 2 digits, such that (for example) 9 becomes 09; if the Track_IDs go into the 100's, change "%02d" to "03d"
Track_IDs <- sprintf("%02d", Track_IDs)

# make a list of all the GraphIds
Graph_IDs <- uniqueID %>% pull(GraphID)

#  Define the number of tracks in the set ; the subtraction  of one represents taking away the background track
N <- (length(Track_IDs) - 1)

# Make a subset of tdfxz that does not have the TrackID with the largest ID value, which is (SHOULD BE!!) the background track

NoBackground_tdfxz <- tdfxz[tdfxz$TRACK_ID != MaxTrackID,]


# You're almost there!
# Make the plot of all the tracks, color coded showing their cell cycle state 
# at any given time point, with their Track_ID number on the left axis
# open pdf to save the ensuing plot

titlecycle <- paste0(condition, " Black Fluorscore-Based Cell Cycle Graph 2 ", "  |  n = ", N)
filename = paste0(condition," HBRv = ", HBRv, " HBRc = ", HBRc, " Black Fluorscore-Based Cell Cycle Graph 2 "," n = ", N, ".pdf")


pdf(filename, width = 10, height = 8)
NoBackground_tdfxz$CellCycleScore_1 <- factor(NoBackground_tdfxz$CellCycleScore_1, levels = c("0", "1", "2", NA))
cellcycleplot1 <-  ggplot(NoBackground_tdfxz, aes(Adjusted_Frame, reorder(GraphID, CellCycleSum), color = CellCycleScore_1, fill = CellCycleScore_1)) +
  geom_tile(color = "white", size = 0.15) +
  scale_color_manual(
    drop = FALSE,
    name = "Cell Cycle Phase",
    values = cols,
    aesthetics = c("color", "fill"),
    labels = c(G1G, S1S, G2G, Unknown)) +
  labs(x = "Tracked Frame", y = "Cell Tracks") +
  theme_classic(base_size = 28) +
  theme(
    # Remove panel border
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "white", size = 0.25),
    panel.background = element_blank(),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(size= 10),
    legend.text.align = 0,
    axis.text.y.left = element_text(size = 1.5),
    axis.line.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank()) +
  scale_y_discrete(breaks = Graph_IDs, labels = Track_IDs) +  # Using a discrete scale for ordered factors
  ggtitle(titlecycle) + theme(plot.title = element_text(hjust = 0.5, size = 5))

print(cellcycleplot1)

dev.off()

# Plotting code
titlecycle <- paste0(condition, " Large Smooth Cell Cycle Graph 2", "  |  n = ", N)
filename <- paste0(condition, " HBRv = ", HBRv, " HBRc = ", HBRc, " mspan = ", mspanbig, " v span = ", vspanbig, " Black Large Smooth Cell Cycle Graph", " n = ", N, ".pdf")

pdf(filename, width = 9, height = 8.5)
NoBackground_tdfxz$CellCycleScore_2 <- factor(NoBackground_tdfxz$CellCycleScore_2, levels = c("0", "1", "2", NA))
cellcycleplot2 <- ggplot(NoBackground_tdfxz, aes(Adjusted_Frame, reorder(GraphID, CellCycleSum), color = CellCycleScore_2, fill = CellCycleScore_2)) +
  geom_tile(color = "white", size = 0.15) +
  scale_color_manual(
    drop = FALSE,
    name = "Cell Cycle Phase",
    values = cols,
    aesthetics = c("color", "fill"),
    labels = c(G1G, S1S, G2G, Unknown)) +
  labs(x = "Tracked Frame", y = "Cell Tracks") +
  theme_classic(base_size = 28) +
  theme(
    # Remove panel border
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "white", size = 0.25),
    panel.background = element_blank(),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(size= 10),
    legend.text.align = 0,
    axis.text.y.left = element_text(size = 1.5),
    axis.line.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank()) +
  scale_y_discrete(breaks = Graph_IDs, labels = Track_IDs) +  # Using a discrete scale for ordered factors
  ggtitle(titlecycle) + theme(plot.title = element_text(hjust = 0.5, size = 5))

print(cellcycleplot2)


dev.off()


titlecycle <- paste0(condition, " Medium Smooth Cell Cycle Graph 2", "  |  n = ", N)
filename <- paste0(condition, " HBRv = ", HBRv, " HBRc = ", HBRc, " mspan = ", mspanmed, " v span = ", vspanmed, " Black Medium Smooth Cell Cycle Graph", " n = ", N, ".pdf")

pdf(filename, width = 10, height = 8)
NoBackground_tdfxz$CellCycleScore_3 <- factor(NoBackground_tdfxz$CellCycleScore_3, levels = c("0", "1", "2", NA))
cellcycleplot3 <- ggplot(NoBackground_tdfxz, aes(Adjusted_Frame, reorder(GraphID, CellCycleSum), color = CellCycleScore_3, fill = CellCycleScore_3)) +
  geom_tile(color = "white", size = 0.15) +
  scale_color_manual(
    drop = FALSE,
    name = "Cell Cycle Phase",
    values = cols,
    aesthetics = c("color", "fill"),
    labels = c(G1G, S1S, G2G, Unknown)) +
  labs(x = "Tracked Frame", y = "Cell Tracks") +
  theme_classic(base_size = 28) +
  theme(
    # Remove panel border
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "white", size = 0.25),
    panel.background = element_blank(),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(size= 10),
    legend.text.align = 0,
    axis.text.y.left = element_text(size = 1.5),
    axis.line.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank()) +
  scale_y_discrete(breaks = Graph_IDs, labels = Track_IDs) +  # Using a discrete scale for ordered factors
  ggtitle(titlecycle) + theme(plot.title = element_text(hjust = 0.5, size = 10))

print(cellcycleplot3)

dev.off()

titlecycle <- paste0(condition, " Small Smooth Cell Cycle Graph 2", "  |  n = ", N)
filename <- paste0(condition, " HBRv = ", HBRv, " HBRc = ", HBRc, " mspan = ", mspansmall, " v span = ", vspansmall, " Black Small Smooth Cell Cycle Graph", " n = ", N, ".pdf")

pdf(filename, width = 10, height = 8)
NoBackground_tdfxz$CellCycleScore_4 <- factor(NoBackground_tdfxz$CellCycleScore_4, levels = c("0", "1", "2", NA))
cellcycleplot4 <- ggplot(NoBackground_tdfxz, aes(Adjusted_Frame, reorder(GraphID, CellCycleSum), color = CellCycleScore_4, fill = CellCycleScore_4)) +
  geom_tile(color = "white", size = 0.15) +
  scale_color_manual(
    drop = FALSE,
    name = "Cell Cycle Phase",
    values = cols,
    aesthetics = c("color", "fill"),
    labels = c(G1G, S1S, G2G, Unknown)) +
  labs(x = "Tracked Frame", y = "Cell Tracks") +
  theme_classic(base_size = 28) +
  theme(
    # Remove panel border
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change axis line
    axis.line = element_line(colour = "white", size = 0.25),
    panel.background = element_blank(),
    legend.title = element_text(color = "black", size = 10),
    legend.text = element_text(size= 10),
    legend.text.align = 0,
    axis.text.y.left = element_text(size = 1.5),
    axis.line.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank()) +
  scale_y_discrete(breaks = Graph_IDs, labels = Track_IDs) +  # Using a discrete scale for ordered factors
  ggtitle(titlecycle) + theme(plot.title = element_text(hjust = 0.5, size = 10))

print(cellcycleplot4)

dev.off()

# I hope this was useful! Thanks for giving it a go!

