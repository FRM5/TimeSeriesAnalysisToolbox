# PRECIPITATION DATA ANALYSIS
# Creating IDF curves from tipping bucket precipitation data
# Created by Gaby Grundemann, EMFRM
# November 2016

### General overview ###
# This script analyses tipping bucket precipitation data
# The raw data is derived from tipping points per 0.1 mm intervals
# The data is then computed into a full timeseries dataframe
# Hereafter, moving averages for different durations are computed
# These values will be aggregated per year
# The result is the maximum precipitation depth per duration 
# Then the maximum intensity for the different years will be computed
# There are two distributions fitted to the data, GEV and Gumbel
# The built-in packages are used to do this
# The one with the lowest error will be used for computing the IDF curves
# The script is set up in such a way that it runs automatically 
# The only settings which are needed in order to run the script are the following:
#     - Durations of the precipitation events, in minutes
#     - Return periods for computing the IDF curves, in years
#     - The working directory
#     - The raw data, which consists of two columns:
#           - Date
#           - Time

### Required libraries ###
library(xts)          # For time-series analysis
library(plyr)         # For handling data, in particular for the function "rename"
library(hydroTSM)     # For computing the moving average
library(lubridate)    # For advance working with dates, in particular for the function "year" 
library(extRemes)     # For performing the extreme value analysis, incl. sensitivity
library(reshape2)     # For reshaping the table as input for the graph
library(ggplot2)      # For displaying the IDF graphs

### Settings, adjust these according to your preference ###
Duration=c(5, 15, 30, 60, 120, 360, 720, 1440)    # duration in minutes 
Return_Period=c(2, 5, 10, 50, 100, 200)           # return periods in years

### Load the data in .txt format ###
# Set the working directory which holds the data file you wish to read
# Read the data. The columns are "date" and "time"
# Convert the data and time stamps to a format recognized by R 
# In order to reduce the size, the interval will be one minute
# For more accuracy, add the seconds as well. But it will take much longer to run
# An extra column with the depth is added
setwd("~/Documents/UPC/Assignment")    
My_Data=read.table("21.txt", header=TRUE)        
My_Data <- data.frame(as.POSIXct(paste(My_Data$date, My_Data$time), format="%m/%d/%Y %H:%M"))   
colnames(My_Data) <- c("DateTime")
My_Data$Depth <- 0.1

### Convert to timeseries format recognized by R ###
# Convert to a XTS format and aggregate the same timesteps 
# Create an empty time series file, from the beginning to the end of the data with 1 minute intervals
# Add the precipitation data to the empty time series
# Convert the empty values to 0
# Convert to a dataframe format
TS <- xts(My_Data[,-1], order.by=My_Data[,1])
TS_Agg <- aggregate(TS, identity, sum)
TS_Empty <- zoo(, seq(start(TS_Agg), end(TS_Agg), by="1 min"))
TS_Merged <- merge(TS_Agg, TS_Empty)
TS_Merged[is.na(TS_Merged)] <- 0
TS_DF <- data.frame(date=index(TS_Merged), coredata(TS_Merged))

### Compute the moving averages ###
# There are two ways to compute the moving averages and annual maximums:
#      The first one is with the "ma" function from the hydroTSM package
#      The second one is with the "rollapply" function 
# Only the first option is used in this script
# Change if the rollapply is prefered over the moving average function

## FIRST OPTION
# Moving average with ma 
for(i in 1:length(Duration)){
  Mov_Avg <- ma(TS_Merged, win.len = Duration[i], FUN = sum)
  Mov_Avg <- data.frame(DateTime=index(Mov_Avg), "TS_Agg" = coredata(Mov_Avg))
  if (i==1){
    Max_Year <- aggregate(Mov_Avg$TS_Agg ~ year(Mov_Avg$DateTime), FUN = max) 
  }else{
    Mov_Avg<- aggregate(Mov_Avg$TS_Agg ~ year(Mov_Avg$DateTime), FUN = max)
    Max_Year<-merge(Max_Year, Mov_Avg, by = "year(Mov_Avg$DateTime)")
  }
  if(Duration[i]<60){
    Colname=paste0(Duration[i], "min")
  }else{
    Colname=paste0(Duration[i]/60, "hr")
  }
  Max_Year<-rename(Max_Year, replace = c("Mov_Avg$TS_Agg" = Colname))
}
Max_Year<-rename(Max_Year, replace = c("year(Mov_Avg$DateTime)" = "Year"))
Max_Year <- Max_Year[-c(1), ] # remove the first year
# Applicable only if there are too little values in that year

## SECOND OPTION
# Turned off, to turn on delte the #
# Moving average with rollapply 
#for (i in 1:length(Duration)) {
#  Mov_Avg <- TS_DF[1]
#  Mov_Avg$Precipitation <- rollapply(TS_DF[2], as.numeric(Duration[i]), sum, fill = NA, align='left')
#  if (i==1){
#    Max_Year <- aggregate(Mov_Avg$Precipitation ~ year(Mov_Avg$date), FUN = max)
#  }else{
#    Mov_Avg <- aggregate(Mov_Avg$Precipitation ~ year(Mov_Avg$date), FUN = max)
#    Max_Year<-merge(Max_Year, Mov_Avg, by="year(Mov_Avg$date)")
#  }
#  if(Duration[i]<60){
#    Colname=paste0(Duration[i], "min")
#  }else{
#    Colname=paste0(Duration[i]/60, "hr")
#  }
#  Max_Year<-rename(Max_Year, replace = c(TS_Agg = Colname))
#}
#Max_Year<-rename(Max_Year, replace = c("year(Mov_Avg$date)" = "Year"))
#Max_Year <- Max_Year[-c(1), ] # remove the year 1994

### Compute the annual maximum intensities ###
# The maximum intensities are in [mm/hr]
Max_Intensity <- Max_Year
for(i in 1:length(Duration)){
  Max_Intensity[,i+1]=Max_Intensity[,i+1]*(60/Duration[i])
}

### Distributions ###
# Apply the Gumbel distribution to the data
Gumbel <- matrix(data=NA, ncol=length(Return_Period))
for(i in 1:length(Duration)){
  value=as.numeric(Max_Intensity[,1+i])
  Fit_Gum <- fevd(value, Max_Intensity, scale=1, type="Gumbel")
  RetLev <- return.level(Fit_Gum, conf=0.05, return.period=(Return_Period))
  if (Duration[i]==1) {
    Gumbel <- t(data.frame(as.numeric(RetLev)))
  }else {
    Gumbel2 <- t(data.frame(as.numeric(RetLev)))
    Gumbel <- rbind(Gumbel, Gumbel2)
  }
}
Gumbel<-Gumbel[-which(apply(Gumbel,1,function(x)all(is.na(x)))),]
row.names(Gumbel) <- Duration

# Apply the GEV distribution to the data
GEV <- matrix(data=NA, ncol=length(Return_Period))
for(i in 1:length(Duration)){
  value=as.numeric(Max_Intensity[,1+i])
  Fit_GEV <- fevd(value, Max_Intensity, scale=1, type="GEV")
  RetLev <- return.level(Fit_GEV, conf=0.05, return.period=(Return_Period))
  if (Duration[i]==1) {
    GEV <- t(data.frame(as.numeric(RetLev)))
  }else {
    GEV2 <- t(data.frame(as.numeric(RetLev)))
    GEV <- rbind(GEV, GEV2)
  }
}
GEV<-GEV[-which(apply(GEV,1,function(x)all(is.na(x)))),]
row.names(GEV) <- Duration

### Sensitivity analysis ###
# Determine if the Gumbel or the GEV is prefered
# Lower AIC and BIC values are better
AICBIC_GEV <- summary(Fit_GEV, silent=TRUE)
AICBIC_Gum <- summary(Fit_Gum, silen=TRUE)
if (AICBIC_Gum$AIC & AICBIC_Gum$BIC < AICBIC_GEV$AIC & AICBIC_GEV$BIC) {
  FINAL <- Gumbel
} else {
  FINAL <- GEV
}

### IDF Curves ###
# Create a graph displaying the IDF curves
# First copy the table with the values of the Gumbel or GEV distribution
# This final table was determined by the previously executed sensitivity analysis
IDF_Table <- data.frame(FINAL) 
Col_Names_RP <- matrix(data=NA) # Create an empty data frame to host the names of the return periods
# That data frame will contain the return periods, needed for the legend
for (i in 1:length(Return_Period)){
  Col_Names_RP[i] <- paste0(Return_Period[i], " years")
}
colnames(IDF_Table) <- c(Col_Names_RP) # Copy the names to the IDF table
Hours <- as.data.frame(c(as.numeric(as.character(t(Duration))))/60) 
# This converts the final durations to hours instead of minutes for the X axis
colnames(Hours) <-"Duration [hrs]" # change column name
IDF_Table <- cbind(Hours, IDF_Table) # add the duration table to the final IDF table.
IDF_Plot<- melt(IDF_Table, id.vars="Duration [hrs]") # order the table by duration, needed for the plotting
IDF_Plot <- ggplot(IDF_Plot, aes(x=`Duration [hrs]`, y=value, color=variable)) + geom_line() + geom_point() + 
  ylab("Intensity [mm/hr]") + labs(title="Intensity Duration Frequency curves") + labs(color="Return Periods:") +
  theme(legend.position = c(.8, .7))
plot(IDF_Plot) 
