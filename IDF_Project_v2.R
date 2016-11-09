#Name: Heather Murdock
#IDF curve generation using R

#STEP 1: load libraries used

library(hydroTSM) #hydroTSM - Management, analysis, interpolation and plot of hydrological time series, with focus on hydrological modelling
library(lubridate) #Lubridate provides tools that make it easier to parse and manipulate dates.
library(xts) #Easily convert one of R's many time-series (and non-time-series) classes to a true time-based object which inherits all of zoo's methods, while allowing for new time-based tools where appropriate.
library(zoo) #zoo - zoo is the creator for an S3 class of indexed totally ordered observations which includes irregular time series.
library(ggplot2) #ggplot2 - Create a new ggplot plot.
library(reshape2) #reshape tables

#STEP 2: Load tipping bucket data for dataset 4 - 08/05/1996 to 08/30/2016
setwd("~/Documents/R-Projects")
mydata <- read.table(file.choose(TRUE), header = TRUE)

#Durations = 30mins, 1hr, 6hr, 12hr, 24hr // Durations in minutes = 30, 60, 360, 720, 1440
#Return periods = 5, 10, 25, 50, 100 years

YearP <- c(1996:2016)
S_durs <- data.frame(30, 60, 360, 720, 1440)
RPs <- data.frame(5, 10, 25, 50, 100)

#label columns
colnames(mydata)[1]= "DateTime"
colnames(mydata)[2]= "Time"

#convert formate of columns and set appropriate timezone (seconds not needed)
mydata$DateTime <- as.POSIXct(paste(mydata$DateTime, mydata$Time), format="%m/%d/%Y %H:%M")

#each tip is 0.1mm - do I need to add a column for Depth??
mydata$Depth <- 0.1
mydata$Time <- NULL #remove time column

#finaldata <-subset(mydata,diff(mydata$Date) >0) #removes overlaping times - clicks in the same second

#STEP 3: Time Series Analysis

#set up data as a time series
TS <- xts(mydata[,-1], order.by= mydata[,1]) #create timeseries for depth ordered by Date and Time?
TS_p <- aggregate(TS, identity, sum) #will not include 1st and last year but functions!
#TS_p <- period.apply(TS, endpoints(TS, "minutes",1), sum) #in a given period - minutes (all data included with period.apply)
TS_em <- zoo( , seq(start(TS_p), end(TS_p), by="1 min")) #one minute uniform timeseries
TS_merged <- merge(TS_p, TS_em)
TS_merged[is.na(TS_merged)] <- 0 #If there is no data then depth is 0 ie: no tipping/no rain
TS.df <- data.frame(date=index(TS_merged), coredata(TS_merged)) #data frame

colnames(TS.df) <- c("DateTime", "Depth")

#STEP 4: Calculate Annual Maximum intensiies

message("
==============================================================
 This code uses durations of 30mins, 1hr, 6hr, 12hr, 24hr (Canadian IDF Standards)
===============================================================")

#max durations in minutes 30, 60, 360, 720, 1440
#rollapply() #moving window
#aggregate() #maximum for the year

dmax30 <- TS.df[1]
dmax30$Depth <- rollapply(TS.df[2], width=30, sum, fill=NA, align='center')
dmax30 <- aggregate((dmax30$Depth) ~ year(dmax30$DateTime), FUN=max)

dmax60 <- TS.df[1]
dmax60$Depth <- rollapply(TS.df[2], width=60, sum, fill=NA, align='center')
dmax60 <- aggregate((dmax60$Depth) ~ year(dmax60$DateTime), FUN=max)

dmax360 <- TS.df[1]
dmax360$Depth <- rollapply(TS.df[2], width=360, sum, fill=NA, align='center')
dmax360 <- aggregate((dmax360$Depth) ~ year(dmax360$DateTime), FUN=max)

dmax720 <- TS.df[1]
dmax720$Depth <- rollapply(TS.df[2], width=720, sum, fill=NA, align='center')
dmax720 <- aggregate((dmax720$Depth) ~ year(dmax720$DateTime), FUN=max)

dmax1440 <- TS.df[1]
dmax1440$Depth <- rollapply(TS.df[2], width=1440, sum, fill=NA, align='center')
dmax1440 <- aggregate((dmax1440$Depth) ~ year(dmax1440$DateTime), FUN=max)

dmaxframe <- data.frame(year = YearP, dmax30=dmax30[,2], dmax60=dmax60[,2], dmax360=dmax360[,2], dmax720=dmax720[,2], dmax1440=dmax1440[,2])

# Possible to import table to complete remaining code - moving average section is the most time consuming
# dmaxframe <- read.table("dmaxframe.txt", header=TRUE)

#STEP 5: Compute GUMBEL for data

message("
        ==============================================================
        Computing GUMBEL
        ===============================================================")
#Factors to compute intensities
I_div <- data.frame(year=1, c1=2, c2=1, c3=1/6, c4=1/12, c5=1/24)

#Create table with correct dimensions to apply factors
Idata <- data.frame(matrix(0,ncol = 6, nrow = length(YearP)))
for (i in 1:length(YearP)) Idata[i,] = I_div

#Apply factors to compute the maximum intensity for each storm duration       
dmaxIframe <- dmaxframe * Idata

#Compute Standard Deviation Table
Sd_p <- data.frame(matrix(0,ncol = length(S_durs), nrow = length(RPs)))
for (i in 1:length(S_durs)) 
  Sd_p[i] <- sd(dmaxIframe[,i])

#Label Standard Deviation Table
colnames(sd_I) <- S_durs
row.names(sd_I) <- RPs

#Compute average intensity table
P_ave <- c(0)
for (i in 1: length(S_durs))
  (P_ave[i] <- mean(dmaxIframe[,i]))

#Compute K values table
K <- data.frame(matrix(0,ncol = 1, nrow = length(RPs)))
for (i in 1: length(RPs))
  K[i] <- (0-(sqrt(6)/pi)*(0.5772+log(log(RPs[i]/((RPs[i])-1)))))

#Label rows and columns of K values table
colnames(K) <- S_durs
row.names(K) <- RPs
  
#STEP 5: Calculating the IDF table

IDF <- P_ave + K*Sd_p
  
#STEP 6: Plot the IDF table

#transpose
IDFt <- t(IDF)
Duration <- t(S_durs)

#use melt function from reshape2 - reshape tables

IDFt2 <- cbind(Duration, IDFt)
IDF_melt <- melt(IDFt, id.vars="Duration")

IDF_Plot <- ggplot(IDF_melt, aes(x=Duration, y=value, color=variable)) + geom_line() + geom_point() + 
  ylab("Intensity [mm/hr]") + labs(title="Intensity Duration Frequency curves") + labs(color="Return Periods:") +
  theme(legend.position = c(.8, .7))

plot(IDF_Plot)
