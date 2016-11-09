#Name: Andre Borin Venturini
#Flash Flood Assigment

#setting working directory
setwd("~/UPC/R_assignment/Andre")
#setting libraries
#eXtensible Time Series: Provide for uniform handling of R's different time-based data classes by extending zoo
library("xts", lib.loc="~/R/win-library/3.3")
#Package to work with h date-times and time-spans
library("lubridate", lib.loc="~/R/win-library/3.3")
#Tools For Highfrequency Data Analysis(irregular TS)
library("highfrequency", lib.loc="~/R/win-library/3.3")
#A Toolbox for Non-Tabular Data Manipulation
library("rlist", lib.loc="~/R/win-library/3.3")
#library for advanced plotting
library("ggplot2", lib.loc="~/R/win-library/3.3")
#library to perform computation on irregular time series
library("zoo", lib.loc="~/R/win-library/3.3")
#library with functions for carrying out analyses on the extreme values of a process
library("extRemes", lib.loc="~/R/win-library/3.3")
#reading the precipitation data with header
precipdata=read.table("3.txt", header=TRUE)
colnames(precipdata)[1]= "Date"
colnames(precipdata)[2]= "Time"
#changing the format of the data
precipdata$DateTime <- as.POSIXct(paste(precipdata$Date, precipdata$Time), format="%m/%d/%Y %H:%M", tz = "UTC")
#assuming 1mm of precipitation
precipdata$Date <- NULL
precipdata$Time <- NULL
precipdata$Depth <- 0.1
#generating a time series with zoo
TS <- xts(precipdata$Depth, order.by=precipdata$DateTime)
#Compute Summary Statistics of zoo Objects(mean all values)
TS_aggreg <- aggregate(TS, identity, sum)
#Merging a set of the time by 1 minute
TS_by1min <- zoo(, seq(start(TS_aggreg), end(TS_aggreg), by = "1 min"))
#Merge both to one
TS_merged <- merge(TS_aggreg, TS_by1min)
#is.nais used If "any" then a row will be regarded as NA if it has any NAs. If "all" then a row will be regarded as NA only if all elements in the row are NA.
TS_merged[is.na(TS_merged)] <- 0
#Use dataframe to coercing objects of class "zoo" to other classes (only need the main values for dataframe)
finalTS.df <- data.frame(date=index(TS_merged), coredata(TS_merged))
#remove objects from global Environment that you will not use anymore
rm(precipdata,TS,TS_aggreg,TS_merged,TS_by1min)
#Apply moving average for user-defined durations using rollapply to get the cumulative depths with the sum function 
#Then find the annual maximum depth by using aggregate to split the data into subsets (each year) using the function max  
#Duration in minutes
tdur <- c(5,10,30,60,120,360,720,1440) 
#controler for creating the data frames
ctr=0
for (i in tdur){
  Depth <- finalTS.df[1]
  Depth$Depth <- rollapply(finalTS.df[2], i, sum, fill = 0, align='center')
  if (ctr < 1){
    MaxDepth <- aggregate( Depth$Depth ~ year(Depth$date), FUN = max)
  }
  else{
    MaxDepth[[ctr+2]] <- aggregate( Depth$Depth ~ year(Depth$date), FUN = max) [,-1]
  }
  ctr=ctr+1
}
colnames(MaxDepth) = c("Year", "15min", "30min", "1hr", "2hr", "4hr", "6hr", "12hr", "24hr")
#Convert to intensity mm/hr
MaxIntensity <- MaxDepth
#units in hours
tdurhr <- c(0.25, 0.5, 1, 2, 4, 6, 12, 24)
ctr=0
for (i in tdurhr){
  MaxIntensity[ctr+2] <- MaxDepth[ctr+2] / i
  ctr=ctr+1
}
#Fit An Extreme Value Distribution function to Data using Gumbel
rp <- c(2,5,10,50,100)
ctr=0
for(i in c(1, 2, 3, 4, 5, 6, 7, 8)){
  val = as.numeric(MaxIntensity[,1+i])
  fitgumbel <- fevd(val, MaxIntensity, type="Gumbel")
  rlev <- return.level(fitgumbel, return.period= rp)
  if (ctr < 1){
    gumbel <- t(data.frame(as.numeric(rlev)))
  }
  else{
    gumbel2  <- t(data.frame(as.numeric(rlev)))
    gumbel <- rbind(gumbel, gumbel2)
  }
  ctr=ctr+1
}
#Fit An Extreme Value Distribution function to Data using GEV
rp <- c(2,5,10,50,100)
ctr=0
for(i in c(1, 2, 3, 4, 5, 6, 7, 8)){
  val = as.numeric(MaxIntensity[,1+i])
  fitGEV <- fevd(val, MaxIntensity, type="GEV")
  rlevG <- return.level(fitGEV, return.period= rp)
  if (ctr < 1){
    GEV <- t(data.frame(as.numeric(rlevG)))
  }
  else{
    GEV2  <- t(data.frame(as.numeric(rlevG)))
    GEV <- rbind(GEV, GEV2)
  }
  ctr=ctr+1
}
#Sensitivity analysis to see if Gumbel or GEV is prefered(lowest value is)
AIC_GEV <- summary(fitGEV, silent=TRUE)
AIC_GEV <- c(AIC_GEV$AIC)
AIC_Gum <- summary(fitgumbel, silent=TRUE)
AIC_Gum <- c(AIC_Gum$AIC)
BIC_GEV <- summary(fitGEV, silent=TRUE)
BIC_GEV <- c(BIC_GEV$BIC)
BIC_Gum <- summary(fitgumbel, silent=TRUE)
BIC_Gum <- c(BIC_Gum$BIC)
if (AIC_Gum & BIC_Gum < AIC_GEV & BIC_GEV) {
  finaldist <- gumbel
} else {
  finaldist <- GEV
}
#Creating the final table
dur = as.data.frame(tdurhr)
finaldata <- data.frame(dur, finaldist)
options(digits=3)
colnames(finaldata) = c("x","2yr","5yr","10yr","50yr","100yr")
rownames(finaldata) = paste(round(tdurhr,digits = 1), "hr")
#Plotting the graph
idfgraph <- 
  ggplot(finaldata, aes(x)) + 
  geom_line(aes(y = `2yr`, colour = "2 years"), size = 0.8) + 
  geom_line(aes(y = `5yr`, colour = "5 years"), size = 0.8) + 
  geom_line(aes(y = `10yr`, colour = "10 years"), size = 0.8) + 
  geom_line(aes(y = `50yr`, colour = "50 years"), size = 0.8) + 
  geom_line(aes(y = `100yr`, colour = "100 years"), size = 0.8) + 
  ylab("Intensity [mm/hour]") +
  xlab("Duration [hours]") +
  labs(title = "IDF curves with precipitation data from 10/1994 to 10/2002") +
  labs(color = "Return Periods, years:") + 
  theme(legend.position = c(0.7, 0.7))
plot(idfgraph)




