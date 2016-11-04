## Assignment on Flas Floods to make IDF curve
## Submitted By: Liton Chandra Mazumder, EMFRM05 Student(2015-2017)

library(lubridate) #Functions to work with date-times and time-spans
library(xlsx) #Functions to read/write/format Excel
library(xts) #Functions for Time Series(TS)
library(highfrequency) #Functions for Time Series(TS)
library(matrixStats) 
library(ggplot2) #Functios for plotting
library(extRemes) #Functions for performing extreme value analysis

setwd("E:/Barcelona/Drought/Practice/1 class")

raindata = read.table("2.txt", header = TRUE)
colnames(raindata)[1]= "Date"
colnames(raindata)[2]= "Time"

# convert to date and time
raindata$DateTime<- as.POSIXct(paste(raindata$Date, raindata$Time), format="%m/%d/%Y %H:%M", tz = "UTC") 

# adding 0.1mm precipitation data in every time steps
raindata$Depth <- 0.1
raindata$Date <- NULL
raindata$Time <- NULL

#creation of the time series
require(xts)
ts <- xts(raindata$Depth, order.by = raindata$DateTime)

ts_aggr <- aggregate(ts, identity, sum)
TS_full <- zoo(, seq(start(ts_aggr), end(ts_aggr), by = "1 min"))
TS_merged <- merge(ts_aggr, TS_full)
TS_merged[is.na(TS_merged)] <- 0
myTS.df <- data.frame(date=index(TS_merged), coredata(TS_merged))

#calculation of the cumulative depth (SUm) and the annual maximum depth
ctr=0
for(i in c(5, 15, 30, 60, 120, 240, 360, 760, 1440)){
  CMDepth <- myTS.df[1]
  CMDepth$Depth <- rollapply(myTS.df[2], i, sum, fill = 0, align='center')
  if (ctr < 1){
    MaxDepth <- aggregate( CMDepth$Depth ~ year(CMDepth$date), FUN = max) 
  }
  else{
    MaxDepth[[ctr+2]] <- aggregate( CMDepth$Depth ~ year(CMDepth$date), FUN = max) [,-1]
  }
  ctr=ctr+1
}

colnames(MaxDepth) = c("Year", "5min", "15min", "30min", "1hr", "2hr", "4hr", "6hr", "12hr", "24hr")

#Calculation of the maximum intensity table for each duration
MaxIntensity <- MaxDepth
MaxIntensity$`5min` <- MaxIntensity$`5min` / 0.0833
MaxIntensity$`15min` <- MaxIntensity$`15min` / 0.25
MaxIntensity$`30min` <- MaxIntensity$`30min` / 0.5
MaxIntensity$`2hr` <- MaxIntensity$`2hr` / 2
MaxIntensity$`6hr` <- MaxIntensity$`6hr` / 6
MaxIntensity$`12hr` <- MaxIntensity$`12hr` / 12
MaxIntensity$`24hr` <- MaxIntensity$`24hr` / 24

# Extreame Analysis by Using Gumbel method for each duration and return period

mlefit <- fevd(MaxIntensity$`5min`, MaxIntensity, type="Gumbel") #Activate for Gumbel method
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_5min <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`15min`, MaxIntensity, type="Gumbel")
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_15min <- as.numeric(rlmle)


mlefit <- fevd(MaxIntensity$`30min`, MaxIntensity, type="Gumbel") 
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_30min <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`2hr`, MaxIntensity, type="Gumbel") 
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_2hr <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`6hr`, MaxIntensity, type="Gumbel")
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_6hr <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`6hr`, MaxIntensity, type="Gumbel") 
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_6hr <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`12hr`, MaxIntensity, type="Gumbel")
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_12hr <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`24hr`, MaxIntensity, type="Gumbel") 
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_24hr <- as.numeric(rlmle)


gumbel<-t(matrix(c(output_5min,output_15min,output_30min,output_2hr,output_6hr,output_12hr,output_24hr),ncol = 7))
colnames(gumbel) <- c(2,5,10,20,50,100,200,500)
gumbel<-as.data.frame(gumbel)

gumbel$t <- c(0.0833, 0.25, 0.5, 2, 6, 12, 24)

#Plotting the graph for IDF Curve
ggplot(gumbel, aes(t, y = value, color = variable)) + 
  geom_line(aes(y = gumbel$`2`, color = "T2")) + geom_point(aes(y = gumbel$`2`, color = "T2")) +
  geom_line(aes(y = gumbel$`5`, color = "T5")) + geom_point(aes(y = gumbel$`5`, color = "T5")) +
  geom_line(aes(y = gumbel$`10`, color = "T10")) + geom_point(aes(y = gumbel$`10`, color = "T10")) +
  geom_line(aes(y = gumbel$`20`, color = "T20")) + geom_point(aes(y = gumbel$`20`, color = "T20")) +
  geom_line(aes(y = gumbel$`50`, color = "T50")) + geom_point(aes(y = gumbel$`50`, color = "T50")) +
  geom_line(aes(y = gumbel$`100`, color = "T100")) + geom_point(aes(y = gumbel$`100`, color = "T100")) +
  geom_line(aes(y = gumbel$`200`, color = "T200")) + geom_point(aes(y = gumbel$`200`, color = "T200")) +
  geom_line(aes(y = gumbel$`500`, color = "T500")) + geom_point(aes(y = gumbel$`500`, color = "T500")) +
  ylab("Intensity [mm/hr]") +
  xlab("Duration [hours]") +
  labs(title = "Intensity-Duration-Frequency Curve") +
  labs(color = "Return Periods")

