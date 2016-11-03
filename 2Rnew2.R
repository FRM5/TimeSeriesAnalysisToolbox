library(lubridate) #library for advance working with dates
library(xlsx) #library for writing data.frame to MS Excel file
library(xts) #library for working with time series
library(highfrequency) #library for handling irregular TS
library(matrixStats) #library for advanced working with matrix
library(ggplot2) #library for advanced plotting
library(extRemes)
setwd("E:/Barcelona/Drought/Practice/1 class")
mydata = read.table("2_cut.txt", header = TRUE)
colnames(mydata)[1]= "Date"
colnames(mydata)[2]= "Time"
mydata$DateTime<- as.POSIXct(paste(mydata$Date, mydata$Time), format="%m/%d/%Y %H:%M", tz = "UTC") # Problem with the Time Zone
mydata$Depth <- 0.1
mydata$Date <- NULL
mydata$Time <- NULL
require(xts)
ts <- xts(mydata$Depth, order.by = mydata$DateTime)
ts_aggr <- aggregate(ts, identity, sum)
TS_full <- zoo(, seq(start(ts_aggr), end(ts_aggr), by = "1 min"))
TS_merged <- merge(ts_aggr, TS_full)
TS_merged[is.na(TS_merged)] <- 0
TS.df <- data.frame(date=index(TS_merged), coredata(TS_merged))

#calculating the cumulative depth and the annual maximum depth
ctr=0
for(i in c(5, 15, 30, 60, 120, 240, 360, 760, 1440)){
  Depth <- TS.df[1]
  Depth$Depth <- rollapply(TS.df[2], i, sum, fill = 0, align='center')
  if (ctr < 1){
    MaxDepth <- aggregate( Depth$Depth ~ year(Depth$date), FUN = max) 
  }
  else{
    MaxDepth[[ctr+2]] <- aggregate( Depth$Depth ~ year(Depth$date), FUN = max) [,-1]
  }
  ctr=ctr+1
}

colnames(MaxDepth) = c("Year", "5min", "15min", "30min", "1hr", "2hr", "4hr", "6hr", "12hr", "24hr")

#generation of max intensity table
MaxIntensity <- MaxDepth
MaxIntensity$`5min` <- MaxIntensity$`5min` / 0.0833
MaxIntensity$`15min` <- MaxIntensity$`15min` / 0.25
MaxIntensity$`30min` <- MaxIntensity$`30min` / 0.5
MaxIntensity$`2hr` <- MaxIntensity$`2hr` / 2
MaxIntensity$`6hr` <- MaxIntensity$`6hr` / 6
MaxIntensity$`12hr` <- MaxIntensity$`12hr` / 12
MaxIntensity$`24hr` <- MaxIntensity$`24hr` / 24

# GEV or Gumbel method

mlefit <- fevd(MaxIntensity$`5min`, MaxIntensity, type="Gumbel") #Activate for Gumbel method
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
outrp_5 <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`15min`, MaxIntensity, type="Gumbel") #Activate for Gumbel method
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
outrp_15 <- as.numeric(rlmle)


mlefit <- fevd(MaxIntensity$`30min`, MaxIntensity, type="Gumbel") #Activate for Gumbel method
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
outrp_30 <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`2hr`, MaxIntensity, type="Gumbel") #Activate for Gumbel method
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
outrp_2hr <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`6hr`, MaxIntensity, type="Gumbel") #Activate for Gumbel method
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
outrp_6hr <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`6hr`, MaxIntensity, type="Gumbel") #Activate for Gumbel method
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
outrp_6hr <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`12hr`, MaxIntensity, type="Gumbel") #Activate for Gumbel method
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
outrp_12hr <- as.numeric(rlmle)

mlefit <- fevd(MaxIntensity$`24hr`, MaxIntensity, type="Gumbel") #Activate for Gumbel method
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
outrp_24hr <- as.numeric(rlmle)


gumbel<-t(matrix(c(outrp_5,outrp_15,outrp_30,outrp_2hr,outrp_6hr,outrp_12hr,outrp_24hr),ncol = 7))
colnames(gumbel) <- c(2,5,10,20,50,100,200,500)
gumbel<-as.data.frame(gumbel)

gumbel$t <- c(0.0833, 0.25, 0.5, 2, 6, 12, 24)

# Plot graph manually
ggplot(gumbel, aes(t, y = value, color = variable)) + 
  geom_line(aes(y = gumbel$`2`, col = "T2")) + geom_point(aes(y = gumbel$`2`, col = "T2")) +
  geom_line(aes(y = gumbel$`5`, col = "T5")) + geom_point(aes(y = gumbel$`5`, col = "T5")) +
  geom_line(aes(y = gumbel$`10`, col = "T10")) + geom_point(aes(y = gumbel$`10`, col = "T10")) +
  geom_line(aes(y = gumbel$`20`, col = "T20")) + geom_point(aes(y = gumbel$`20`, col = "T20")) +
  geom_line(aes(y = gumbel$`50`, col = "T50")) + geom_point(aes(y = gumbel$`50`, col = "T50")) +
  geom_line(aes(y = gumbel$`100`, col = "T100")) + geom_point(aes(y = gumbel$`100`, col = "T100")) +
  geom_line(aes(y = gumbel$`200`, col = "T200")) + geom_point(aes(y = gumbel$`200`, col = "T200")) +
  geom_line(aes(y = gumbel$`500`, col = "T500")) + geom_point(aes(y = gumbel$`500`, col = "T500")) +
  ylab("Intensity [mm/hr]") +
  xlab("Duration [hours]") +
  labs(title = "Intensity-Duration-Frequency Curve") +
  labs(color = "Return Periods")

