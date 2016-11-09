## Assignment on Flash Floods to Develop IDF curve
## Submitted By: Liton Chandra Mazumder, EMFRM05 Student(2015-2017)

library(lubridate) #Functions to work with date-times and time-spans
library(xlsx) #Functions to read/write/format Excel
library(xts) #Functions for Time Series(TS)
library(highfrequency) #Functions for Time Series(TS)
library(matrixStats) 
library(ggplot2) #Functios for plotting
library(extRemes) #Functions for performing extreme value analysis

# At first set working directory
setwd("E:/Barcelona/FlashFloods/Practice/1 class")

#load the data
raindata = read.table("2.txt", header = TRUE)
colnames(raindata)[1]= "Date"
colnames(raindata)[2]= "Time"

# convert to date and time
raindata$DateTime<- as.POSIXct(paste(raindata$Date, raindata$Time), format="%m/%d/%Y %H:%M", tz = "UTC") 

# adding 0.1mm precipitation data in every time steps
raindata$Depth <- 0.1
# removing unecssary columns
raindata$Date <- NULL
raindata$Time <- NULL

#creation of the time series
require(xts)
ts <- xts(raindata$Depth, order.by = raindata$DateTime)
# aggreating the timeseries per minute
ts_aggr <- aggregate(ts, identity, sum)
# create regular timeseries using zoo
TS_full <- zoo(, seq(start(ts_aggr), end(ts_aggr), by = "1 min"))
# To merge irregular and regular timeseries
TS_merged <- merge(ts_aggr, TS_full)
# filling the NA with 0 values
TS_merged[is.na(TS_merged)] <- 0
# making data frame from the final timeseries
myTS.df <- data.frame(date=index(TS_merged), coredata(TS_merged))

#calculation of the cumulative depth (Sum) and the annual maximum depth for various time intervals from 5 minutes to 24 hours
Min5 <- rollapply(myTS.df$ts_aggr,5, sum, fill = 0, align='center')
Min5 <- data.frame(Date = myTS.df$date, Depth = Min5)
Min5 <- xts(Min5$Depth,myTS.df$date)
MaxMin5 <- apply.yearly(Min5, FUN = "max")
colnames(MaxMin5) <- "Min5"

Min15 <- rollapply(myTS.df$ts_aggr,15, sum, fill = 0, align='center')
Min15 <- data.frame(Date = myTS.df$date, Depth = Min15)
Min15 <- xts(Min15$Depth,myTS.df$date)
MaxMin15 <- apply.yearly(Min15, FUN = "max")
colnames(MaxMin15) <- "Min15"

Min30 <- rollapply(myTS.df$ts_aggr,30, sum, fill = 0, align='center')
Min30 <- data.frame(Date = myTS.df$date, Depth = Min30)
Min30 <- xts(Min30$Depth,myTS.df$date)
MaxMin30 <- apply.yearly(Min30, FUN = "max")
colnames(MaxMin30) <- "Min30"

Min45 <- rollapply(myTS.df$ts_aggr,45, sum, fill = 0, align='center')
Min45 <- data.frame(Date = myTS.df$date, Depth = Min45)
Min45 <- xts(Min45$Depth,myTS.df$date)
MaxMin45 <- apply.yearly(Min45, FUN = "max")
colnames(MaxMin45) <- "Min45"

Hour2 <- rollapply(myTS.df$ts_aggr,120, sum, fill = 0, align='center')
Hour2 <- data.frame(Date = myTS.df$date, Depth = Hour2)
Hour2 <- xts(Hour2$Depth,myTS.df$date)
MaxHour2 <- apply.yearly(Hour2, FUN = "max")
colnames(MaxHour2) <- "Hour2"

Hour6 <- rollapply(myTS.df$ts_aggr,360, sum, fill = 0, align='center')
Hour6 <- data.frame(Date = myTS.df$date, Depth = Hour6)
Hour6 <- xts(Hour6$Depth,myTS.df$date)
MaxHour6 <- apply.yearly(Hour6, FUN = "max")
colnames(MaxHour6) <- "Hour6"

Hour12 <- rollapply(myTS.df$ts_aggr,720, sum, fill = 0, align='center')
Hour12 <- data.frame(Date = myTS.df$date, Depth = Hour12)
Hour12 <- xts(Hour12$Depth,myTS.df$date)
MaxHour12 <- apply.yearly(Hour12, FUN = "max")
colnames(MaxHour12) <- "Hour12"

Hour24 <- rollapply(myTS.df$ts_aggr,1440, sum, fill = 0, align='center')
Hour24 <- data.frame(Date = myTS.df$date, Depth = Hour24)
Hour24 <- xts(Hour24$Depth,myTS.df$date)
MaxHour24 <- apply.yearly(Hour24, FUN = "max")
colnames(MaxHour24) <- "Hour24"

#Calculation of the maximum intensity(i = P/t) table for each duration

IntensityMin5 <- MaxMin5/(0.083)
IntensityMin15 <- MaxMin15/(0.25)
IntensityMin30 <- MaxMin30/(0.50)
IntensityMin45 <- MaxMin45/(0.75)
IntensityHour2 <- MaxHour2/(2)
IntensityHour6 <- MaxHour6/(6)
IntensityHour12 <- MaxHour12/(12)
IntensityHour24 <- MaxHour24/(24)

# Extreame Analysis by Using Gumbel method for each duration and return period

mlefit <- fevd(IntensityMin5, MaxMin5, type="Gumbel") #Activate for Gumbel method
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_5min <- as.numeric(rlmle)

mlefit <- fevd(IntensityMin15, MaxMin15, type="Gumbel")
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_15min <- as.numeric(rlmle)


mlefit <- fevd(IntensityMin30, MaxMin30, type="Gumbel") 
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_30min <- as.numeric(rlmle)

mlefit <- fevd(IntensityMin45, MaxMin45, type="Gumbel") 
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_45min <- as.numeric(rlmle)

mlefit <- fevd(IntensityHour2, MaxHour2, type="Gumbel") 
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_2hr <- as.numeric(rlmle)

mlefit <- fevd(IntensityHour6, MaxHour6, type="Gumbel")
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_6hr <- as.numeric(rlmle)


mlefit <- fevd(IntensityHour12, MaxHour12, type="Gumbel")
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_12hr <- as.numeric(rlmle)

mlefit <- fevd(IntensityHour24, MaxHour24, type="Gumbel") 
rlmle <- return.level(mlefit, conf = 0.05, return.period = c(2,5,10,20,50,100,200,500))
output_24hr <- as.numeric(rlmle)


gumbel<-t(matrix(c(output_5min,output_15min,output_30min,output_45min,output_2hr,output_6hr,output_12hr,output_24hr),ncol = 8))
colnames(gumbel) <- c(2,5,10,20,50,100,200,500)
gumbel<-as.data.frame(gumbel)

gumbel$t <- c(0.0833, 0.25, 0.5, 0.75, 2, 6, 12, 24)

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

