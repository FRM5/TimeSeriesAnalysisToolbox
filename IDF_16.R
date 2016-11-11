library(lubridate) 
library(xts)
library(highfrequency)
library(matrixStats)
library(ggplot2)
library(forecast)
data = read.table("16.txt",header = FALSE, skip = 1)
colnames(data)[1]= "Date"
colnames(data)[2]= "Time"
data$DateTime <- as.POSIXct(paste(data$Date, data$Time), format="%m/%d/%Y %H:%M", tz = "GMT") 
data$Date <- NULL
data$Time <- NULL
data$Depth <- 0.1
#Timeseries processing
timeserie <- xts(data[,-1], order.by=data[,1])
aggregation <- aggregate(timeserie, identity, sum)
interval <- zoo(, seq(start(aggregation), end(aggregation), by = "1 min"))
TS_merge <- merge(aggregation,interval)
TS_merge[is.na(TS_merge)] <- 0
TS_frame <- data.frame(date=index(TS_merge), coredata(TS_merge))
#Designed duration
Duration <- c(5,10,30,60,120,360,720,1440)
#Designed Return period
RP <- c(5,10,50,100,200,500)
final <- matrix(data=NA, nrow=length(Duration), ncol=length(RP))
for (i in 1:length(Duration))
{
  TS_new <- TS_frame[1]
  TS_new$Depth <- ma(TS_frame[2],Duration[i], centre = FALSE)
  TS_new$Depth <- TS_new$Depth * 60
  Annually <- aggregate(TS_new$Depth ~ year(TS_new$date), FUN = max)
  Average <- mean(Annually$V1)
  SD <- sd(Annually$V1)
  alpha <- SD * sqrt(6)/pi
  beta <- Average - 0.5772*alpha
  for (j in 1:length(RP))
  {
    freqc <- (-sqrt(6) / pi) * ((0.5772 + log((log(as.numeric(RP[j]) / (as.numeric(RP[j]) - 1))), base = exp(1))))
    intensity <- Average + SD * freqc
    final[i,j] <- intensity
  }
}
final <- as.data.frame(final)
IDFplot <- as.data.frame(final)
IDFplot <- c(t((Duration)))
IDFplot <- as.data.frame(IDFplot)
colnames(IDFplot)[1] = "Duration"
IDFplot <- cbind(IDFplot, final)
colnames(IDFplot)[-1] <- paste(RP, "Years")
melt <- melt(IDFplot, id.vars="Duration")
finalplot <- ggplot(melt, aes(Duration, value, color = variable)) + geom_line() + geom_point() + ylab("Intensity, mm/hour") +  xlab("Duration, hours") +  labs(title = "Intensity-Duration-Frequency curves") + labs(color = "Return Periods, years:")
finalplot
