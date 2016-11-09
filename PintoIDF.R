library(lubridate) #package for working with date-times and time spans
library(xts) #package that provides uniform handling of different time-based data classes
library(highfrequency) #package for handling irregular time series
library(matrixStats) #provides funcions for operating with rows and columns in newun matrix
library(ggplot2) #package for graph creation
library(rlist) #For data manipulation with list objects
library(reshape2) #package for flexibility restructure and aggregating data
#Command for opening file
CodeFile <- file.choose()
TSdata = read.table(CodeFile, header=TRUE, skip = 1)
colnames(TSdata)[1]= "Date"
colnames(TSdata)[2]= "Time"
#Command for handling date and time
TSdata$DateTime <- as.POSIXct(paste(TSdata$Date, TSdata$Time), format="%m/%d/%Y %H:%M", tz = "UTC")
TSdata$Date <- NULL
TSdata$Time <- NULL
TSdata$Depth <- 0.1
RP <- c(2,5,10,50,100) #Return periods in years
duridf <- c(15,30,45,60,120,360,720,1440) #Duration in minutes
#Commands for making uniform TS
TS <- xts(TSdata[,-1], order.by=TSdata[,1])
TS_gr <- aggregate(TS, identity, sum)
TS_om <- zoo(, seq(start(TS_gr), end(TS_gr), by = "1 min"))
TS_merged <- merge(TS_gr, TS_om)
TS_merged[is.na(TS_merged)] <- 0
TS.df <- data.frame(date=index(TS_merged), coredata(TS_merged))
#creation of empty matrix for storing the results
IDFtable <- matrix(data=NA, nrow=length(duridf), ncol=length(RP))
#loop for calculation of moving average for each duration and RP
for (i in 1:length(duridf)) {
  newun <- TS.df[1]
  newun$Depth <- rollapply(TS.df[2], as.numeric(duridf[i]), sum, fill = NA, align='left')
  newun$Depth <- (newun$Depth / as.numeric(duridf[i])) * 60
  newun <- aggregate( newun$Depth ~ year(newun$date), FUN = max)
  aver <- mean(newun$TS_gr)
  standdev <- standdev(newun$TS_gr)
    for (j in 1:length(RP)) {
    freqf <- (-sqrt(6) / pi) * ((0.5772 + log((log(as.numeric(RP[j]) / (as.numeric(RP[j]) - 1))), base = exp(1))))
    intfin <- aver + standdev * freqf
    IDFtable[i,j] <- intfin
  }
}
#converting matrix to table
IDFtable <- as.data.frame(IDFtable)
IDFtablearranged <- c(as.numeric(as.character(t(duridf))))
IDFtablearranged <- as.data.frame(IDFtablearranged)
colnames(IDFtablearranged)[1] = "dur"
IDFtablearranged <- cbind(IDFtablearranged, IDFtable)
#rearranging final table in order to make a plot
mted <- melt(IDFtablearranged, id.vars="dur")
#plotting the results
IDFplot <- ggplot(mted, aes(dur, value, color = variable)) + geom_line() + geom_point()
plot(IDFplot) 
