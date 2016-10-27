library(lubridate) #library for advance working with dates
library(xts) #library for working with time series
library(highfrequency) #library for handling irregular TS
library(matrixStats) #library for advanced working with matrix
library(ggplot2) #library for advanced plotting
library(rlist)
library(reshape2)
askrt <- readline(prompt="Enter return periods for IDF curves computation (coma separated list) :")
askrt <- strsplit(askrt, ",")[[1]] 
askdur <- readline(prompt="Enter durations for IDF curves computation (coma separated list) :")
askdur <- strsplit(askdur, ",")[[1]]
myFile <- file.choose()
mydata = read.table(myFile, header=TRUE)
colnames(mydata)[1]= "Date"
colnames(mydata)[2]= "Time"
mydata$DateTime <- as.POSIXct(paste(mydata$Date, mydata$Time), format="%m/%d/%Y %H:%M", tz = "UTC")
mydata$Date <- NULL
mydata$Time <- NULL
mydata$Depth <- 0.1
myts <- xts(mydata[,-1], order.by=mydata[,1])
myts_gr <- aggregate(myts, identity, sum)
myts_om <- zoo(, seq(start(myts_gr), end(myts_gr), by = "1 min"))
myts_merged <- merge(myts_gr, myts_om)
myts_merged[is.na(myts_merged)] <- 0
myts.df <- data.frame(date=index(myts_merged), coredata(myts_merged))
final <- matrix(data=NA, nrow=length(askdur), ncol=length(askrt))
for (i in 1:length(askdur)) {
  a <- myts.df[1]
  a$Depth <- rollapply(myts.df[2], as.numeric(askdur[i]), sum, fill = NA, align='left')
  a$Depth <- (a$Depth / as.numeric(askdur[i])) * 60
  a <- aggregate( a$Depth ~ year(a$date), FUN = max)
  aver <- mean(a$myts_gr)
  sd <- sd(a$myts_gr)
    for (j in 1:length(askrt)) {
      freqf <- (-sqrt(6) / pi) * ((0.5772 + log((log(as.numeric(askrt[j]) / (as.numeric(askrt[j]) - 1))), base = exp(1))))
      intfin <- aver + sd * freqf
      final[i,j] <- intfin
    }
}
final <- as.data.frame(final)
finaln <- c(as.numeric(as.character(t(askdur))))
finaln <- as.data.frame(finaln)
colnames(finaln)[1] = "dur"
finaln <- cbind(finaln, final)
melted <- melt(finaln, id.vars="dur")
finalplot <- ggplot(melted, aes(dur, value, color = variable)) + geom_line() + geom_point()
answerplot <- readline(prompt="Do you want to plot current IDF graph? (type Y for yes and N for n) :")
if (answerplot == "Y"){
  plot(finalplot) 
}





