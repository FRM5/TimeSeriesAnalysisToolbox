# Flash Flood Assignment - Developing IDF Curve
# Adibtya Asyhari - Erasmus Mundus Flood Risk Management 2015/2017

library(xts)
library(lubridate)
library(highfrequency)
library(extRemes)
library(ggplot2)
library(zoo)
library(hydroTSM)

myFile <- file.choose()
mydata = read.table(myFile, header=TRUE)
colnames(mydata)[1]= "Date"
colnames(mydata)[2]= "Time"
mydata$DateTime <- as.POSIXct(paste(mydata$Date, mydata$Time), format="%m/%d/%Y %H:%M", tz = "UTC")
mydata$Depth <- 0.1

TS <- xts(mydata$Depth, order.by=mydata$DateTime)
TS_gr <- aggregate(TS, identity, sum)
TS_om <- zoo(, seq(start(TS_gr), end(TS_gr), by = "1 min"))
TS_merged <- merge(TS_gr, TS_om)
TS_merged[is.na(TS_merged)] <- 0
TS.df <- data.frame(date=index(TS_merged), coredata(TS_merged))

#write.table(myts.df, "D:/ON_GOING_PROJECTS/Assignment_FFDF/data.txt", sep="\t")

#calculating the cumulative depth and the annual maximum depth
ctr=0
for(i in c(15, 30, 60, 120, 240, 360, 760, 1440)){
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

colnames(MaxDepth) = c("Year", "15min", "30min", "1hr", "2hr", "4hr", "6hr", "12hr", "24hr")

#Calculating the annual maximum intensity table
MaxIntensity <- MaxDepth
ctr=0
for(i in c(0.25, 0.5, 1, 2, 4, 6, 12, 24)){
  MaxIntensity[ctr+2] <- MaxDepth[ctr+2] / i
  ctr=ctr+1
}

#Extremes analysis
ctr=0
for(i in c(1, 2, 3, 4, 5, 6, 7, 8)){
  val = as.numeric(MaxIntensity[,1+i])
  mlefit <- fevd(val, MaxIntensity, type="Gumbel")
  rlmle <- return.level(mlefit, conf = 0.05, return.period= c(2,5,10,50,100))
  if (ctr < 1){
    gumbel <- t(data.frame(as.numeric(rlmle)))
  }
  else{
    new_gumbel  <- t(data.frame(as.numeric(rlmle)))
    gumbel <- rbind(gumbel, new_gumbel)
  }
  ctr=ctr+1
}

#Creating the final table
dur = as.data.frame(c(0.25, 0.5, 1, 2, 4, 6, 12, 24))
output <- data.frame(dur, gumbel)
options(digits=3)
colnames(output) = c("x","2yr","5yr","10yr","50yr","100yr")
rownames(output) = c("15min", "30min", "1hr", "2hr", "4hr", "6hr", "12hr", "24hr")

#Plotting the graph
idfgraph <- 
  ggplot(output, aes(x)) + 
  geom_line(aes(y = `2yr`, colour = "2 years"), size = 1) + 
  geom_point(aes(y = `2yr`, colour = "2 years"), size = 1.6) + 
  geom_line(aes(y = `5yr`, colour = "5 years"), size = 1) + 
  geom_point(aes(y = `5yr`, colour = "5 years"), size = 1.6) + 
  geom_line(aes(y = `10yr`, colour = "10 years"), size = 1) + 
  geom_point(aes(y = `10yr`, colour = "10 years"), size = 1.6) + 
  geom_line(aes(y = `50yr`, colour = "50 years"), size = 1) + 
  geom_point(aes(y = `50yr`, colour = "50 years"), size = 1.6) + 
  geom_line(aes(y = `100yr`, colour = "100 years"), size = 1) + 
  geom_point(aes(y = `100yr`, colour = "100 years"), size = 1.6) +
  ylab("Intensity [mm/hour]") +
  xlab("Duration [hours]") +
  labs(title = "Intensity-Duration-Frequency Curves") +
  theme(panel.background = element_rect(fill = "white")) +
  theme(panel.grid.major = element_line(colour = "black", linetype = "dotted")) +
  theme(panel.grid.minor = element_line(colour = "black", linetype = "dotted")) +
  theme(legend.background = element_rect(fill= "white", size=0.5, 
                                         linetype="solid", colour ="black")) +
  labs(color = "Return Periods, years:") + 
  theme(legend.position = c(0.8, 0.68))

plot(idfgraph)