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

#Creating the uniform time series
TS <- xts(mydata$Depth, order.by=mydata$DateTime)
TS_gr <- aggregate(TS, identity, sum)
TS_om <- zoo(, seq(start(TS_gr), end(TS_gr), by = "1 min"))
TS_merged <- merge(TS_gr, TS_om)
TS_merged[is.na(TS_merged)] <- 0
TS.df <- data.frame(date=index(TS_merged), coredata(TS_merged))

#Exporting the data
#write.table(myts.df, "D:/ON_GOING_PROJECTS/Assignment_FFDF/data.txt", sep="\t")

#Calculating the cumulative depth and the annual maximum depth
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

#Calculating the annual maximum intensity
MaxIntensity <- MaxDepth
ctr=0
for(i in c(0.25, 0.5, 1, 2, 4, 6, 12, 24)){
  MaxIntensity[ctr+2] <- MaxDepth[ctr+2] / i
  ctr=ctr+1
}

#Goodness of fit test and the extremes analysis
ctr=0
for(i in c(1, 2, 3, 4, 5, 6, 7, 8)){
  val = as.numeric(MaxIntensity[,1+i])
  fitgev <- fevd(val, MaxIntensity, type="GEV")
  fitgmb <- fevd(val, MaxIntensity, type="Gumbel")
  fit1 <- summary(fitgev)
  fit2 <- summary(fitgmb)
  #Choosing the best fit by comparing the value of BIC, the lower value the better
  if (fit1$BIC < fit2$BIC){
    mlefit <- fitgev
  }
  else{
    mlefit <- fitgmb
  }
  rlmle <- return.level(mlefit, conf = 0.05, return.period= c(2,5,10,50,100))
  if (ctr < 1){
    extreme <- t(data.frame(as.numeric(rlmle)))
  }
  else{
    new_extreme  <- t(data.frame(as.numeric(rlmle)))
    extreme <- rbind(extreme, new_extreme)
  }
  ctr=ctr+1
}

#Creating the final table
dur = as.data.frame(c(0.25, 0.5, 1, 2, 4, 6, 12, 24))
output <- data.frame(dur, extreme)
options(digits=3)
colnames(output) = c("x","2yr","5yr","10yr","50yr","100yr")
rownames(output) = c("15min", "30min", "1hr", "2hr", "4hr", "6hr", "12hr", "24hr")

#Plotting the IDF curve
idfgraph <- 
  ggplot(output, aes(x)) + 
  geom_line(aes(y = `2yr`, colour = "2 years")) + 
  geom_point(aes(y = `2yr`, colour = "2 years")) + 
  geom_line(aes(y = `5yr`, colour = "5 years")) + 
  geom_point(aes(y = `5yr`, colour = "5 years")) + 
  geom_line(aes(y = `10yr`, colour = "10 years")) + 
  geom_point(aes(y = `10yr`, colour = "10 years")) + 
  geom_line(aes(y = `50yr`, colour = "50 years")) + 
  geom_point(aes(y = `50yr`, colour = "50 years")) + 
  geom_line(aes(y = `100yr`, colour = "100 years")) + 
  geom_point(aes(y = `100yr`, colour = "100 years")) +
  ylab("Intensity [mm/hour]") +
  xlab("Duration [hours]") +
  labs(title = "Intensity-Duration-Frequency Curves")

plot(idfgraph)

#Creating the synthetic precipitation with the return period of 2 years
idf2 <- data.frame(output$x, output$`2yr`)
colnames(idf2) = c("x","y")
idf <- nls(y ~ c/(f+x^e), start=list(c=1, f=1, e=1),data=idf2)
int <-predict(idf, newdata = data.frame(x=c(seq(from = 0.25, to = 6, by = 0.25))))
synth <- data.frame(c(seq(from = 0.25, to = 6, by = 0.25)))
synth[2] <- int
synth[3] <- synth[1] * synth[2]
synth[4] <- c(synth[1,3],diff(as.matrix(synth[3])))

P <- as.matrix(unlist(synth[4]))
med <- max(P)
val <- which.max(P)
Q <- P[-val]

if (length(P) %% 2 != 0){ #check if function is odd
  t <- sort(which(index(Q) %% 2 == 0),decreasing = TRUE)
  h <- sort(which(index(Q) %% 2 != 0),decreasing = FALSE)
} else {
  h <- sort(which(index(Q) %% 2 == 0),decreasing = TRUE)
  t <- sort(which(index(Q) %% 2 != 0),decreasing = FALSE)
}

synth[5] <- c(Q[h],med,Q[t])

colnames(synth) = c("Duration","Intensity","Cumulative Depth","Incremental Depth","Precipitation")

#Plotting the hyetograph
e <- as.numeric(as.matrix(synth[5]))
barplot(e, main="2 Years - 6 hours Synthetic Precipitation", xlab="Time", ylab="Precipitation")