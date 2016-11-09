library(lubridate) #library for advance working with dates
library(xts) #library for working with time series
library(highfrequency) #library for handling irregular TS
library(matrixStats) #library for advanced working with matrix
library(ggplot2) #library for advanced plotting
library(rlist) #toolbox for Non-Tabular Data Manipulation
library(reshape2) #flexibly Reshape Data
library(forecast) #methods and tools for displaying and analysing univariate time series
library(gridExtra) #library for multiple plotting
myFile <- file.choose() #useful function without setting working directory
message("Reading data from the file")
mydata = read.table(myFile, header=FALSE, skip = 1)
colnames(mydata)[1]= "Date"
colnames(mydata)[2]= "Time"
mydata$DateTime <- as.POSIXct(paste(mydata$Date, mydata$Time), format="%m/%d/%Y %H:%M") #merging date and time
mydata$Date <- NULL
mydata$Time <- NULL
mydata$Depth <- 0.1
message("Creating regular time series from provided data")
#in the following part of code we create a regular TS with time step of 1 min.
myts <- xts(mydata[,-1], order.by=mydata[,1])
myts_gr <- aggregate(myts, identity, sum) #if we have records for same minute, we will sum it, in order to have
#one record for one minute. Extremely important for creation of regular TS
myts_om <- zoo(, seq(start(myts_gr), end(myts_gr), by = "1 min"))
myts_merged <- merge(myts_gr, myts_om) #merging irrefular TS with depth data and regular TS with proper time steps
myts_merged[is.na(myts_merged)] <- 0
myts.df <- data.frame(date=index(myts_merged), coredata(myts_merged)) #converting regular TS to data.frame
#asking user for RP and Durations. All the calculations and plotting will be done automatically
askrt <- readline(prompt="Enter return periods for IDF curves computation (coma separated list) :")
askdur <- readline(prompt="Enter durations for IDF curves computation (coma separated list) :")
askrt <- strsplit(askrt, ",")[[1]] #converting user input into kind of Vector with nubers
askdur <- strsplit(askdur, ",")[[1]]
final <- matrix(data=NA, nrow=length(askdur), ncol=length(askrt))
kscheck <- matrix(data=NA, nrow=1, ncol=length(askdur))
askslks <- readline(prompt="Enter significance level for Kolmogorov-Smirnov test (one value with dot separator, like 0.5) :")
askslks <- as.numeric(askslks)
message("Performing annual extremes analysis")
#the following line of the code contains for loop cycle for automatization of work
#this code can calculate as much durations and RP as user wants. 
#Finally MA function was used for improving the performance.
#It's pointless to explain every line here as it is simple nested loop
for (i in 1:length(askdur)) {
  anextr <- myts.df[1]
  anextr$Depth <- ma(myts.df[2], as.numeric(askdur[i]), centre = FALSE)
  anextr$Depth <- anextr$Depth * 60
  anextr <- aggregate( anextr$Depth ~ year(anextr$date), FUN = max)
  aver <- mean(anextr$V1)
  sd <- sd(anextr$V1)
  alpha <- sd * sqrt(6) / pi
  beta <- aver - 0.5772 * alpha
  for (j in 1:length(askrt)) {
    freqf <- (-sqrt(6) / pi) * ((0.5772 + log((log(as.numeric(askrt[j]) / (as.numeric(askrt[j]) - 1))), base = exp(1))))
    intfin <- aver + sd * freqf
    final[i,j] <- intfin
  }
  #creation of empty matrix for results from K-S test
  fitmat <- matrix(data = NA, nrow = nrow(anextr), ncol = 7)
  fitmat[,2] <- as.matrix(anextr[2])
  fitmat <- fitmat[order(fitmat[,2]),]
  fitmat[,1] <- 1:nrow(anextr)
  for (k in 1:nrow(fitmat)){
    fitmat[k,3] <- (fitmat[k,2] - aver)^2
    fitmat[k,4] <- (fitmat[k,1] - 0.44) / (nrow(fitmat) + 0.12)
    fitmat[k,5] <- (exp(-exp(-(fitmat[k,2] - beta)/alpha)))
    fitmat[k,6] <- abs(fitmat[k,5] - fitmat[k,4])
    if (k > 1){
      fitmat[k,7] <- abs(fitmat[k,5] - fitmat[k-1,4])
    }else{
      fitmat[k,7] <- 0
    }
    kscheck[,i] <- max(fitmat[,6], fitmat[,7])
  }
  if (kscheck[,i] < sqrt(-0.5*log(askslks/2))/sqrt(nrow(anextr))){ #equation for critical value of K-S test
    message("Data for duration ", askdur[i],"min pass Kolmogorov-Smirnov test")
  }else{
    message("Data for duration ", askdur[i],"min DID NOT pass Kolmogorov-Smirnov test")
  }
}
#filling final table with the IDF results
final <- as.data.frame(final)
finaln <- c(as.numeric(as.character(t(askdur))))
finaln <- as.data.frame(finaln)
colnames(finaln)[1] = "dur"
finaln <- cbind(finaln, final)
colnames(finaln)[1] = "Duration"
colnames(finaln)[-1] <- paste(askrt, "Years")
melted <- melt(finaln, id.vars="Duration")
#the following part include calculations of A-B hyetograph
#moving average was used again for calculation of cumulative intensities
answercalc <- readline(prompt="Do you want to calculate rainfall hyetograph? (type Y for yes and N for no) :")
if (answercalc == "Y"){
  rt <- readline(prompt="Enter return period (years): ")
  dur <- readline(prompt="Enter duration (minutes): ")
  inc <- readline(prompt="Enter increments (minutes): ")
  rt <- as.numeric(rt)
  dur <- as.numeric(dur)
  inc <- as.numeric(inc)
  durhg <- matrix(data = NA, nrow = 1, ncol = (dur/inc))
  for (i in 1:(dur/inc)){
    durhg[,i] <- inc * i
  }
  finalhg <- matrix(data = NA, nrow = length(durhg), ncol = 4)
  finalhg[,1] <- t(durhg)
  freqf2 <- (-sqrt(6) / pi) * ((0.5772 + log((log(rt / (rt - 1))), base = exp(1))))
  hyetoplot <- matrix(data = NA, nrow = length(durhg), ncol = 2)
  hyetoplot[1,1] <- round(length(durhg)/2)
  for (i in 1:length(durhg)) {
    hyetdur <- myts.df[1]
    hyetdur$Depth <- ma(myts.df[2], as.numeric(durhg[i]), centre = FALSE)
    hyetdur$Depth <- hyetdur$Depth * 60
    hyetdur <- aggregate( hyetdur$Depth ~ year(hyetdur$date), FUN = max)
    aver2 <- mean(hyetdur$V1)
    sd2 <- sd(hyetdur$V1)
    intfinhg <- aver2 + sd2 * freqf2
    finalhg[i,2] <- intfinhg
    finalhg[i,3] <- finalhg[i,2] / (60 / finalhg[i,1])
    if (i > 1){
      finalhg[i,4] <- finalhg[i,3] - finalhg[i-1,3]
    }else{
      finalhg[i,4] <- finalhg[i,3]
    }
    if (i == i){
      message("Intensity for increment has been calculated (", durhg[i],"min)")
    }
  }
  hyetoplot[,2] <- sort(finalhg[,4], decreasing = TRUE)
  for ( i in 2:length(durhg)){
    if (i%%2 == 0){
      hyetoplot[i,1] <- hyetoplot[1,1] + i/2
    }else{
      hyetoplot[i,1] <- hyetoplot[1,1] - (i-1)/2
    }
  }
  hyetoplot <- hyetoplot[order(hyetoplot[,1]),]
  hyetoplot <- as.data.frame(hyetoplot)
  
}
hyetoplot[,1] <- finalhg[,1]
colnames(hyetoplot) = c("Duration", "Depth, mm")
#plotting IDF and A-B hyetograph graphs
finalplot <- ggplot(melted, aes(Duration, value, color = variable)) + geom_line() + geom_point() + ylab("Intensity, mm/hour") +  xlab("Duration, hours") +  labs(title = "Intensity-Duration-Frequency curves") + labs(color = "Return Periods, years:")
ggplhyeto <- ggplot(data=hyetoplot, aes(x=Duration, y=`Depth, mm`, fill = `Depth, mm`)) +  geom_bar(stat="identity") + scale_fill_continuous(low="blue", high="green") + ylab("Depth, mm") +  xlab("Increments, minutes") +  labs(title = "Alternate Block Hyetograph")
answerplot <- readline(prompt="Do you want to plot IDF graph and Alternating Block hyetograph? (type Y for yes and N for n) :")
if (answerplot == "Y"){
  plotgrphs <- grid.arrange(finalplot, ggplhyeto)
}
answerfinish <- readline(prompt="Do you want to clean Environment for better reading the tables? (type Y for yes and N for n) :")
if (answerfinish == "Y"){
  rm(finalhg, durhg, anextr, sd2, sd, rt,myts_om, myts_merged, myts_gr, myFile, k, j, intfinhg, intfin, i, freqf2, freqf, dur, beta, aver2, aver, askslks, answercalc, answerplot, answerfinish, alpha, myts, melted, hyetdur, final)
  colnames(fitmat) = c("Ranking", "Intensity", "Squared res.", "Dringorten par.", "CDF Gumbel", "KS d1", "KS d2")
}
answerfinish2 <- readline(prompt="Do you want to clear Environment? (type Y for yes and N for n) :")
if (answerfinish2 == "Y"){
  rm(list = ls())
}
message("Program finished executing, please re-run the program")
#Finally the end. Actually only for computing IDF curves with user input I can do this application in 50 lines
#But as I include Kolmogorov-Smirnov test and Synthetic precipitation...and some comments, now I have 158 :D