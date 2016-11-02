library(zoo) # timeseries
library(xts) # timeseries
library(highfrequency) # timeseries
library(lubridate)
library(ggplot2) # plotting
library(reshape2)

setwd("D:/OneDrive/FlashFloods/Assignment") # set working directory
MyData = read.table("9.txt", header = TRUE) # read the file

# formating the data
colnames(MyData)[1] = "Date" # change column name
colnames(MyData)[2] = "Time" # change column name
MyData$Date_Time <- as.POSIXct(paste(MyData$Date, MyData$Time), format="%m/%d/%Y %H:%M") # change class to date
MyData$Depth <- 0.1 # create column with depth
MyData$Date <- NULL # delete unnecessary column
MyData$Time <- NULL # delete unnecessary column
ts <- xts(MyData[,-1], order.by=MyData[,1]) # create timeseries
ts_aggr <- aggregate(ts, identity, FUN= sum) # aggreate timeseries per minute
ts_full <- zoo(, seq(start(ts_aggr), end(ts_aggr), by = "1 min")) # create regular timeseries
z <- merge(ts_aggr, ts_full) # merge irregular and regular timeseries
z[is.na(z)] <- 0 # fill the NA with 0 values
ts.frame <- data.frame(date=index(z), coredata(z)) # make data frame from the final timeseries
colnames(ts.frame) <- c("Date_Time", "Depth") # change column names

dur <- c(5, 15, 30, 60, 120, 360, 720, 1440, 2880) # pick the duration of the time intervals (windows)
rp <- c(2, 5, 10, 20, 25, 50, 100, 500) # pick the return periods
# For Duration, proper analysis needs to be done for the related catchments. The durations need to fit the catchment's concentration time.
# Also small values, below ~4 hours can be used for pluvial floods in urban areas, since they affect the sewage and drainage system 
# For Return Periods, the maximum value that can be givem is based n the number of years of the data set. For the case of ~25 years, for return
# periods above 25 years the error is extremely high, but they are produced just for research purposes.

IDF <- matrix(data = NA, nrow = length(dur), ncol = length(rp)) # create emply IDF matrix
ctr <- 0 # controler for creating the below data frames

# loop for creation of IDF and Maximum annual intesities. First loop is for each time window, second loop for each return period
for (i in 1:length(dur)) {   
  a <- ts.frame[1]
  a$Depth <- rollapply(ts.frame[2], as.numeric(dur[i]), sum, fill = NA, align='left')
  a$Depth <- (a$Depth / as.numeric(dur[i])) * 60
  a <- aggregate( a$Depth ~ year(a$Date_Time), FUN = max)
  a_sort <- a[with(a, order(Depth)),]
  aver <- mean(a$Depth)
  sd <- sd(a$Depth)
  if (ctr < 1){
    Maxima <- a
    Maxima_sort <- data.frame(a_sort[,-1])
    Koltest <- data.frame(rank=index(a))
    Koltest$Gr <- (Koltest$rank-0.44)/(length(Koltest$rank)+0.12)
    Koltest[[ctr+3]] <- exp(-exp(-(a_sort$Depth-(aver-0.5772*sd*sqrt(6)/pi))/(sd*sqrt(6)/pi))) 
  }
  else{
    Maxima[[ctr+2]] <- a[,-1]
    Maxima_sort[[ctr+1]] <- a_sort[,-1]
    Koltest[[ctr+3]] <- exp(-exp(-(a_sort$Depth-(aver-0.5772*sd*sqrt(6)/pi))/(sd*sqrt(6)/pi)))
  }
  ctr=ctr+1
  for (j in 1:length(rp)) {
    IDF[i,j] <- -log(-log(1-1/rp[j]))*sd*sqrt(6)/pi+aver-0.5772*sd*sqrt(6)/pi
  }
}
colnames(Maxima)[1] <- "Year"
colnames(Maxima)[-1] <- paste(dur, "Min")
colnames(Maxima_sort) <- paste(dur, "Min")
colnames(IDF) <- rp
IDF <- as.data.frame(IDF)
IDF_final <- cbind(as.data.frame(dur/60), IDF)
colnames(IDF_final)[1] = "Duration"
colnames(IDF_final)[-1] = paste(rp, "Years")
melted <- melt(IDF_final, id.vars="Duration")
IDFplot <- ggplot(melted, aes(Duration, value, color = variable)) + geom_line() + geom_point() +
  ylab("Intensity [mm/hr]") + xlab("Duration [hr]") + labs(title = "IDF") + labs(color = "Return Periods")
plot(IDFplot) 

# For goodtness of fit, Kolmogorov test was used. The related created data are: Koltest, Koltest2, Test, Test2, MaxTest, MaxKol and the final
# one is GootFit. For this table, the values need to be checked manually, with the maximum allowed values of the test, based on the significance
# and the number of years of the initial data set. Based on the wbove, the proposed Gumbel distribution fits the data properly. Koltest is the 
# differences of the Gringorten and the Gumbel values of the same line. Koltest2 is the differences of Gringorten with the Gumbel values of one
# row below. From these 2 tables the maximum difference is assigned to table Goodfit

colnames(Koltest)[1] <- "Index"
colnames(Koltest)[2] <- "Gringorten"
colnames(Koltest)[-(1:2)] <- dur
Koltest2 <- data.frame(Koltest[-1,-2])
Test <- data.frame(Koltest$Index)
Test2 <- data.frame(Koltest2$Index)
ctr <- 2  
for (i in 1:length(dur)){ 
  Test[[ctr]] <- abs(Koltest$Gringorten - Koltest[i+2])
  Test2[[ctr]] <- abs(Koltest$Gringorten - Koltest2[i+1])
  ctr <- ctr+1
} 
Test <- Test[-1]
Test2 <- Test2[-1]
colnames(Test) <- paste(dur, "Min")
colnames(Test2) <- paste(dur, "Min")
MaxTest2 <- apply(Test2-0, 2, max)
MaxTest <- apply(Test-0, 2, max)
MaxKol <- rbind(MaxTest2, MaxTest)
Goodfit <- t(as.data.frame(apply(MaxKol, 2, max)))
rownames(Goodfit) <- "MaxDif"
colnames(Goodfit) <- paste(dur, "Min")

# rm(IDF, a, a_sort, ctr, i, j, melted, aver, sd, z,ts_aggr, ts_full, MaxTest, MaxTest2, MaxKol, Koltest, Koltest2)