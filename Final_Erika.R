#Intensity-Duration-Frequency (IDF) code
#Created by: Erika Roxana Landaverde. EMFRM 15/16
#IMPORTANT:Column 1= Date, Column 2 = Time
#Durations and return periods are already established as it can take a long time for short durations

require(zoo)
require(xts)
require(hydroTSM)
require(lubridate)
require(ggplot2)
require(reshape)

message("Intensity-Duration-Frequency (IDF) code
================================================================================")
message("**Info: This code consideres durations of 15min, 30 min, 1 hr, 3 hrs, 8 hrs, 12 hrs and 1 day.")
message("Please keep in mind that the input data should be a .txt, with column1=date 
and column2=time (Including headers)**")
invisible(readline(prompt="Press [enter] to continue"))
#Read data from given text file
message("Please choose rainfall data file")
mydata <- read.table(file.choose(TRUE),header=TRUE)
message("*Wait for data to be processed*")

#Setting up headers for tables
colnames(mydata)[2]= "Time"
colnames(mydata)[1]= "Date"

#Creating a column with date and time 
mydata$Date <- strptime(paste(mydata$Date,mydata$Time), "%m/%d/%Y %H:%M" , tz = "UTC")

#Adding a depth column to the table and eliminating the extra column
mydata$depth <- 0.1
mydata$Time <- NULL
final_data <-subset(mydata,diff(mydata$Date) >0) #Removes overlapping times

#Convert into time series
ts <- xts(mydata [,-1], order.by=mydata[,1])
#Convert into 1 minute uniform time series
new_ts <- period.apply(ts, endpoints(ts, "minutes", 1), sum)
new2_ts <- zoo(, seq(start(new_ts), end(new_ts), by = "1 min"))
final_ts <- merge(new_ts, new2_ts)
final_ts[is.na(final_ts)] <- 0
final <- data.frame(date=index(final_ts), coredata(final_ts)) 
colnames(final) <- c("Date", "Depth") 

message("============================================================================================
Calculating Annual Maximum precipitations. This process might take several minutes...
============================================================================================")

#NOTE: The following lines (moving average, gumbel and so on) can be improved by doing a loops...I'm just not familiar with the function.
dur_15max <-final[1]
dur_15max$Depth <- rollapply(final[2],width=15, FUN=sum, fill=NA, align = 'center')
dur_15max <- aggregate(dur_15max$Depth ~ year(dur_15max$Date), FUN=max)

dur_30max <-final[1]
dur_30max$Depth <- rollapply(final[2],width=30, FUN=sum, fill=NA, align = 'center')
dur_30max <- aggregate(dur_30max$Depth ~ year(dur_30max$Date), FUN=max)

dur_60max <-final[1]
dur_60max$Depth <- rollapply(final[2],width=60, FUN=sum, fill=NA, align = 'center')
dur_60max <- aggregate(dur_60max$Depth ~ year(dur_60max$Date), FUN=max)

dur_180max <-final[1]
dur_180max$Depth <- rollapply(final[2],width=180, FUN=sum, fill=NA, align = 'center')
dur_180max <- aggregate(dur_180max$Depth ~ year(dur_180max$Date), FUN=max)

dur_480max <-final[1]
dur_480max$Depth <- rollapply(final[2],width=480, FUN=sum, fill=NA, align = 'center')
dur_480max <- aggregate(dur_480max$Depth ~ year(dur_480max$Date), FUN=max)

dur_720max <-final[1]
dur_720max$Depth<- rollapply(final[2],width=720, FUN=sum, fill=NA, align = 'center')
dur_720max <- aggregate(dur_720max$Depth ~ year(dur_720max$Date), FUN=max)

dur_1440max <-final[1]
dur_1440max$Depth <- rollapply(final[2],width=1440, FUN=sum, fill=NA, align = 'center')
dur_1440max <- aggregate(dur_1440max$Depth ~ year(dur_1440max$Date), FUN=max)

message("================================================================================
Calculating IDFs...
================================================================================")
#Gumbel Method.
#For this method the following data is needed
  #Maximum Precipitation Table
  dur <- c("15 min", "30 min", "60 min", "180 min", "480 min", "720 min", "1440 min")
  RT <- c("5 years", "10 years", "25 years", "50 years", "100 years")
  rt <- c(5,10,25,50,100)
  dura <- c(15,30,60,180,480,720,1440)
  Max_Pre <- data.frame(dur_15max$Depth, dur_30max$Depth, dur_60max$Depth, dur_180max$Depth, dur_480max$Depth, dur_720max$Depth, dur_1440max$Depth)
  colnames(Max_Pre) <- dur
  rownames(Max_Pre) <- dur_15max[,1]
  #Gumbel Parameters
  Gumbel <- data.frame(matrix(data=NA,ncol=length(dur), nrow=2))
  colnames(Gumbel) <- dur
  rownames(Gumbel) <- c("standard deviation", "Annual Average")
  #Standard deviation
  Gumbel[1,1] <- sd(Max_Pre$`15 min`)
  Gumbel[1,2] <- sd(Max_Pre$`30 min`)
  Gumbel[1,3] <- sd(Max_Pre$`60 min`)
  Gumbel[1,4] <- sd(Max_Pre$`180 min`)
  Gumbel[1,5] <- sd(Max_Pre$`480 min`)
  Gumbel[1,6] <- sd(Max_Pre$`720 min`)
  Gumbel[1,7] <- sd(Max_Pre$`1440 min`)
  #Annual Average
  Gumbel[2,] <- colMeans(Max_Pre, na.rm = FALSE, dims =1)
  #K parameter
  K <- matrix(data=NA, nrow=1, ncol=5)
  colnames(K) <- RT
  K[1,] <- (-sqrt(6) / pi) * ((0.5772 + log(log(rt/(rt-1))))) 
  #Intensities
  FP = data.frame(matrix(data=NA, ncol=length(RT), nrow=length(dur)))
  colnames(FP) <- RT
  rownames(FP) <- dur
  #For 5 years
  FP[1,1] <- (Gumbel[2,1]+Gumbel[1,1]*K[1,1])/0.25
  FP[2,1] <- (Gumbel[2,2]+Gumbel[1,2]*K[1,1])/0.5
  FP[3,1] <- Gumbel[2,3]+Gumbel[1,3]*K[1,1]
  FP[4,1] <- (Gumbel[2,4]+Gumbel[1,4]*K[1,1])/3
  FP[5,1] <- (Gumbel[2,5]+Gumbel[1,5]*K[1,1])/8
  FP[6,1] <- (Gumbel[2,6]+Gumbel[1,6]*K[1,1])/12
  FP[7,1] <- (Gumbel[2,7]+Gumbel[1,7]*K[1,1])/24
  #For 10 years
  FP[1,2] <- (Gumbel[2,1]+Gumbel[1,1]*K[1,2])/0.25
  FP[2,2] <- (Gumbel[2,2]+Gumbel[1,2]*K[1,2])/0.5
  FP[3,2] <- Gumbel[2,3]+Gumbel[1,3]*K[1,2]
  FP[4,2] <- (Gumbel[2,4]+Gumbel[1,4]*K[1,2])/3
  FP[5,2] <- (Gumbel[2,5]+Gumbel[1,5]*K[1,2])/8
  FP[6,2] <- (Gumbel[2,6]+Gumbel[1,6]*K[1,2])/12
  FP[7,2] <- (Gumbel[2,7]+Gumbel[1,7]*K[1,2])/24
  #For 25 years
  FP[1,3] <- (Gumbel[2,1]+Gumbel[1,1]*K[1,3])/0.25
  FP[2,3] <- (Gumbel[2,2]+Gumbel[1,2]*K[1,3])/0.5
  FP[3,3] <- (Gumbel[2,3]+Gumbel[1,3]*K[1,3])
  FP[4,3] <- (Gumbel[2,4]+Gumbel[1,4]*K[1,3])/3
  FP[5,3] <- (Gumbel[2,5]+Gumbel[1,5]*K[1,3])/8
  FP[6,3] <- (Gumbel[2,6]+Gumbel[1,6]*K[1,3])/12
  FP[7,3] <- (Gumbel[2,7]+Gumbel[1,7]*K[1,3])/24
  #For 50 years
  FP[1,4] <- (Gumbel[2,1]+Gumbel[1,1]*K[1,4])/0.25
  FP[2,4] <- (Gumbel[2,2]+Gumbel[1,2]*K[1,4])/0.5
  FP[3,4] <- (Gumbel[2,3]+Gumbel[1,3]*K[1,4])
  FP[4,4] <- (Gumbel[2,4]+Gumbel[1,4]*K[1,4])/3
  FP[5,4] <- (Gumbel[2,5]+Gumbel[1,5]*K[1,4])/8
  FP[6,4] <- (Gumbel[2,6]+Gumbel[1,6]*K[1,4])/12
  FP[7,4] <- (Gumbel[2,7]+Gumbel[1,7]*K[1,4])/24
  #For 100 years
  FP[1,5] <- (Gumbel[2,1]+Gumbel[1,1]*K[1,5])/0.25
  FP[2,5] <- (Gumbel[2,2]+Gumbel[1,2]*K[1,5])/0.5
  FP[3,5] <- (Gumbel[2,3]+Gumbel[1,3]*K[1,5])
  FP[4,5] <- (Gumbel[2,4]+Gumbel[1,4]*K[1,5])/3
  FP[5,5] <- (Gumbel[2,5]+Gumbel[1,5]*K[1,5])/8
  FP[6,5] <- (Gumbel[2,6]+Gumbel[1,6]*K[1,5])/12
  FP[7,5] <- (Gumbel[2,7]+Gumbel[1,7]*K[1,5])/24
#Result summary
invisible(readline(prompt="Press [enter] to see results"))
message("================================================================================
Maximum Annual Precipitation")
invisible(print(Max_Pre))
message("================================================================================
Standard Deviation and Annual Average Precipitation")
invisible(print(Gumbel))
message("================================================================================
Frequency Precipitation values")
invisible(print(FP))
message("================================================================================")
invisible(readline(prompt="Press [enter] to plot IDF curves"))
#Plot
IDF <- ggplot(FP) +
  geom_line(aes(x=dura, y =FP$`5 years`),
            colour = 'red', size = 1) +
  geom_line(aes(x=dura, y =FP$`10 years`),
            colour = 'blue', size = 1) +
  geom_line(aes(x=dura, y =FP$`25 years`),
            colour = 'green', size = 1) +
  geom_line(aes(x=dura, y =FP$`50 years`),
            colour = 'black', size = 1) +
  geom_line(aes(x=dura, y =FP$`100 years`),
             colour = 'yellow', size = 1) + xlab("Duration (min)") + ylab("Intensity (mm/hr)") + ggtitle("Intensity-Duration-Frequency curve (IDF)") 
print(IDF)
message("================================================================================
END OF CODE")