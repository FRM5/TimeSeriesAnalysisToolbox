setwd("~/Documents/EMFRM 2015-2017/UPC Barcelona/Drought - RStudio")
# Load libraries
library(xts)
library(lubridate)
library(highfrequency)
library(extRemes)
library(ggplot2)
library(zoo)
library(hydroTSM)

# Load data & Change column names
mydata = read.table("8.txt", header = TRUE)
colnames(mydata)[1] = "Date"
colnames(mydata)[2] = "Time"

# Change Date format
mydata$DateTime <- as.POSIXct(paste(mydata$Date, mydata$Time), format = "%m/%d/%Y %H:%M" , tz = "UTC" )

# Add column "Depth" with value 0.1
mydata$Depth <- 0.1

# Returning NULL value
mydata$Date <- NULL
mydata$Time <- NULL

# Create Time series with the ZOO package
TS <- xts(mydata[,-1], order.by = mydata[,1])
# create 1 min uniform time series
TS_aggr <- aggregate(TS, identity, sum)
TS_1min <- zoo(, seq(start(TS_aggr), end(TS_aggr), by = "1 min"))
TS_merg <- merge(TS_aggr, TS_1min)
TS_merg[is.na(TS_merg)] <- 0
# Create Data Frame
TS_DF <- data.frame(date = index(TS_merg), coredata(TS_merg))
colnames(TS_DF) <- c("Date_Time", "Depth")

# Creating Time series for different time intervalls
# 15 min
min15 <- TS_DF[1]
min15$Depth <- rollapply(TS_DF[2], width = 15, FUN = sum, fill = NA, align = 'center' )
min15 <- aggregate(min15$Depth ~ year(min15$Date_Time), FUN = max )
min15$Depth <- (((min15$Depth)/15)*60)

# 30 min
min30 <- TS_DF[1]
min30$Depth <- rollapply(TS_DF[2], width = 30, FUN = sum, fill = NA, align = 'center' )
min30 <- aggregate(min30$Depth ~ year(min30$Date_Time), FUN = max )
min30$Depth <- (((min30$Depth)/30)*60)

# 45 min
min45 <- TS_DF[1]
min45$Depth <- rollapply(TS_DF[2], width = 45, FUN = sum, fill = NA, align = 'center' )
min45 <- aggregate(min45$Depth ~ year(min45$Date_Time), FUN = max )
min45$Depth <- (((min45$Depth)/45)*60)

# 1 hour (60 min)
min60 <- TS_DF[1]
min60$Depth <- rollapply(TS_DF[2], width = 60, FUN = sum, fill = NA, align = 'center' )
min60 <- aggregate(min60$Depth ~ year(min60$Date_Time), FUN = max )
min60$Depth <- (((min60$Depth)/60)*60)

# 12 hour (720 min)
min720 <- TS_DF[1]
min720$Depth <- rollapply(TS_DF[2], width = 720, FUN = sum, fill = NA, align = 'center' )
min720 <- aggregate(min720$Depth ~ year(min720$Date_Time), FUN = max )
min720$Depth <- (((min720$Depth)/720)*60)

# 24 hour (1440 min)
min1440 <- TS_DF[1]
min1440$Depth <- rollapply(TS_DF[2], width = 1440, FUN = sum, fill = NA, align = 'center' )
min1440 <- aggregate(min1440$Depth ~ year(min1440$Date_Time), FUN = max )
min1440$Depth <- (((min1440$Depth)/1440)*60)

# Add Data Frame with max. annual intensity
MaxIntensity <- data.frame(min15$Depth, min30$Depth, min45$Depth, min60$Depth, min720$Depth, min1440$Depth)
names(MaxIntensity) <- c("15min", "30min", "45min", "1hr", "12hr", "24hr")
row.names(MaxIntensity) <- c("1994":"2016")

# Gumbel distribution and Annual Average
RT <- c(2,5,10,50,100) #Return Periods
DUR <- c("15min", "30min", "45min", "1hr", "12hr", "24hr")
Gumbel <- data.frame(matrix(data = NA, ncol = length(DUR), nrow = 2))
colnames(Gumbel) <- DUR
rownames(Gumbel) <- c("SD", "Annual Aver.")
Gumbel[1,1] <- sd(MaxIntensity$`15min`)
Gumbel[1,2] <- sd(MaxIntensity$`30min`)
Gumbel[1,3] <- sd(MaxIntensity$`45min`)
Gumbel[1,4] <- sd(MaxIntensity$`1hr`)
Gumbel[1,5] <- sd(MaxIntensity$`12hr`)
Gumbel[1,6] <- sd(MaxIntensity$`24hr`)
Gumbel[2,] <- colMeans(MaxIntensity, na.rm = FALSE, dims = 1)

# Rainfall Intensities per duration
# 15 min
Int15min <- data.frame(matrix(data = NA, ncol = 3, nrow = 5))
colnames(Int15min) <- c("K", "Freq. Prec.", "I")
rownames(Int15min) <- c("2y", "5y", "10y", "50y", "100y")
Int15min$K <- (-sqrt(6)/pi)*(0.5772 + log(log(RT/(RT-1))) )
Int15min$`Freq. Prec.` <- Gumbel[2,1]+(Int15min$K*Gumbel[1,1])
Int15min$I <- Int15min$`Freq. Prec.`/0.25

# 30 min
Int30min <- data.frame(matrix(data = NA, ncol = 3, nrow = 5))
colnames(Int30min) <- c("K", "Freq. Prec.", "I")
rownames(Int30min) <- c("2y", "5y", "10y", "50y", "100y")
Int30min$K <- (-sqrt(6)/pi)*(0.5772 + log(log(RT/(RT-1))) )
Int30min$`Freq. Prec.` <- Gumbel[2,2]+(Int30min$K*Gumbel[1,2])
Int30min$I <- Int30min$`Freq. Prec.`/0.5

# 45 min
Int45min <- data.frame(matrix(data = NA, ncol = 3, nrow = 5))
colnames(Int45min) <- c("K", "Freq. Prec.", "I")
rownames(Int45min) <- c("2y", "5y", "10y", "50y", "100y")
Int45min$K <- (-sqrt(6)/pi)*(0.5772 + log(log(RT/(RT-1))) )
Int45min$`Freq. Prec.` <- Gumbel[2,3]+(Int45min$K*Gumbel[1,3])
Int45min$I <- Int45min$`Freq. Prec.`/0.75

# 1 hour
Int1hr <- data.frame(matrix(data = NA, ncol = 3, nrow = 5))
colnames(Int1hr) <- c("K", "Freq. Prec.", "I")
rownames(Int1hr) <- c("2y", "5y", "10y", "50y", "100y")
Int1hr$K <- (-sqrt(6)/pi)*(0.5772 + log(log(RT/(RT-1))) )
Int1hr$`Freq. Prec.` <- Gumbel[2,4]+(Int1hr$K*Gumbel[1,4])
Int1hr$I <- Int1hr$`Freq. Prec.`/1

# 12 hour
Int12hr <- data.frame(matrix(data = NA, ncol = 3, nrow = 5))
colnames(Int12hr) <- c("K", "Freq. Prec.", "I")
rownames(Int12hr) <- c("2y", "5y", "10y", "50y", "100y")
Int12hr$K <- (-sqrt(6)/pi)*(0.5772 + log(log(RT/(RT-1))) )
Int12hr$`Freq. Prec.` <- Gumbel[2,5]+(Int12hr$K*Gumbel[1,5])
Int12hr$I <- Int12hr$`Freq. Prec.`/12

# 24 hour
Int24hr <- data.frame(matrix(data = NA, ncol = 3, nrow = 5))
colnames(Int24hr) <- c("K", "Freq. Prec.", "I")
rownames(Int24hr) <- c("2y", "5y", "10y", "50y", "100y")
Int24hr$K <- (-sqrt(6)/pi)*(0.5772 + log(log(RT/(RT-1))) )
Int24hr$`Freq. Prec.` <- Gumbel[2,6]+(Int24hr$K*Gumbel[1,6])
Int24hr$I <- Int24hr$`Freq. Prec.`/24

# Merge Intensities into IDF
I_total <- data.frame(matrix(data = NA, ncol = 5, nrow = 6))
colnames(I_total) <- c("2y", "5y", "10y", "50y", "100y")
rownames(I_total) <- DUR
I_total[1,] <- Int15min$I
I_total[2,] <- Int30min$I
I_total[3,] <- Int45min$I
I_total[4,] <- Int1hr$I
I_total[5,] <- Int12hr$I
I_total[6,] <- Int24hr$I

#print the IDF
dur <- c(15, 30, 45, 60, 720, 1440)
IDF <- ggplot(I_total) + 
  geom_line(aes(x = dur, y = I_total$`2y`), 
            colour = 'red', size = 1) + 
  geom_line(aes(x = dur, y = I_total$`5y`), 
            colour = 'blue', size = 1) + 
  geom_line(aes(x = dur, y = I_total$`10y`), 
            colour = 'yellow', size = 1) + 
  geom_line(aes(x = dur, y = I_total$`50y`), 
            colour = 'black', size = 1) + 
  geom_line(aes(x = dur, y = I_total$`100y`), 
            colour = 'orange', size = 1) +
  xlab("Duration (min)") + ylab("Intensity (mm/h)") + ggtitle("IDF Curves")

print(IDF)
