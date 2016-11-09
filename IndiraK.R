library(zoo)
library(xts)
library(lubridate)
library(highfrequency)
library(matrixStats)
library(ggplot2)
#set working directory
mydata = read.table("17.txt", row.names=NULL)
ls()
names(mydata)
mydata$Datetime
#change the column names
colnames(mydata) [1]="Date"
colnames(mydata) [2]="Time"
#to put dates in TS into dates format
mydata$DateTime <- as.POSIXct(paste(mydata$Date, mydata$Time),format="%m/%d/%Y %H:%M",tz= "Europe/paris")
mydata$Date <- NULL
mydata$Time <- NULL
mydata$Depth <- 0.1
#xts = extensible time series to create new objectives
myTS <- xts(mydata[,-1], order.by=mydata[,1])
#aggregate is used to compute mean of all values with the same time series index
myts_gr <- aggregate(myTS, identity, sum)
#to check if memory is enough, if not then set the new limit, remove from Environment tab myts_omin
memory.limit()
memory.limit(size = 7000)
rm(myts_omin)
#creating sequence of myts_omin
myts_omin <- zoo(, seq(start(myts_gr), end(myts_gr), by = "1 min"))
myts_omin <- xts(myts_omin)
myts_merged <- merge(myts_gr, myts_omin)
#NA values to turn into 0
myts_merged[is.na(myts_merged)] <- 0
#create a data framework with main values of merged information
myts.df <- data.frame(date=index(myts_merged), coredata(myts_merged))
#same for 15 min subset
min15 <- myts.df[1]
#apply rolling function
min15$Depth <- rollapply(myts.df[2], 15, sum, fill = NA, align='left')
min15$Depth <- min15$Depth / 15
min15 <- aggregate( min15$Depth ~ year(min15$date), FUN = max)
names(min15) = c("date", "Depth")
MaxIntensity <- min15$date
#to get 1 h max intencity of 15 min for every year
MaxIntensity <- as.data.frame(MaxIntensity)
MaxIntensity$min15 <- min15$Depth * 60
names(MaxIntensity) = c("Year", "15min")
#do same for 1h, 2h, 6, 12, 24 h
min60 <- myts.df[1]
min60$Depth <- rollapply(myts.df[2], 60, sum, fill = NA, align='left')
min60$Depth <- min60$Depth / 60
min60 <- aggregate( min60$Depth ~ year(min60$date), FUN = max)
names(min60) = c("date", "Depth")
MaxIntensity$min60 <- min60$Depth * 60
min360 <- myts.df[1]
min360$Depth <- rollapply(myts.df[2], 360, sum, fill = NA, align='left')
min360$Depth <- min360$Depth / 360
min360 <- aggregate( min360$Depth ~ year(min360$date), FUN = max)
names(min360) = c("date", "Depth")
MaxIntensity$min360 <- min360$Depth * 60
min720 <- myts.df[1]
min720$Depth <- rollapply(myts.df[2], 720, sum, fill = NA, align='left')
min720$Depth <- min720$Depth / 720
min720 <- aggregate( min720$Depth ~ year(min720$date), FUN = max)
names(min720) = c("date", "Depth")
MaxIntensity$min720 <- min720$Depth * 60
min1440 <- myts.df[1]
min1440$Depth <- rollapply(myts.df[2], 1440, sum, fill = NA, align='left')
min1440$Depth <- min1440$Depth / 1440
min1440 <- aggregate( min1440$Depth ~ year(min1440$date), FUN = max)
names(min1440) = c("date", "Depth")
MaxIntensity$min1440 <- min1440$Depth * 60
#Make new file statdata per column you want the mean -1 first row remove the NA = false because not necessary and the dimension is always 1
#the statedata is to create the mean value first for each duration
statdata <- colMeans(MaxIntensity[,-1], na.rm = FALSE, dims = 1)
statdata <- data.frame(t(statdata))
names(statdata) = c("min15", "min60", "min360")
#standard deviation of the mean values on each column
statdata[2,1] <- sd(MaxIntensity$`15min`)
statdata[2,4] <- sd(MaxIntensity$min60)
statdata[2,5] <- sd(MaxIntensity$min120)
statdata[2,6] <- sd(MaxIntensity$min360)
statdata[2,7] <- sd(MaxIntensity$min720)
statdata[2,8] <- sd(MaxIntensity$min1440)
#equation gumble inundation + duration * frequency
#Xi = mean + SD * Return period
#calculate the numbers of gumbel
gumbel <- matrix(c(2,-0.164,5,0.719,10,1.305,50,2.592,100,3.137), nrow = 2, ncol = 5)
options(digits=3)
gumbel <- as.data.frame(gumbel)
names(gumbel) = c("2", "5", "10", "50", "100")
#create new table with intensity and durations
intdur <- matrix(data = NA, nrow = 8, ncol = 5)
intdur[1,1] = statdata[1,1] + gumbel[2,1] * statdata[2,1]
intdur[2,1] = statdata[1,2] + gumbel[2,1] * statdata[2,2]
intdur[3,1] = statdata[1,3] + gumbel[2,1] * statdata[2,3]
intdur[4,1] = statdata[1,4] + gumbel[2,1] * statdata[2,4]
intdur[5,1] = statdata[1,5] + gumbel[2,1] * statdata[2,5]
intdur[6,1] = statdata[1,6] + gumbel[2,1] * statdata[2,6]
intdur[7,1] = statdata[1,7] + gumbel[2,1] * statdata[2,7]
intdur[8,1] = statdata[1,8] + gumbel[2,1] * statdata[2,8]
intdur[1,2] = statdata[1,1] + gumbel[2,2] * statdata[2,1]
intdur[2,2] = statdata[1,2] + gumbel[2,2] * statdata[2,2]
intdur[3,2] = statdata[1,3] + gumbel[2,2] * statdata[2,3]
intdur[4,2] = statdata[1,4] + gumbel[2,2] * statdata[2,4]
intdur[5,2] = statdata[1,5] + gumbel[2,2] * statdata[2,5]
intdur[6,2] = statdata[1,6] + gumbel[2,2] * statdata[2,6]
intdur[7,2] = statdata[1,7] + gumbel[2,2] * statdata[2,7]
intdur[8,2] = statdata[1,8] + gumbel[2,2] * statdata[2,8]
intdur[1,3] = statdata[1,1] + gumbel[2,3] * statdata[2,1]
intdur[2,3] = statdata[1,2] + gumbel[2,3] * statdata[2,2]
intdur[3,3] = statdata[1,3] + gumbel[2,3] * statdata[2,3]
intdur[4,3] = statdata[1,4] + gumbel[2,3] * statdata[2,4]
intdur[5,3] = statdata[1,5] + gumbel[2,3] * statdata[2,5]
intdur[6,3] = statdata[1,6] + gumbel[2,3] * statdata[2,6]
intdur[7,3] = statdata[1,7] + gumbel[2,3] * statdata[2,7]
intdur[8,3] = statdata[1,8] + gumbel[2,3] * statdata[2,8]
intdur[1,4] = statdata[1,1] + gumbel[2,4] * statdata[2,1]
intdur[2,4] = statdata[1,2] + gumbel[2,4] * statdata[2,2]
intdur[3,4] = statdata[1,3] + gumbel[2,4] * statdata[2,3]
intdur[4,4] = statdata[1,4] + gumbel[2,4] * statdata[2,4]
intdur[5,4] = statdata[1,5] + gumbel[2,4] * statdata[2,5]
intdur[6,4] = statdata[1,6] + gumbel[2,4] * statdata[2,6]
intdur[7,4] = statdata[1,7] + gumbel[2,4] * statdata[2,7]
intdur[8,4] = statdata[1,8] + gumbel[2,4] * statdata[2,8]
intdur[1,5] = statdata[1,1] + gumbel[2,5] * statdata[2,1]
intdur[2,5] = statdata[1,2] + gumbel[2,5] * statdata[2,2]
intdur[3,5] = statdata[1,3] + gumbel[2,5] * statdata[2,3]
intdur[4,5] = statdata[1,4] + gumbel[2,5] * statdata[2,4]
intdur[5,5] = statdata[1,5] + gumbel[2,5] * statdata[2,5]
intdur[6,5] = statdata[1,6] + gumbel[2,5] * statdata[2,6]
intdur[7,5] = statdata[1,7] + gumbel[2,5] * statdata[2,7]
intdur[8,5] = statdata[1,8] + gumbel[2,5] * statdata[2,8]
#put the internal duration in a new data frama to see the table
intdur <- as.data.frame(intdur)
names(intdur) = c("2", "5", "10", "50", "100")
#plotting the new data frame
dataplot <- matrix(data = NA, nrow = 8, ncol = 1)
dataplot[1,1] = 0.25
dataplot[4,1] = 1
dataplot[5,1] = 2
dataplot[6,1] = 6
dataplot[7,1] = 12
dataplot[8,1] = 24
dataplot <- as.data.frame(dataplot)
dataplot$`2yr` <- intdur$`2`
dataplot$`5yr` <- intdur$`5`
dataplot$`10yr` <- intdur$`10`
dataplot$`50yr` <- intdur$`50`
dataplot$`100yr` <- intdur$`100`
colnames(dataplot)[1] = "dur"
#IDF graph
idfgraph <-
ggplot(dataplot, aes(dur)) +
geom_line(aes(y = `2yr`, colour = "2 years",), size = 1) +
geom_point(aes(y = `2yr`, colour = "2 years"), size = 1.5) +
geom_line(aes(y = `5yr`, colour = "5 years"), size = 1) +
geom_point(aes(y = `5yr`, colour = "5 years"), size = 1.5) +
geom_line(aes(y = `10yr`, colour = "10 years"), size = 1) +
geom_point(aes(y = `10yr`, colour = "10 years"), size = 1.5) +
geom_line(aes(y = `50yr`, colour = "50 years"), size = 1) +
geom_point(aes(y = `50yr`, colour = "50 years"), size = 1.5) +
geom_line(aes(y = `100yr`, colour = "100 years"), size = 1) +
geom_point(aes(y = `100yr`, colour = "100 years"), size = 1.5) +
ylab("Intensity, mm/hour") +
xlab("Duration, hours") +
labs(title = "Intensity-Duration-Frequency curves") +
theme(panel.background = element_rect(fill = "white smoke")) +
theme(panel.grid.major = element_line(colour = "black", linetype = "dotted")) +
theme(panel.grid.minor = element_line(colour = "white smoke", linetype = "dotted")) +
theme(legend.background = element_rect(fill="white", size=0.5,
linetype="solid", colour ="darkblue")) +
labs(color = "Return Periods, years:") +






