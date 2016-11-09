library(zoo) # to perform computation on irregular time series
library(ggplot2) # for creating elegant and complex plots
library(xts) #library for working with time series
library(TTR) # Functions and data to construct technical trading rules with R
library(plyr) # tool for splitting time series
library(tidyr) # used to tidy up data

# reading dataset 
dat <- read.table("7.txt", header = TRUE)
dat$DateTime<- as.POSIXct(paste(dat$date, dat$time), format="%m/%d/%Y %H:%M", tz = "UTC")
dat$date = NULL
dat$time = NULL
dat$Depth <- 0.1

# Forming a [zoo] data frame for the irregular time series
dagg<-aggregate(dat[,2], list(dat$DateTime), sum ) # aggregates the rainfall data by date and time
dat.zoo<-zoo(dagg[,2],dagg[,1]) # converting to zoo 
dmrg <- merge(dat.zoo, zoo(,seq(start(dat.zoo),end(dat.zoo),by="min")), all=TRUE)
dmrg[is.na(dmrg)] <- 0
dframe <- data.frame(dmrg)

# Making a 1 minute time series to create a regular time series
from <- as.POSIXct("1994-02-08 01:09:00", tz = "UTC")
to <- as.POSIXct("2016-08-30 14:30:00", tz = "UTC")
ts.1min <- data.frame(seq(from, to, by = ("min")))
colnames(ts.1min)[1] = "DateTime"
d1min <- cbind.data.frame(ts.1min[1], dframe[1])
colnames(d1min)[1]="DateTime"
colnames(d1min)[2]="Depth"
#plot(d1min, ylab="Depth", xlab="Year")


# Aggregating the 1 minute time series to 1hr, 2hr, 4hr, 6hr, 12hr and 24hr
d1zoo <- zoo(d1min$Depth, order.by = d1min$DateTime)
d60s <- rollapply(d1zoo, 60, sum, fill=NA, align='right')
d2hs <- rollapply(d1zoo, 120, sum, fill=NA, align='right')
d4hs <- rollapply(d1zoo, 240, sum, fill=NA, align='right')
d6hs <- rollapply(d1zoo, 360, sum, fill=NA, align='right')
d12s <- rollapply(d1zoo, 720, sum, fill=NA, align='right')
d24s <- rollapply(d1zoo, 1440, sum, fill=NA, align='right')

# Creating a data frame for further computation
d1h <- data.frame(d60s)
d2h <- data.frame(d2hs/2)
d4h <- data.frame(d4hs/4)
d6h <- data.frame(d6hs/6)
d12h <- data.frame(d12s/12)
d24h <- data.frame(d24s/24)

# Combining timestamp with the data frames for further computation
df1 <- cbind.data.frame(ts.1min[1], d1h[1])
df2 <- cbind.data.frame(ts.1min[1], d2h[1])
df4 <- cbind.data.frame(ts.1min[1], d4h[1])
df6 <- cbind.data.frame(ts.1min[1], d6h[1])
df12 <- cbind.data.frame(ts.1min[1], d12h[1])
df24 <- cbind.data.frame(ts.1min[1], d24h[1])

# Substituting NA values with 0
df1[is.na(df1)] <- 0
df2[is.na(df2)] <- 0
df4[is.na(df4)] <- 0
df6[is.na(df6)] <- 0
df12[is.na(df12)] <- 0
df24[is.na(df24)] <- 0


colnames(df1)[1] = "DateTime"

# Splitting the DateTime format to Year for further computation
df1$year <- strftime(df1$DateTime, "%Y")
df2$year <- strftime(df2$DateTime, "%Y")
df4$year <- strftime(df4$DateTime, "%Y")
df6$year <- strftime(df6$DateTime, "%Y")
df12$year <- strftime(df12$DateTime, "%Y")
df24$year <- strftime(df24$DateTime, "%Y")

# Aggregating 1hr, 2hr, 4hr, 6hr, 12hr, and 24hr maximum values per year
d1mx <- aggregate(df1[, 2], list(df1$year), max)
d2mx <- aggregate(df2[, 2], list(df2$year), max)
d4mx <- aggregate(df4[, 2], list(df4$year), max)
d6mx <- aggregate(df6[, 2], list(df6$year), max)
d12mx <- aggregate(df12[, 2], list(df12$year), max)
d24mx <- aggregate(df24[, 2], list(df24$year), max)

# Computation of Rainfall Intensity
Tyears <- c(2,5,10,50,100) # to create a matrix of 2, 5, 10, 50, 100 years
duration <- c(1,2, 4, 6, 12, 24) # to create a matrix of 1hr, 2hr, 4hr, 6hr, 12hr, and 24hr
Kt <- c(-0.164, 0.719, 1.305, 2.592, 3.137) # to introduce Gumbel distribution constants 
Gumb <- data.frame(Tyears, Kt) 
mean1 <- c(mean(d1mx$x),mean(d2mx$x),mean(d4mx$x),mean(d6mx$x),mean(d12mx$x),mean(d24mx$x)) # to compute mean of the maximum aggregates 
sd1 <- c(sd(d1mx$x),sd(d2mx$x),sd(d4mx$x),sd(d6mx$x),sd(d12mx$x),sd(d24mx$x)) # to compute standard deviation of the maximum aggregates
Int.tab <- data.frame(duration, mean1, sd1)

#Calculation of Rainfall Intensity using Gumbel Distribution 
Int.tab$Yr2 <- mean1 + sd1*Kt[1]
Int.tab$Yr5 <- mean1 + sd1*Kt[2]
Int.tab$Yr10 <- mean1 + sd1*Kt[3]
Int.tab$Yr50 <- mean1 + sd1*Kt[4]
Int.tab$Yr100 <- mean1 + sd1*Kt[5]

# Combining different duraration intensities on a single data frame 
d_comp <- data.frame(d1mx$Group.1, d1mx$x, d2mx$x, d4mx$x, d6mx$x, d12mx$x, d24mx$x)
colnames(d_comp) = c("Year", "1hr", "2hr", "4hr", "6hr", "12hr", "24hr")


# Drawing IDF curves
lables <- c("Yr2", "Yr5", "Yr10", "Yr50", "Yr100") 
plot(Int.tab$duration, Int.tab$Yr100, xlab="Duration",ylab="Intensity (mm/hr)", type = "b", pch=2, lty = 1, lwd=1,  col = "blue")
legend("topright", title = "Legend", lables, lwd = 1, cex = 0.75, col = c("454","139","575","133","blue"))
lines(Int.tab$duration, Int.tab$Yr50, type = "b", pch=2, lty = 1, lwd=1,  col = "133")
lines(Int.tab$duration, Int.tab$Yr10, type = "b", pch=2, lty = 1, lwd=1,  col = "575")
lines(Int.tab$duration, Int.tab$Yr5, type = "b", pch=2, lty = 1, lwd=1,  col = "139")
lines(Int.tab$duration, Int.tab$Yr2, type = "b", pch=2, lty = 1, lwd=1,  col = "454")