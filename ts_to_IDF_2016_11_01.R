#################   SETUP

# Set the working directory and require the packages
# Not sure all these packages are still required - can scrap the unused ones after all is sorted
# setwd("C:/DebrisFlowFlashFlood/timeseries_ex/data/15")
library(zoo)#for working with non-uniform ts
library(highfrequency) #for aggregating the non-uniform time series
library(xts) #for working with ts
library(lmomco) #for extreme value analysis
library(reshape2) # for melting the data frame before plotting with ggplot2
library(ggplot2) #for plotting the IDF
library(hydroTSM) #for the rolling mean/sum
library(data.table) # for conversion to uniform ts
library(chron) # yet another date/time package...
library(lubridate) # for conversion to uniform ts
library(QRM) #for gumbel and gof
library(goftest) # for goodness of fit

############################### CREATE SUPPORTING VARIABLES

# Need to edit this if your data has different years (should really make this flow from the read-in process)
yrlist <- c(2005:2016)
n = length(yrlist)
dur = c(0.25,0.5,1,2,3,4,5,6,12,24)
T <- c(2,5,10,20,25,50,100,250,500)

################################ READ DATA IN AND CONVERT TO TS

# Read the table (without header) and change column names
# FYI - in this dataset, each entry (each "tick") represents 0.1 mm of water (like a tipping bucket)
mydata = read.table("15nh.txt")
colnames(mydata)[1] = "Date"
colnames(mydata)[2] = "Time"

# Convert mydata to time series: concatenate date and time, then stores as a POSIXlt POSIXt object...in central european summer time (not going to worry about time zones for now)
ticks <- as.POSIXct(paste(mydata$Date, mydata$Time, sep = " "), format="%m/%d/%Y %H:%M:%S")

############################## CONVERT TO UNIFORM TS

ticksframe <- data.frame(timestamp = ticks, count = 0.1)
# A data.table is a special type of data frame. It does everything a dataframe does, plus more functions.
tickstable <- data.table(ticksframe)

# Group tickstable by the minute, and for each set of obs in the same minute, sum the 0.1mms.
minutely <- tickstable[, sum(count), by=floor_date(timestamp, 'minute')]
names(minutely) <- c('time', 'precip_mm')

# Create empty data frame of the proper length (based on ticks length - should fix to have this determined from the read-in process)
ts_empty <- trunc(seq(from = as.POSIXct("2005-06-21 16:24:34", format="%Y-%m-%d %H:%M:%S"), to = as.POSIXct("2016-08-30 05:24:39", format="%Y-%m-%d %H:%M:%S"), by=60),units=c("mins"))
ts_empty <- data.frame(time=ts_empty) # merge only works on dataframes.

# Merge the empty ts with the "minutely" data, based on common column "time," it's like a sql join
complete_ts <- merge(ts_empty, minutely, by='time', all.x=TRUE)
# * merge names its first two arguments 'x' and 'y'. So all.x says return a row ...
#      (possibly NA -- no data) for every time in the argument I gave ts_empty.
#      If there were a time in 'y' (minutely) that didn't match one in ts_empty, it
#      would be thrown away. Because there were a few NA times in the original dataset, these
#      are discarded and don't exist in the complete_ts (original sum was 5520.5, now it is 5519.5
#      b/c there were 10 NAs).

# Now complete_ts has NA for long stretches of no precip. Set to 0.
complete_ts$precip_mm[is.na(complete_ts$precip_mm)] <- 0

######################## CALCULATE MAXIMUM ANNUAL INTENSITY BY DURATION

# Convert to zoo, perform moving window sum for each duration (align center default)
myTS <-zoo(complete_ts[,-1], order.by=complete_ts[,1])
msum0015 <- ma(myTS, win.len = 15, FUN = sum)
msum0030 <- ma(myTS, win.len = 30, FUN = sum)
msum0100 <- ma(myTS, win.len = 60, FUN = sum)
msum0200 <- ma(myTS, win.len = 120, FUN = sum)
msum0300 <- ma(myTS, win.len = 180, FUN = sum)
msum0400 <- ma(myTS, win.len = 240, FUN = sum)
msum0500 <- ma(myTS, win.len = 300, FUN = sum)
msum0600 <- ma(myTS, win.len = 360, FUN = sum)
msum1200 <- ma(myTS, win.len = 720, FUN = sum)
msum2400 <- ma(myTS, win.len = 1440, FUN = sum)

# This is the table with max sum precip for each year (rows) and duration (columns)
annomax <- data.frame(year=yrlist, max0015 = apply.yearly(msum0015,FUN="max")[,1], max0030 = apply.yearly(msum0030,FUN="max")[,1], max0100 = apply.yearly(msum0100,FUN="max")[,1], max0200 = apply.yearly(msum0200,FUN="max")[,1], max0300 = apply.yearly(msum0300,FUN="max")[,1], max0400 = apply.yearly(msum0400,FUN="max")[,1], max0500 = apply.yearly(msum0500,FUN="max")[,1], max0600 = apply.yearly(msum0600,FUN="max")[,1], max1200 = apply.yearly(msum1200,FUN="max")[,1], max2400 = apply.yearly(msum2400,FUN="max")[,1])
# This is converted to hourly maximums, mm per hour
annomax_mmph <- data.frame(year=yrlist,max0015 = annomax$max0015*4, max0030 = annomax$max0030*2, max0100 = annomax$max0100, max0200 = annomax$max0200/2, max0300 = annomax$max0300/3, max0400 = annomax$max0400/4, max0500 = annomax$max0500/5, max0600 = annomax$max0600/6, max1200 = annomax$max1200/12, max2400 = annomax$max2400/24)
annomax_mmph$year = NULL
# This is a transposition of the annomax_mmph df used to plot against the final IDF curve
annomax_bydur = t(annomax_mmph)

#################### PERFORM GOODNESS OF FIT TESTS

# Note - goodness of fit tests should usually be performed on samples greater than what we have available (12 years per duration)
  # This could be helped by using the POT method, but it is not currently included in this exercise.

# Goodness of fit tests using Kolmogorov-Smirnov
  # Note - also attempted chi square test, but the answers appear to be unreliable, likely because sample size is so small (?)
  # Note - depending on the tool, different parameters have different names. For Gumbel, xi (other places known as mu) is location, alpha (other places known as beta or sigma) is for scale

  # Gumbel
    # Estimate parameters for each duration
sdgof <- c(0)
for (i in 1:length(dur)) sdgof[i] = sd(annomax_mmph[,i])
Pavegof <- c(0)
for (i in 1:length(dur)) Pavegof[i] = sum(annomax_mmph[i]/n)
beta <- c(0)
for (i in 1:length(dur)) beta[i] = sdgof[i]*sqrt(6)/pi
mu <- c(0)
for (i in 1:length(dur)) mu[i] = Pavegof[i] - 0.5772*beta[i]
    # Generate random values from the gumbel distribution of those parameters for each duration
ygum <- data.frame(matrix(ncol = length(dur), nrow = n))
for (i in 1:length(dur)) ygum[,i] = rGumbel(n,mu=mu[i],sigma=beta[i])
    # Test the randomly generated numbers against observed for each duration using Kolmogorov-Smirnov (and attempted chi square)
      # Note - tried to make a table of this, not working -- gumres_ks = data.frame(matrix(ncol = length(dur), nrow = n))
for (i in 1:length(dur)) print(ks.test(annomax_mmph[i],ygum[,i]))
    # I think chi square (below) is not working - they all come out the same, and the same value for p with log pearson...
      # for (i in 1:length(dur)) print(chisq.test(annomax_mmph[i],ygum[,i]))

  # Log-Pearson III
    # Estimate parameters for LP3 using lmomco2 functions for each duration
lp3mu <- c(0)
for (i in 1:length(dur)) lp3mu[i] = FrequencyAnalysis(annomax_mmph[,i],"lp3", nep = nonexceeds())$distribution$parameters$para[1]
lp3sigma <- c(0)
for (i in 1:length(dur)) lp3sigma[i] = FrequencyAnalysis(annomax_mmph[,i],"lp3", nep = nonexceeds())$distribution$parameters$para[2]
lp3gamma <-c(0)
for (i in 1:length(dur)) lp3gamma[i] = FrequencyAnalysis(annomax_mmph[,i],"lp3", nep = nonexceeds())$distribution$parameters$para[3]
    # Create a log pearson 3 function using its relationship to gamma distribution, recommended here: http://stats.stackexchange.com/questions/55563/estimating-parameters-of-log-pearson-iii-distribution-in-r
rlogpearson <- function(m,a,b,c) return( exp(rgamma(m,shape=a,rate=b) - c) )
    # Generate random numbers from the LP3 distribution with the parameters for each duration
ylp3 <- data.frame(matrix(ncol = length(dur), nrow = n))
for (i in 1:length(dur)) ylp3[,i] = rlogpearson(n,lp3mu[i],lp3sigma[i],lp3gamma[i])
    # Test the randomly generated numbers against observed for each duration using Kolmogorov-Smirnov (and attempted chi square)
for (i in 1:length(dur)) print(ks.test(annomax_mmph[i],ylp3[,i]))
    # Chi square is still not working, gives the same answer for Gumbel and LP3 for all durations...
      # for (i in 1:length(dur)) print(chisq.test(annomax_mmph[i],ylp3[,i]))

################# FREQUENCY ANALYSIS USING GUMBEL

# Note - could use the values provided by lmomco2 above, but I think Vicente said to run it ourselves?

# Calculate the standard dev for each duration
sd <- c(0)
for (i in 1:length(dur)) sd[i+1] <- sd(annomax_mmph[,i])
# Calculate the K for each return period
K <- c(0)
for (i in 1:length(T)) K[i] <- ((0-sqrt(6))/pi)*(0.5772+log(log(T[i]/(T[i]-1))))
# Create a table for K values
IDF_K <- data.frame(ReturnPd=T, mn15 = c(0), mn30 = c(0), mn60=c(0),hr2=c(0),hr3=c(0),hr4=c(0),hr5=c(0),hr6=c(0),hr12=c(0),hr24=c(0))
for (i in 1:length(dur)) IDF_K[i+1]= K
# Create a table for sd values
IDF_S <- data.frame(ReturnPd=T, mn15 = c(0), mn30 = c(0), mn60=c(0),hr2=c(0),hr3=c(0),hr4=c(0),hr5=c(0),hr6=c(0),hr12=c(0),hr24=c(0))
for (i in 1:length(T)) IDF_S[i,]= sd
# Create a table for the average for each duration
Pave <- c(0)
for (i in 1:length(dur)) Pave[i+1] <- sum(annomax_mmph[i]/n)
IDF_Pave <- data.frame(ReturnPd=T, mn15 = c(0), mn30 = c(0), mn60=c(0),hr2=c(0),hr3=c(0),hr4=c(0),hr5=c(0),hr6=c(0),hr12=c(0),hr24=c(0))
for (i in 1:length(T)) IDF_Pave[i,]= Pave

# Combine using the formula Intensity (for each return pd)=Pave+K*S
IDF = IDF_Pave + IDF_K*IDF_S

# Transpose the IDF :/
IDF$ReturnPd = NULL
rownames(IDF) = T
IDF = t(IDF)

# Add duration column for plotting, and re-set the column names
IDF_toplot = data.frame(Duration=dur, IDF[,1],IDF[,2],IDF[,3],IDF[,4],IDF[,5],IDF[,6],IDF[,7],IDF[,8],IDF[,9])
for (i in 1:length(T)) colnames(IDF_toplot)[i+1] = T[i]

######################### PLOT IDF CURVE AND ADD OBSERVED DATA FOR VISUAL COMPARISON

# Plot the IDF_toplot
plot(IDF_toplot$Duration, IDF_toplot$`2`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="Duration", ylab="Intensity", pch=18, col=26)
par(new=TRUE)
plot(IDF_toplot$Duration, IDF_toplot$`5`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", axes=F, pch=18, col=28)
par(new=TRUE)
plot(IDF_toplot$Duration, IDF_toplot$`10`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", axes=F, pch=18, col=30)
par(new=TRUE)
plot(IDF_toplot$Duration, IDF_toplot$`20`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", axes=F, pch=18, col=32)
par(new=TRUE)
plot(IDF_toplot$Duration, IDF_toplot$`25`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", axes=F, pch=18, col=34)
par(new=TRUE)
plot(IDF_toplot$Duration, IDF_toplot$`50`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", axes=F, pch=18, col=36)
par(new=TRUE)
plot(IDF_toplot$Duration, IDF_toplot$`100`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", axes=F, pch=18, col=38)
par(new=TRUE)
plot(IDF_toplot$Duration, IDF_toplot$`250`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", axes=F, pch=18, col=81)
par(new=TRUE)
plot(IDF_toplot$Duration, IDF_toplot$`500`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", axes=F, pch=18, col=48)

# Flip annomax for visual comparison to observed data
annomax_toplot = data.frame(Duration=c(.25,.5,1,2,3,4,5,6,12,24), annomax_bydur[,1],annomax_bydur[,2],annomax_bydur[,3],annomax_bydur[,4],annomax_bydur[,5],annomax_bydur[,6],annomax_bydur[,7],annomax_bydur[,8],annomax_bydur[,9],annomax_bydur[,10],annomax_bydur[,11],annomax_bydur[,12])
years = c(2005:2016)
for (i in 1:12) colnames(annomax_toplot)[i+1] = years[i]

# Add observed data to the plot (annual maximum intensities)
# I'm commenting this out for now, so the IDF curve comes out clean the first time

#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2005`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2006`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2007`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2008`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2009`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2010`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2011`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2012`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2013`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2014`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2015`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")
#par(new=TRUE)
#plot(annomax_toplot$Duration, annomax_toplot$`2016`, type='p', xlim=c(0,25), ylim=c(0,170), xlab="", ylab="", pch=1, col="black")

###################### GENERATE DESIGN STORM HYETOGRAPH (Alternating Blocks)

# This example is for 500 yr storm, 1hr intervals, 6-hour storm

# Extract durations from 1-6 hours for 500 yr flood
Dstorm_i <- IDF[,9]
Dstorm_i = Dstorm_i[3:8]

# Generate cumulative and incremental precip
cumulative = Dstorm_i * c(1:6)
incremental = c(cumulative[1],0)
for (i in 1:5) incremental[i+1] = cumulative[i+1]-sum(incremental[1:i])
# not finished - need to find a way to alternate the incremental series...
# also, is it possible that incremental hr 4 is higher than hr 3??? how?

# Plot reference https://www.r-bloggers.com/mastering-r-plot-part-1-colors-legends-and-lines/
# Plot reference http://www.sixhat.net/plotting-multiple-data-series-in-r.html

###############NEXT STEPS
# Use S-shape method or alternating blocks to get a sample precip event for one of the durations
# Things to clean up: dates in annomax, all the table flipping with IDFs..., some of the names are unclear..., unnecessary objects, automate things like the beginning/end dates (based on the original data)
                                                              