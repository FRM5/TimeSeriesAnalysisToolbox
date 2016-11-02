# library(lubridate)
# library(forecast)
# library(xts)
# library(googleVis)

#load data
raindata = read.table("2004.txt", header = TRUE)

#rename two columns
colnames(raindata)[1] = "dts"
colnames(raindata)[2] = "tms"

# convert to date and time
dates <-strptime(paste(raindata$dts, raindata$tms), "%m/%d/%Y %H:%M:%S", tz = "UTC")
raindata$dts <- dates

# add 0.1mm precipitation
raindata["precip"] <-  0.1
raindata$tms <- NULL

#create time series with zoo
TSrain <- xts(raindata[,-1], order.by=raindata[,1])
data1m <- period.apply(TSrain, endpoints(TSrain, 'minutes', 1), sum)
min_ts <- merge(data1m, zoo(, seq(start(data1m), end(data1m), by = "1 min")))
min_ts[is.na(min_ts)] <- 0
rm(dates, raindata, TSrain, data1m)

# create intensity for durations
data_30min <- rollapply(min_ts, 30, sum, align = 'left')/0.5
min_30 <- tapply(data_30min$data1m, year(index(data_30min)), max)
rm(data_30min)
data_1h <- rollapply(min_ts, 60, sum, align = 'left')
hour_1 <- tapply(data_1h$data1m, year(index(data_1h)), max)
rm(data_1h)
data_3h <- na.omit(rollapply(min_ts, 180, sum, align = 'left'))/3
hour_3 <- tapply(data_3h$data1m, year(index(data_3h)), max)
rm(data_3h)
data_6h <- na.omit(rollapply(min_ts, 360, sum, align = 'left'))/6
hour_6 <- tapply(data_6h$data1m, year(index(data_6h)), max)
rm(data_6h)
data_12h <- na.omit(rollapply(min_ts, 720, sum, align = 'left'))/12
hour_12 <- tapply(data_12h$data1m, year(index(data_12h)), max)
rm(data_12h)
data_18h <- na.omit(rollapply(min_ts, 1080, sum, align = 'left'))/18
hour_18 <- tapply(data_18h$data1m, year(index(data_18h)), max)
rm(data_18h)
data_24h <- na.omit(rollapply(min_ts, 1440, sum, align = 'left'))/24
hour_24 <- tapply(data_24h$data1m, year(index(data_24h)), max)
rm(min_ts, data_24h)

# create return periods
RP = c(5, 10, 25, 50, 100, 200)

# calculate rain intensity for each duration and return period 
RP_30min <- (-log(-log(1-(1/RP)))) * ((sd(min_30) * sqrt(6))/pi) + (mean(min_30)-0.5772 * ((sd(min_30) * sqrt(6))/pi))
RP_1h <- (-log(-log(1-(1/RP)))) * ((sd(hour_1) * sqrt(6))/pi) + (mean(hour_1)-0.5772 * ((sd(hour_1) * sqrt(6))/pi))
RP_3h <- (-log(-log(1-(1/RP)))) * ((sd(hour_3) * sqrt(6))/pi) + (mean(hour_3)-0.5772 * ((sd(hour_3) * sqrt(6))/pi))
RP_6h <- (-log(-log(1-(1/RP)))) * ((sd(hour_6) * sqrt(6))/pi) + (mean(hour_6)-0.5772 * ((sd(hour_6) * sqrt(6))/pi))
RP_12h <- (-log(-log(1-(1/RP)))) * ((sd(hour_12) * sqrt(6))/pi) + (mean(hour_12)-0.5772 * ((sd(hour_12) * sqrt(6))/pi))
RP_18h <- (-log(-log(1-(1/RP)))) * ((sd(hour_18) * sqrt(6))/pi) + (mean(hour_18)-0.5772 * ((sd(hour_18) * sqrt(6))/pi))
RP_24h <- (-log(-log(1-(1/RP)))) * ((sd(hour_24) * sqrt(6))/pi) + (mean(hour_24)-0.5772 * ((sd(hour_24) * sqrt(6))/pi))

# join results in one table
all_data <- data.frame(t(data.frame(RP_30min, RP_1h, RP_3h, RP_6h, RP_12h, RP_18h, RP_24h)))
rm(RP_30min, RP_1h, RP_3h, RP_6h, RP_12h, RP_18h, RP_24h)
all_data["duration"] <- c(0.5, 1, 3, 6, 12, 18, 24)
colnames(all_data)[1] = "5 year"
colnames(all_data)[2] = "10 year"
colnames(all_data)[3] = "25 year"
colnames(all_data)[4] = "50 year"
colnames(all_data)[5] = "100 year"
colnames(all_data)[6] = "200 year"

# plot results
plot(gvisLineChart(all_data, xvar = "duration", yvar = c("5 year", "10 year", "25 year", "50 year", "100 year", "200 year"), options = list(gvis.editor = "Edit chart", title = "IDF curve for February 1994 to August 2016 rainfall event", vAxes="[{title:'intensity (mm/hour)'}]", hAxes="[{title:'duration (hour)'}]", curveType = "function", width = 1300, height = 700)))