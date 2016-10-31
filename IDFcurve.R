#load data
raindata = read.table("rain.txt", header = TRUE)

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
min_ts <- zoo(, seq(start(data1m), end(data1m), by = "1 min"))
merged_ts <- merge(data1m, min_ts)
merged_ts[is.na(merged_ts)] <- 0
rm(dates, raindata, TSrain, data1m, min_ts)

# create intensity for durations
data_30min <- (rollapply(merged_ts, 30, sum))/0.5
min_30 <- data_30min[30:length(data_30min)]
data_1h <- rollapply(merged_ts, 60, sum)
hour_1 <- data_1h[60:length(data_1h)]
data_3h <- (rollapply(merged_ts, 180, sum))/3
hour_3 <- data_3h[180:length(data_3h)]
data_6h <- (rollapply(merged_ts, 360, sum))/6
hour_6 <- data_6h[360:length(data_6h)]
data_12h <- (rollapply(merged_ts, 720, sum))/12
hour_12 <- data_12h[720:length(data_12h)]
data_18h <- (rollapply(merged_ts, 1080, sum))/18
hour_18 <- data_12h[1080:length(data_12h)]
data_24h <- (rollapply(merged_ts, 1440, sum))/24
hour_24 <- data_24h[1440:length(data_24h)]
rm(merged_ts, data_30min, data_1h, data_3h, data_6h, data_12h, data_18h, data_24h)

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
rm(min_30, hour_1, hour_3, hour_6, hour_12, hour_18, hour_24)

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