#Time Series Assignment for dataset "5.txt"
#by: Ahmed Nasr
#Set working directory
setwd("~/")
setwd("E:/EMFRM5/UPC/Courses/Debris Flow and Flash Floods/TimeSeries/Oct24/Nov9")

#Block1
#Required Libraries/Packages
require(chron)
require(forecast)
require(ggplot2)
require(highfrequency)
require(hydroTSM)
require(matrixStats)
require(stats)
require(timeDate)
require(xts)
require(tseries)
require(TTR)
require(zoo)
require(dplyr)
require(tidyr)
require(RcppRoll)
require(gumbel)
require(lubridate)

#Block2
#read the two columns of date and time separately
mydata=read.table("5.txt", header=TRUE)
mydata$date=as.character(mydata$date)
mydata$time=as.character(mydata$time)

#convert date and time into datetime to be understood by R
datetime<-as.POSIXct(paste(mydata$date, mydata$time), format="%m/%d/%Y %H:%M", tz = "")
years <- year(datetime)
mydata <- data.frame(mydata, datetime,years, depth=0.1)

#Delete date and time columns
mydata$date <- NULL
mydata$time <- NULL

#Convert dataframe "mydata" into a uniform time series of 1 minute interval
dat.xts <- xts(x = mydata$depth, mydata$datetime)
colnames(dat.xts)[1] <- "depth"

ts_aggr <- aggregate(dat.xts, identity, sum)
ts_full <- zoo(, seq(start(ts_aggr), end(ts_aggr), by = "1 min"))
tsfinal <- merge(ts_aggr, ts_full)
tsfinal[is.na(tsfinal)] <- 0
tsfinal <- data.frame(date=(index(tsfinal)), depth=coredata(tsfinal))
colnames(tsfinal)[2]<-"depth"


#Block3
#Obtain max Intensity per year for different durations in y in min
y <- c(5,15,30,60,180,360,720,1440,2880,4320)
for (i in 1:length(y))
{
  #i=1
  ts<-data.frame(date=tsfinal$date ,Int=60*ma(tsfinal$depth, order = y[i], centre=TRUE))
  colnames(ts)[2]<-"Int"
  ts[is.na(ts)] <- 0
  ts1 <- aggregate( ts$Int ~ year(ts$date), FUN = max)
  colnames(ts1)<-c("year","Int")
  if (i < 2){
    Intensitymax <- ts1
  }
  else{
    ts_DF1 <- ts1
    Intensitymax[1+i]<-ts_DF1[2] 
  }
  rm(ts)
  rm(ts1)
}
colnames(Intensitymax)<- c("Year", "5min", "15min", "30min" , "1hr", "3hr", "6hr", "12hr", "24hr", "48hr", "72hr")
Intensitymax$Year <-NULL
#to remove years(any year) with no data (0 values) in case you have gaps of a year or more in your data
Intensitymax<-Intensitymax[!rowSums(Intensitymax[-c(1:2)] == 0) >= 1,]

#Block4
#statistics for Gumbel dist (mean, stdev, alpha, and beta)
Intensitymax1<-as.matrix.data.frame(Intensitymax)
mean<- colMeans(Intensitymax)
stdev<- colSds(Intensitymax1)
mean<-data.frame(mean)
stdev<-data.frame(stdev)
alpha<-stdev*(sqrt(6))/pi
beta<-mean-0.5772*alpha
statistics1<-t(data.frame(mean=mean, stdev=stdev, alpha=alpha, beta=beta))
rownames(statistics1)<- c("mean", "stdev", "alpha", "beta")
rm(Intensitymax1)

#Return Periods vector (in years)
x <- c(2,5,10,25,50,100,200,500)
gmblk1<-x
#final table with intensities and durations
intdurfreqx <- matrix(data = NA, nrow = length(y), ncol = length(x))
rownames(intdurfreqx) <- as.character(y)
colnames(intdurfreqx) <- as.character(t(x))
for (i in 1:length(y))#loop on rows (Durations)
{
  for (j in 1:length(x))#loop on columns(ReturnPeriods)
  {
    gmblk1[j] <-(-sqrt(6) / pi) * ((0.5772 + log((log(x[j] / (x[j] - 1))), base = exp(1))))
    intdurfreqx[i,j]= statistics1[1,i]+gmblk1[j]*statistics1[2,i]
  }
}

intdurfreq1x <- data.frame(duration=y/1,intdurfreqx)
colnames(intdurfreq1x) <- c("duration", "2yr", "5yr", "10yr", "25yr", "50yr", "100yr", "200yr", "500yr")
rownames(intdurfreq1x) <- c("5min", "15min", "30min" , "1hr", "3hr", "6hr", "12hr", "24hr", "48hr", "72hr")

#Block5
#Goodness of Fit using K-S Test
KS_Critical<-matrix(c(.995,.929,.829,.734,.669,.617,.576,.542,.513,.489,+
                      .468,.450,.432,.418,.404,.392,.381,.371,.361,.352,+
                      .344,.337,.330,.323,.317,.311,.305,.300,.295,.290,+
                                                .285,.281,.277,.273,.269))
K_out<-matrix(y)
K_calc<-matrix(y)

for (i in 1:length(y))
{
  Intensitymaxgf<-Intensitymax[i]
  colnames(Intensitymaxgf)[1]<-"Int"
  Intensitymaxgf<-data.frame(Int=Intensitymaxgf[with(Intensitymaxgf, order(Int)), ])
  I2<-matrix(t(Intensitymaxgf))
  for (j in 1:length(I2))
  {
    I2[j]=((j-0.44)/(length(I2)+0.12))
  }
  I3<-matrix(t(Intensitymaxgf))
  I3<-exp(-exp(-(I3-statistics1[4,i])/statistics1[3,i]))
  I4<-abs(I3-I2)
  I5<-I4
  I5[1]<-0
  for (j in 2:length(I5))
  {
    I5[j]=abs(I3[j]-I2[j-1])
  }
  K_calc[i]<-max(max(I4),max(I5))
  if (K_calc[i]<KS_Critical[length(I5)]){
    K_out[i]<-"K-S Test Passed"
  }else
  {
    K_out[i]<-"K-S Test Failed"
  }
}
View(K_out)

#Block6
#IDF CURVE PLOTTING
idfcurve<-
  ggplot(intdurfreq1x, aes(duration)) + 
  geom_line(aes(y = `2yr`, colour = "2 years",), size = 1) + 
  geom_point(aes(y = `2yr`, colour = "2 years"), size = 1.5) + 
  geom_line(aes(y = `5yr`, colour = "5 years"), size = 1) + 
  geom_point(aes(y = `5yr`, colour = "5 years"), size = 1.5) + 
  geom_line(aes(y = `10yr`, colour = "10 years"), size = 1) + 
  geom_point(aes(y = `10yr`, colour = "10 years"), size = 1.5) + 
  geom_line(aes(y = `25yr`, colour = "25 years"), size = 1) + 
  geom_point(aes(y = `25yr`, colour = "25 years"), size = 1.5) + 
  geom_line(aes(y = `50yr`, colour = "50 years"), size = 1) + 
  geom_point(aes(y = `50yr`, colour = "50 years"), size = 1.5) +
  geom_line(aes(y = `100yr`, colour = "100 years"), size = 1) + 
  geom_point(aes(y = `100yr`, colour = "100 years"), size = 1.5) + 
  geom_line(aes(y = `200yr`, colour = "200 years"), size = 1) + 
  geom_point(aes(y = `200yr`, colour = "200 years"), size = 1.5) + 
  geom_line(aes(y = `500yr`, colour = "500 years"), size = 1) + 
  geom_point(aes(y = `500yr`, colour = "500 years"), size = 1.5) +
  ylab("Intensity, mm/hour") +
  xlab("Duration, minutes") + 
  labs(title = "IDF Curves") +
  theme(panel.background = element_rect(fill = "white smoke")) +
  theme(panel.grid.major = element_line(colour = "black", linetype = "dotted")) +
  theme(panel.grid.minor = element_line(colour = "white smoke", linetype = "dotted")) +
  theme(legend.background = element_rect(fill="white", size=0.5, 
                                         linetype="solid", colour ="darkblue")) +
  labs(color = "Return Periods, years:") + 
  theme(legend.position = c(.92, .8))

plot(idfcurve)

write.table(intdurfreq1x,"intdurfreq1x.txt", append = FALSE, sep="\t", row.names = FALSE)

#Block7
#synthetic precipitation generation based on one of the calculated return period
#and storm duration(multiple of 5 minutes) using method of Alternating Blocks
#A
answercalc <- readline(prompt="Do you want to calculate intensity for different return period and duration? (type Y for yes and N for n) :")
#B
if (answercalc == "Y"){
  rt <- readline(prompt="Select one fo the following return periods: type: 2 for 2yr, 3 for 5yr, 4 for 10yr, 5 for 25yr, 6 for 50yr, 7 for 100yr, 8 for 200yr, 9 for 500yr")
  durr <- readline(prompt="Enter duration in minutes as multiple of 5 minutes : ")
  rt1 <- as.numeric(rt)
  dur <- as.numeric(durr)
  duration<-data.frame(duration=(intdurfreq1x[1]))
  Int<- data.frame(Int=(intdurfreq1x[rt1]))
  Synth<-data.frame(duration,Int=Int)
  colnames(Synth)[2]<- "Int"
  #fitting IDF curve to a power law in teh following step
  idf=nls(Int~c/(f+duration^e),start=list(c=1000,f=10, e=1),data=Synth)
  duration<-data.frame(duration=seq(from = 5, to = dur, by = 5))
  Int1<-data.frame(Int1=predict(idf,newdata=duration))
  Pcum<-data.frame(Pcum=Int1*duration)
  Pcum1<-as.matrix(Pcum)
  Pdis<- Pcum1
  Pdis[2:length(Pcum1)]<-diff(Pcum1)
  Palt<-Pdis*0
  
  if (length(Palt)%%2==0){ 
    ctr=0
    for (j in 1:(length(Palt)/2))
    {
      Palt[1-j+(length(Palt)/2)]<-Pdis[j+ctr]
      Palt[j+(length(Palt)/2)]<-Pdis[j+ctr+1]
      ctr=ctr+1
    }
    
  }else
  { 
    Palt[((length(Palt)+1)/2)]<-Pdis[1]
    ctr=1
    for (j in 1:((length(Palt)-1)/2))
    {
      Palt[j+((length(Palt)+1)/2)]<-Pdis[j+ctr]
      Palt[1-j+((length(Palt)-1)/2)]<-Pdis[j+ctr+1]
      ctr=ctr+1
    }
  }
  Psynth <- data.frame(duration,Int1,Pcum,Pdis,Palt, row.names = NULL)
  colnames(Psynth) <- c("duration", "Intensity", "Pcum", "P", "P_alternating")
  plot(Psynth$duration,Psynth$P_alternating, xlab="Duration (mins)", ylab="Precipitation depth(mm)", type="S",main="Synthetic Precipitation")
  write.table(Psynth, "Psynth.txt", append = FALSE, sep="\t", row.names = FALSE)
  print("Program finished executing, please re-run the program")
 }else
{
  print("Program finished executing, please re-run the program")  
}

 