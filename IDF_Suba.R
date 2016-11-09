#Suba Madurangani Hewavidana
#Assignement-Flash Flood
#Flood Risk Management

#Calling packages
require(xts)
require(zoo)
library(lubridate)
library(timeDate)
library(extRemes)
library(ggplot2)

#Calling data
mydata=read.csv("13.csv",header=TRUE)

#Making data into date-time format
mydata$datetime <-as.POSIXct(mydata$datetime, format = "%m/%d/%Y %H:%M", usetz = FALSE)
mydata$rocord<- 0.1 #Assigning depth

a <- mydata$datetime[1] #calling first date in the series
b <- mydata$datetime[nrow(mydata)] #Calling last date in the series

#Summing the records per one minute
sum <- tapply(mydata$rocord, as.POSIXct(mydata$datetime, "minute"), sum)
ts_sum <- data.frame(sum)

# making the series in one minute intervall for the time span
Emty <- seq(a,b,"min") 
zero <- xts(rep(0,length(Emty)), order.by=Emty)

#Making a time series with available record
Merge <-merge(x = zero, y = ts_sum,  all.x =   TRUE, fill=0)[,2]

rm(zero,mydata,ts_sum)#remove unwanted data to make the space for the next step

#Moving Average to get the maximum intensity
ctr=0
dur <- c(5,10,30,60,120,360,720,1440)#in minute
time <- dur/60#in hour

for(i in dur){
  mam <- data.frame(date=index(Merge))
  
  mam$depth <- rollapply(coredata(Merge),i,sum,fill=0)
  mam$depth <- mam$depth/time[ctr+1]
  
  if (ctr < 1)#extracting maximums
  {
    maxd <- aggregate(mam$depth ~ year(mam$date),mam, FUN=max)
  }
  else
  {
    maxd[[ctr+2]] <- aggregate(mam$depth ~ year(mam$date),mam, FUN=max)[,-1]
  }
  ctr=ctr+1
}


colnames(maxd) <- c("Year",paste(round(time,digits = 1), "hour"))

#Choose between GEV or Gumbel method
RP <- c(5,10,50,100,200,500) 
output <- data.frame(matrix(ncol = length(dur), nrow = length(RP)))

for (i in 1:length(dur)){
  val = as.numeric(maxd[,1+i])
  fitgev <- fevd(val, maxd, type="GEV")
  fitgmb <- fevd(val, maxd, type="Gumbel")
  fit1 <- summary(fitgev)
  fit2 <- summary(fitgmb)
  
  #Choosing the best fit using AIC and BIC
  if (fit1$AIC & fit1$BIC < fit2$AIC & fit2$BIC){
    mlefit <- fitgev}
  else {mlefit <- fitgmb}
  rlmle <- return.level(mlefit, conf = 0.05, return.period = RP)
  output[i] <- as.numeric(rlmle)
} 
output <- data.frame(t(output))
colnames(output) <- RP
rownames(output) <- paste(round(time,digits = 1), "hour")

#Plotting the grapgh
ggplot(output, aes(time, y = value, color = variable)) + 
  geom_line(aes(y = output$`5`, col = "T5"),size=1) +
  geom_line(aes(y = output$`10`, col = "T10"),size=1) +
  geom_line(aes(y = output$`50`, col = "T50"),size=1) +
  geom_line(aes(y = output$`100`, col = "T100"),size=1) + 
  geom_line(aes(y = output$`200`, col = "T200"),size=1) +
  geom_line(aes(y = output$`500`, col = "T500"),size=1) +
  ylab("Intensity [mm/hour]") +
  xlab("Duration [hours]") +
  labs(title = "Intensity-Duration-Frequency Curve") +
  labs(color = "Return Periods")+
  theme(panel.background = element_rect(fill = "white"))+
  theme(legend.background = element_rect(fill= "white", size=0.5, 
                                         linetype="solid", colour ="black")) 

