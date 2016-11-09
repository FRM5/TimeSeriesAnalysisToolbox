## Name      : Hani Maisarah Mohamad
## Assignment: IDF Curve for Flash Flood (2016)

library(lubridate)
library(xts)
library(timeDate)
library(extRemes)
library(ggplot2)

mydata=read.table("18.txt",header=TRUE,col.names=c("Date","Time"))
mydata$dt <- as.POSIXct(paste(mydata$Date,mydata$Time),format = "%m/%d/%Y %H:%M")
mydata$depth <- 0.1
mydata$Date <- NULL
mydata$Time <- NULL

#Create uniform time series with 1 min intervals
grpdata <- tapply(mydata$depth, mydata$dt, sum)
mydata_ts <- xts(grpdata,as.POSIXct(rownames(grpdata)))
bymin <- seq(start(mydata_ts), end(mydata_ts), "min")
dummycell <- xts(rep(0,length(bymin)), order.by=bymin)
new_TS <- merge(dummycell, mydata_ts, fill=0, all.x = TRUE)[,2]

rm(mydata,grpdata,mydata_ts,bymin,dummycell) #Remove to free-up memory space

#Apply moving average for user-defined durations and find annual maximum
tdur <- c(5,10,30,60,120,360,720,1440) #Duration in minutes
t <- tdur/60
dma <- data.frame(date = index(new_TS))
ctr=0
for (i in tdur){
  dma$depth <- rollapply(coredata(new_TS),i,sum,fill=0)
  #assign(paste("ma", i, sep = ""), dma) #Activate to display new variables
  dma$depth <- dma$depth/t[ctr+1] #Convert to intensity mm/hr
  if (ctr < 1){
    maxd <- aggregate(dma$depth ~ year(dma$date),dma, FUN=max)}
  else{
    maxd[[ctr+2]] <- aggregate(dma$depth ~ year(dma$date),dma, FUN=max)[,-1]}
  ctr=ctr+1
}
colnames(maxd) <- c("Year",paste(round(t,digits = 1), "hour"))

#GEV or Gumbel method
rp <- c(5,10,50,100,200,500) # Specify return period
outrp <- data.frame(matrix(ncol = length(tdur), nrow = length(rp)))
for (i in 1:length(tdur)){
  val = as.numeric(maxd[,1+i])
  fitgev <- fevd(val, maxd, type="GEV")
  fitgmb <- fevd(val, maxd, type="Gumbel")
  fit1 <- summary(fitgev)
  fit2 <- summary(fitgmb)
  #Choose the best fit using AIC and BIC
  if (fit1$AIC & fit1$BIC < fit2$AIC & fit2$BIC){
    mlefit <- fitgev}
  else {mlefit <- fitgmb}
  rlmle <- return.level(mlefit, conf = 0.05, return.period = rp)
  outrp[i] <- as.numeric(rlmle)
} 
outrp <- data.frame(t(outrp))
colnames(outrp) <- rp
rownames(outrp) <- paste(round(t,digits = 1), "hour")

#Plot graph manually
ggplot(outrp, aes(t, y = value, color = variable)) + 
  geom_line(aes(y = outrp$`5`, col = "T5")) +
  geom_line(aes(y = outrp$`10`, col = "T10")) +
  geom_line(aes(y = outrp$`50`, col = "T50")) +
  geom_line(aes(y = outrp$`100`, col = "T100")) + 
  geom_line(aes(y = outrp$`200`, col = "T200")) +
  geom_line(aes(y = outrp$`500`, col = "T500")) +
  ylab("Intensity [mm/hr]") +
  xlab("Duration [hours]") +
  labs(title = "Intensity-Duration-Frequency Curve") +
  labs(color = "Return Periods")

#Creating synthetic rainfall
idfdata <- data.frame(t, outrp$`100`)
colnames(idfdata) = c("Td","i")
idffit <- nls(i ~ c/(f+Td^e), start=list(c=1, f=1, e=1), data=idfdata)
fitdat <-predict(idffit, newdata = data.frame(Td=c(seq(1,24))))
designR <- data.frame(c(seq(1,24)))
designR[2] <- fitdat
designR[3] <- designR[1] * designR[2]
designR[4] <- c(designR[1,3],diff(as.matrix(designR[3])))

rain <- as.matrix(unlist(designR[4]))
med <- max(rain)
num <- which.max(rain)
Q <- rain[-num]

if (length(rain) %% 2 != 0){ #check if function is odd
  tl <- sort(which(index(Q) %% 2 == 0),decreasing = TRUE)
  hd <- sort(which(index(Q) %% 2 != 0),decreasing = FALSE)
} else {
  hd <- sort(which(index(Q) %% 2 == 0),decreasing = TRUE)
  tl <- sort(which(index(Q) %% 2 != 0),decreasing = FALSE)
}
designR[5] <- c(Q[hd],med,Q[tl])
colnames(designR) = c("Duration","Intensity","Cum.Depth","Inc.Depth","Precipitation")

e <- as.numeric(as.matrix(designR[5]))
barplot(e, main="T100 - 24 hours Synthetic Precipitation", xlab="Time (hr)", ylab="Precipitation (mm)")
