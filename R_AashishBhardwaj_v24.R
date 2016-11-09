rm(list=ls(all=TRUE))
setwd("C:/temp_AB_upc/001_R_Studio_161014")
library(dplyr)
library(forecast)
library(zoo)
library(ggplot2)
df = read.table("10.1.txt", header = TRUE)
df$TimeStamp <- as.POSIXct(paste(df$Date, df$Time), format="%m/%d/%Y %H:%M")
df$Date <- NULL
df$Time <- NULL
df$Depth <- 0.1
ag<-aggregate(df[,2], list(df$TimeStamp), sum )
colnames(ag)[2]= "DP"
colnames(ag)[1]= "DT"
ts2<-data.frame(seq(ag$DT[1],ag$DT[length(ag$DT)] , by="min"))
colnames(ts2)[1]= "DT"
ts2$DP2 <- 0
ts3<-data.frame(left_join(ts2, ag))
ts3$DP2 <- NULL
ts3[is.na(ts3)]<-0
ts3$m15<-(60*as.numeric(ma(ts3$DP, 15, centre= FALSE)))
ts3$m30<-(60*as.numeric(ma(ts3$DP, 30, centre= FALSE)))
ts3$hr1<-(60*as.numeric(ma(ts3$DP, 60, centre= FALSE)))
ts3$hr2<-(60*as.numeric(ma(ts3$DP, 120, centre= FALSE)))
ts3$hr6<-(60*as.numeric(ma(ts3$DP, 360, centre= FALSE)))
ts3$hr12<-(60*as.numeric(ma(ts3$DP, 720, centre= FALSE)))
ts3$hr24<-(60*as.numeric(ma(ts3$DP, 1440, centre= FALSE)))
ts3[is.na(ts3)]<-0
ts3$year<- as.numeric(strftime(ts3$DT, "%Y"))
ts3<- as.matrix(as.data.frame(lapply(ts3, as.numeric)))
idf<-data.frame(as.numeric(seq(1994, 2016)))
colnames(idf)[1]= "year"
idf$m15<-unlist(aggregate(ts3[,3], list(ts3[,10]),max)[2])
idf$m30<-unlist(aggregate(ts3[,4], list(ts3[,10]),max)[2])
idf$hr1<-unlist(aggregate(ts3[,5], list(ts3[,10]),max)[2])
idf$hr2<-unlist(aggregate(ts3[,6], list(ts3[,10]),max)[2])
idf$hr6<-unlist(aggregate(ts3[,7], list(ts3[,10]),max)[2])
idf$hr12<-unlist(aggregate(ts3[,8], list(ts3[,10]),max)[2])
idf$hr24<-unlist(aggregate(ts3[,9], list(ts3[,10]),max)[2])
t2<-data.frame(c(0.25,0.5,1,2,6,12,24))
colnames(t2)[1]= "hr"
t2$mean<-colMeans(idf)[-1]
t2$sd<-c(sd(idf$m15),sd(idf$m30),sd(idf$hr1),sd(idf$hr2),sd(idf$hr6),sd(idf$hr12),sd(idf$hr24))
t2$rp2 <- (t2$mean + (t2$sd * -0.164))
t2$rp5 <- (t2$mean + (t2$sd * 0.719))
t2$rp10 <- (t2$mean + (t2$sd * 1.305))
t2$rp50 <- (t2$mean + (t2$sd * 2.529))
t2$rp100 <- (t2$mean + (t2$sd * 3.137))
idfgraph <- 
  ggplot(t2, aes(hr)) +
  geom_line(aes(y = rp2, colour = "2 years",), size = 1) + 
  geom_point(aes(y = rp2, colour = "2 years"), size = 1.5) +
  geom_line(aes(y = rp5, colour = "5 years"), size = 1)  +
  geom_point(aes(y = rp5, colour = "5 years"), size = 1.5) + 
  geom_line(aes(y = rp10, colour = "10 years"), size = 1)  +
  geom_point(aes(y = rp10, colour = "10 years"), size = 1.5) + 
  geom_line(aes(y = rp50, colour = "50 years"), size = 1) +
  geom_point(aes(y = rp50, colour = "50 years"), size = 1.5) + 
  geom_line(aes(y = rp100 , colour = "100 years"), size = 1)  +
  geom_point(aes(y = rp100, colour = "100 years"), size = 1.5) +
  ylab("Intensity, mm/hour") +
  xlab("Duration, hours") +
  labs(title = "Intensity-Duration-Frequency curves") +
  theme(panel.background = element_rect(fill = "white smoke"))+
  theme(panel.grid.major = element_line(colour = "black", linetype = "dotted"))+
  theme(panel.grid.minor = element_line(colour = "white smoke", linetype = "dotted"))+
  theme(legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="darkblue")) +
  labs(color = "Return Periods, years:") +
  theme(legend.position = c(.92, .8))
###################################################################################
#To do the goodness of fit using K-S test, reference to http://www.mathematik.uni-kl.de/~schwaar/Exercises/Tabellen/table_kolmogorov.pdf
gf<-idf
com<-data.frame(c(0.25,0.5,1,2,6,12,24)) # this dataframe will store the max value for d1 and d2 for all durations
colnames(com)[1] = "Duration"
for (m in 2:(length(gf))){
  gd <- data.frame(gf[1],gf[m])
  gd <- gd[order(gd[,2]),]
  colnames(gd)[2] = "id"
  colnames(gd)[1] = "DT"
  men<- mean(gd[,2])
  sd<- sd(gd[,2])

    for(n in 1:nrow(gd))
    {
      a<- (sd*sqrt(6)/3.14)
      b<- (men - (a* 0.5772))
      gd[n,3] <- ((gd[n,2]-men)^2)
      gd[n,4] <- ((n-0.44)/(nrow(gd)+0.12)) # calculating gr
      gd[n,5] <- exp(-exp(-(gd[n,2]-b)/a)) # calculating gu
      gd[n,6] <- abs ( gd[n,5] - gd[n,4])#d1
      
        if(n>1){
          gd[n,7] <- abs (gd[n,5] - gd[(n-1),4])#d2
        }else{
          gd[n,7]<- 0
        }
    }
  
  
  com[(m-1),2] <- max(max(gd[6]),max(gd[7]))
  rm(gd)
  if((com[(m-1),2]) < 0.247){
    print("The curve fits good for hour value :  ") 
    print(com[(m-1),1] )
  }else{
    print("The curve do not fits good for hour value :  ") 
    print(com[(m-1),1] )
    }
}
######################################################################################
# displaying goodness of fit graphically
gfmn15 <- data.frame(idf$year,idf$m15)
gfmn15 <- gfmn15[order(gfmn15$idf.m15),]
colnames(gfmn15)[2]= "id"
colnames(gfmn15)[1]= "DT"
mean15 <- mean(gfmn15$id)
sd15<- sd(gfmn15$id)
gfmn15$r <- ((gfmn15$id-mean15)^2)
for (i in 1:length(gfmn15$id))
{ a<- (sd15*sqrt(6)/3.14)
b<- (mean15 - (a*0.5772))
gfmn15$gr[i] <- ((i-0.44)/(length(gfmn15$id) + 0.12))
gfmn15$gu[i] <- exp(-exp(-(gfmn15$id[i]-b)/a))
gfmn15$d1[i] <- abs(gfmn15$gu[i]-gfmn15$gr[i])
if(i>1){
  gfmn15$d2[i] <- abs(gfmn15$gu[i]-gfmn15$gr[i-1])
}else{
  gfmn15$d2[i]<-0
}
}
max <- max(max(gfmn15$d2),max(gfmn15$d1))
gof <- 
  ggplot(gfmn15, aes(id)) +
  geom_line(aes(y = gr , colour = "Gringorten",), size = 1) + 
  geom_point(aes(y = gr, colour = "Gringorten"), size = 1.5) +
  geom_line(aes(y = gu, colour = "Gumbelverteilung"), size = 1)  +
  geom_point(aes(y = gu, colour = "Gumbelverteilung"), size = 1.5) 
