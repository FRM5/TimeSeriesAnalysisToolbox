## Name      : Hani Maisarah Mohamad
## Assignment: IDF Curve for Flash Flood (2016)

library(extRemes)
library(ggplot2)

# Read data from a .txt file
mydata=read.table("18.txt",header=TRUE,col.names=c("Date", "Time")) 

# Merge and change column format
dt <-strptime(paste(mydata$Date,mydata$Time), "%m/%d/%Y %H:%M:%S")

# Create bins, annual maximum intensity based on user-defined durations in minutes
# Warning: Long run-time for small durations (i.e. 5 minutes)!
tint <- c("10 mins","30 mins","60 mins","120 mins","360 mins","720 mins","1440 mins")
t <- (as.numeric(gsub("([0-9]*).*","\\1",tint))/60)
ctr=0
for (i in tint){
  n <- data.frame(table(cut(dt, breaks = i)))
  #assign(paste("tint", which(tint == i), sep = ""), n) #Activate to check
  n[2] <- n[2]*0.1/t[ctr+1]
  if (ctr < 1){
    x <- aggregate(n$Freq ~ format(as.Date(n$Var1), "%Y"), n, max)}
  else{
    x[[ctr+2]] <- aggregate(n$Freq ~ format(as.Date(n$Var1), "%Y"), n, max)[,-1]}
  ctr=ctr+1
}
colnames(x) <- c("Year",paste(round(t,digits = 1), "hour"))

# GEV or Gumbel method
rp <- c(5,10,50,100,500) # Specify return period
outrp <- data.frame(matrix(ncol = length(tint), nrow = length(rp)))
for (i in 1:length(tint)){
  val = as.numeric(x[,1+i])
  mlefit <- fevd(val, method = "GMLE", type="GEV")
  #mlefit <- fevd(val, x, type="Gumbel") #Activate for Gumbel method
  rlmle <- return.level(mlefit, conf = 0.05, return.period = rp)
  outrp[i] <- as.numeric(rlmle)
} 
outrp <- data.frame(t(outrp))
colnames(outrp) <- rp
rownames(outrp) <- paste(round(t,digits = 1), "hour")

# Plot graph manually
ggplot(outrp, aes(t, y = value, color = variable)) + 
  geom_line(aes(y = outrp$`5`, col = "T5")) +
  geom_line(aes(y = outrp$`10`, col = "T10")) +
  geom_line(aes(y = outrp$`50`, col = "T50")) +
  geom_line(aes(y = outrp$`100`, col = "T100")) + 
  geom_line(aes(y = outrp$`500`, col = "T500")) +
  ylab("Intensity [mm/hr]") +
  xlab("Duration [hours]") +
  labs(title = "Intensity-Duration-Frequency Curve") +
  labs(color = "Return Periods")
