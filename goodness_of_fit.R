# calculation of goodness of fit by Kolmogorov-Smirnov method for 30 minute duration
# order yearly maximum rainfall intensity by decreasing
fit_30 <- data.frame(min_30[order(min_30, decreasing = FALSE)])
colnames(fit_30)[1] <- "rain"
# Empirical Likelihood Gringorten
fit_30["Gringorten"] <- (index(fit_30)-0.44)/(12 + 0.12)
# Theoretical distribution Gumbel
fit_30["Gumbel"] <- exp(-exp((-fit_30$rain+(mean(fit_30$rain)-0.5772 * ((sd(fit_30$rain) * sqrt(6))/pi)))/((sd(fit_30$rain) * sqrt(6))/pi)))
# calculation of d1, d2, alfa and c values
alfa30 <- (sd(fit_30$rain) * sqrt(6))/pi
fit_30["d1"] <- fit_30$Gumbel - fit_30$Gringorten
fit_30["d2"] <- c(NA, fit_30[2:12,3]-fit_30[1:11,2]) 
dmax_30min <- max(abs(fit_30$d1), abs(na.omit(fit_30$d2)))
# plotting
plot(fit_30$rain, fit_30$Gringorten, xlab = 'rain intensity (mm/hour)', ylab = 'probability', main = '30 minutes')
lines(fit_30$rain, fit_30$Gumbel)

# calculation of goodness of fit by Kolmogorov-Smirnov method for 1 hour duration
# order yearly maximum rainfall intensity by decreasing
fit_1 <- data.frame(hour_1[order(hour_1, decreasing = FALSE)])
colnames(fit_1)[1] <- "rain"
# Empirical Likelihood Gringorten
fit_1["Gringorten"] <- (index(fit_1)-0.44)/(12 + 0.12)
# Theoretical distribution Gumbel
fit_1["Gumbel"] <- exp(-exp((-fit_1$rain+(mean(fit_1$rain)-0.5772 * ((sd(fit_1$rain) * sqrt(6))/pi)))/((sd(fit_1$rain) * sqrt(6))/pi)))
# calculation of d1, d2, alfa and c values
alfa1 <- (sd(fit_1$rain) * sqrt(6))/pi
fit_1["d1"] <- fit_1$Gumbel - fit_1$Gringorten
fit_1["d2"] <- c(NA, fit_1[2:12,3]-fit_1[1:11,2]) 
dmax_1h <- max(abs(fit_1$d1), abs(na.omit(fit_1$d2)))
plot(fit_1$rain, fit_1$Gringorten,  xlab = 'rain intensity (mm/hour)', ylab = 'probability', main = '1 hour')
lines(fit_1$rain, fit_1$Gumbel)


# calculation of goodness of fit by Kolmogorov-Smirnov method for 3 hours duration
# order yearly maximum rainfall intensity by decreasing
fit_3 <- data.frame(hour_3[order(hour_3, decreasing = FALSE)])
colnames(fit_3)[1] <- "rain"
# Empirical Likelihood Gringorten
fit_3["Gringorten"] <- (index(fit_3)-0.44)/(12 + 0.12)
# Theoretical distribution Gumbel
fit_3["Gumbel"] <- exp(-exp((-fit_3$rain+(mean(fit_3$rain)-0.5772 * ((sd(fit_3$rain) * sqrt(6))/pi)))/((sd(fit_3$rain) * sqrt(6))/pi)))
# calculation of d1, d2, alfa and c values
alfa3 <- (sd(fit_3$rain) * sqrt(6))/pi
fit_3["d1"] <- fit_3$Gumbel - fit_3$Gringorten
fit_3["d2"] <- c(NA, fit_3[2:12,3]-fit_3[1:11,2]) 
dmax_3h <- max(abs(fit_3$d1), abs(na.omit(fit_3$d2)))
plot(fit_3$rain, fit_3$Gringorten,  xlab = 'rain intensity (mm/hour)', ylab = 'probability', main = '3 hours')
lines(fit_3$rain, fit_3$Gumbel)

# calculation of goodness of fit by Kolmogorov-Smirnov method for 6 hours duration
# order yearly maximum rainfall intensity by decreasing
fit_6 <- data.frame(hour_6[order(hour_6, decreasing = FALSE)])
colnames(fit_6)[1] <- "rain"
# Empirical Likelihood Gringorten
fit_6["Gringorten"] <- (index(fit_6)-0.44)/(12 + 0.12)
# Theoretical distribution Gumbel
fit_6["Gumbel"] <- exp(-exp((-fit_6$rain+(mean(fit_6$rain)-0.5772 * ((sd(fit_6$rain) * sqrt(6))/pi)))/((sd(fit_6$rain) * sqrt(6))/pi)))
# calculation of d1, d2, alfa and c values
alfa6 <- (sd(fit_6$rain) * sqrt(6))/pi
fit_6["d1"] <- fit_6$Gumbel - fit_6$Gringorten
fit_6["d2"] <- c(NA, fit_6[2:12,3]-fit_6[1:11,2]) 
dmax_6h <- max(abs(fit_6$d1), abs(na.omit(fit_6$d2)))
plot(fit_6$rain, fit_6$Gringorten,  xlab = 'rain intensity (mm/hour)', ylab = 'probability', main = '6 hours')
lines(fit_6$rain, fit_6$Gumbel)

# calculation of goodness of fit by Kolmogorov-Smirnov method for 12 hours duration
# order yearly maximum rainfall intensity by decreasing
fit_12 <- data.frame(hour_12[order(hour_12, decreasing = FALSE)])
colnames(fit_12)[1] <- "rain"
# Empirical Likelihood Gringorten
fit_12["Gringorten"] <- (index(fit_12)-0.44)/(12 + 0.12)
# Theoretical distribution Gumbel
fit_12["Gumbel"] <- exp(-exp((-fit_12$rain+(mean(fit_12$rain)-0.5772 * ((sd(fit_12$rain) * sqrt(6))/pi)))/((sd(fit_12$rain) * sqrt(6))/pi)))
# calculation of d1, d2, alfa and c values
alfa12 <- (sd(fit_12$rain) * sqrt(6))/pi
fit_12["d1"] <- fit_12$Gumbel - fit_12$Gringorten
fit_12["d2"] <- c(NA, fit_12[2:12,3]-fit_12[1:11,2]) 
dmax_12h <- max(abs(fit_12$d1), abs(na.omit(fit_12$d2)))
plot(fit_12$rain, fit_12$Gringorten,  xlab = 'rain intensity (mm/hour)', ylab = 'probability', main = '12 hours')
lines(fit_12$rain, fit_12$Gumbel)

# calculation of goodness of fit by Kolmogorov-Smirnov method for 18 hours duration
# order yearly maximum rainfall intensity by decreasing
fit_18 <- data.frame(hour_18[order(hour_18, decreasing = FALSE)])
colnames(fit_18)[1] <- "rain"
# Empirical Likelihood Gringorten
fit_18["Gringorten"] <- (index(fit_18)-0.44)/(12 + 0.12)
# Theoretical distribution Gumbel
fit_18["Gumbel"] <- exp(-exp((-fit_18$rain+(mean(fit_18$rain)-0.5772 * ((sd(fit_18$rain) * sqrt(6))/pi)))/((sd(fit_18$rain) * sqrt(6))/pi)))
# calculation of d1, d2, alfa and c values
alfa18 <- (sd(fit_18$rain) * sqrt(6))/pi
fit_18["d1"] <- fit_18$Gumbel - fit_18$Gringorten
fit_18["d2"] <- c(NA, fit_18[2:12,3]-fit_18[1:11,2]) 
dmax_18h <- max(abs(fit_18$d1), abs(na.omit(fit_18$d2)))
plot(fit_18$rain, fit_18$Gringorten,  xlab = 'rain intensity (mm/hour)', ylab = 'probability', main = '18 hours')
lines(fit_18$rain, fit_18$Gumbel)

# calculation of goodness of fit by Kolmogorov-Smirnov method for 24 hours duration
# order yearly maximum rainfall intensity by decreasing
fit_24 <- data.frame(hour_24[order(hour_24, decreasing = FALSE)])
colnames(fit_24)[1] <- "rain"
# Empirical Likelihood Gringorten
fit_24["Gringorten"] <- (index(fit_24)-0.44)/(12 + 0.12)
# Theoretical distribution Gumbel
fit_24["Gumbel"] <- exp(-exp((-fit_24$rain+(mean(fit_24$rain)-0.5772 * ((sd(fit_24$rain) * sqrt(6))/pi)))/((sd(fit_24$rain) * sqrt(6))/pi)))
# calculation of d1, d2, alfa and c values
alfa24 <- (sd(fit_24$rain) * sqrt(6))/pi
fit_24["d1"] <- fit_24$Gumbel - fit_24$Gringorten
fit_24["d2"] <- c(NA, fit_24[2:12,3]-fit_24[1:11,2]) 
plot(fit_24$rain, fit_24$Gringorten,  xlab = 'rain intensity (mm/hour)', ylab = 'probability', main = '24 hours')
lines(fit_24$rain, fit_24$Gumbel)

# combining dmax and alfa results
dmax_24h <- max(abs(fit_24$d1), abs(na.omit(fit_24$d2)))
rm(fit_30,fit_1, fit_3,fit_6,fit_12,fit_18,fit_24)
all_dmax <- data.frame(dmax_30min, dmax_1h,dmax_3h,dmax_6h,dmax_12h,dmax_18h,dmax_24h)
all_alfa <- data.frame(alfa30, alfa1,alfa3,alfa6,alfa12,alfa18,alfa24)
rm(dmax_30min, dmax_1h,dmax_3h,dmax_6h,dmax_12h,dmax_18h,dmax_24h)
rm(alfa30, alfa1,alfa3,alfa6,alfa12,alfa18,alfa24)