load (file = "../Processed Data/WSDetections.RData")

library (mgcv)
library(lubridate)
library (maptools)
library (ggplot2)

# Merge with shark info
DET.WS <- merge (DET.WS, WS.TAGS)
# Convert raw sensor values to meters
DET.WS$DEPTH <- (DET.WS$SENSOR1 - 2) *3
DET.WS$DATETIME.NUM <- (as.numeric(DET.WS$DATETIME) - as.numeric (min (DET.WS$DATETIME)))/3600/24
DET.WS.DEPTH <- DET.WS[!is.na(DET.WS$DEPTH), ]
# Explore using factores
dmodel <- lme (DEPTH ~ as.factor(month(DATETIME)) + as.factor(hour(DATETIME)) + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH)
# Residuals are not autocorrelated. Except slightly for TZ-032 (9th level)
residuals <- resid(dmodel, type = "normalized")
tapply (residuals, DET.WS.DEPTH$ECOCEAN, acf)
# There is a difference between day/night the changes occur somewhere between 6
# and 7 and between 18 and 19. Sunrise and sunset times should do good.
qplot(row.names(as.data.frame(dmodel$coefficients$fixed[grepl("hour", names(dmodel$coefficients$fixed))])), 
      as.data.frame(dmodel$coefficients$fixed[grepl("hour", names(dmodel$coefficients$fixed))])[,1], 
      geom = "bar", stat = "identity") + coord_flip()
# The change between months is more progressive. Dec, Jan and Feb are virtually
# the same. Then it grows quickly in March and April (max) and then decreasses
# progressively back until november, december. For the final model we should use categories or a smoother 
qplot(row.names(as.data.frame(dmodel$coefficients$fixed[grepl("month", names(dmodel$coefficients$fixed))])), 
      as.data.frame(dmodel$coefficients$fixed[grepl("month", names(dmodel$coefficients$fixed))])[,1], 
      geom = "bar", stat = "identity") + coord_flip()
# We check that the random factor is neccessary. We conclude that it's much
# better to consider the shark as a random factor
dmodel.2 <- lm (DEPTH ~ as.factor(month(DATETIME)) + as.factor(hour(DATETIME)) + SEX + SIZE, data = DET.WS.DEPTH)
anova(dmodel, dmodel.2)
# We compare residuals with predictors and responses. They are all fine. No
# different variances required for different strata :)
qplot (as.factor (month(DET.WS.DEPTH$DATETIME)), residuals, geom = "boxplot")
qplot (as.factor (hour(DET.WS.DEPTH$DATETIME)), residuals, geom = "boxplot")

m.coord <- matrix (39.6, -7.9, nrow = 1, ncol = 2)
# We calculate day night categories
pb <- txtProgressBar(min = 0, max = nrow(DET.WS.DEPTH), style = 3)
for (i in 1:nrow(DET.WS.DEPTH)){
  sunrise <- sunriset(m.coord, DET.WS.DEPTH$DATETIME[i], POSIXct.out = TRUE, direction = "sunrise")$time
  sunset <- sunriset(m.coord, DET.WS.DEPTH$DATETIME[i], POSIXct.out = TRUE, direction = "sunset")$time
  c.dawn <- crepuscule (m.coord, DET.WS.DEPTH$DATETIME[i], POSIXct.out = TRUE, direction = "dawn", solarDep = 6)$time
  c.dusk <- crepuscule (m.coord, DET.WS.DEPTH$DATETIME[i], POSIXct.out = TRUE, direction = "dusk", solarDep = 6)$time
  if (DET.WS.DEPTH$DATETIME[i] < sunrise | DET.WS.DEPTH$DATETIME[i] > sunset) {
    DET.WS.DEPTH$DAY.NIGHT.S[i] <- "NIGHT" 
  } else {
    DET.WS.DEPTH$DAY.NIGHT.S[i] <- "DAY"
  }
  if (DET.WS.DEPTH$DATETIME[i] < c.dawn | DET.WS.DEPTH$DATETIME[i] > c.dusk) {
    DET.WS.DEPTH$DAY.NIGHT.C[i] <- "NIGHT" 
  } else {
    DET.WS.DEPTH$DAY.NIGHT.C[i] <- "DAY"
  }
  setTxtProgressBar(pb, value = i)
}
close (pb)

# Because it takes a lot to calculate this shit we'll save the object
# save (DET.WS.DEPTH, file = "../Processed Data/DepthInfo.RData")

# Create new models with the new factors and test their performance, Crepuscule is better
dmodel.1 <- lme (DEPTH ~ as.factor(month(DATETIME)) + as.factor(hour(DATETIME)) + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
dmodel.3 <- lme (DEPTH ~ as.factor(month(DATETIME)) + DAY.NIGHT.S + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
dmodel.4 <- lme (DEPTH ~ as.factor(month(DATETIME)) + DAY.NIGHT.C + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
anova (dmodel.1, dmodel.3, dmodel.4)
# Is there any interaction? Yes, it is significant but we are not really
# interested in the interaction.. also. The depth at night is always the same or
# higher (the pattern never reverses), so we sould be fine ignoring it
dmodel.5 <- lme (DEPTH ~ as.factor(month(DATETIME)) + DAY.NIGHT.C +as.factor(month(DATETIME)):DAY.NIGHT.C + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
ggplot(DET.WS.DEPTH, aes (x = as.factor(month(DATETIME)))) + geom_boxplot(aes (y = DEPTH, fill = DAY.NIGHT.C))


# See how the gamm with a smoother for month works
DET.WS.DEPTH$YDAY <- yday (DET.WS.DEPTH$DATETIME)
DET.WS.DEPTH$WEEK <- week (DET.WS.DEPTH$DATETIME)
DET.WS.DEPTH$MONTH <- month (DET.WS.DEPTH$DATETIME)
DET.WS.DEPTH$HOUR <- hour (DET.WS.DEPTH$DATETIME) + minute (DET.WS.DEPTH$DATETIME) / 60

# It's good...
dmodel.6 <- gamm (DEPTH ~ s(YDAY) + DAY.NIGHT.C + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH)
dmodel.7 <- gamm (DEPTH ~ s(WEEK) + DAY.NIGHT.C + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH)
dmodel.8 <- gamm (DEPTH ~ s(WEEK) + DAY.NIGHT.C + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH)
dmodel.9 <- gamm (DEPTH ~ s(WEEK) + s(HOUR) + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH)
dmodel.10 <- gamm (DEPTH ~ s(WEEK) + s(HOUR), random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH)
dmodel.11 <- gamm (DEPTH ~ s(WEEK) + DAY.NIGHT.C, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH)
# To easy work use it as a number
DET.WS.DEPTH$DNN <- as.numeric (as.factor(DET.WS.DEPTH$DAY.NIGHT.C))-1
dmodel.12 <- gamm (DEPTH ~ s(WEEK) + DNN, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH)

#Twelve is the one I want, now I have to make them cycle
DET.WS.C.1 <- DET.WS.DEPTH[DET.WS.DEPTH$WEEK <= 5,]
DET.WS.C.1$WEEK <- DET.WS.C.1$WEEK + 53
DET.WS.C.2 <- DET.WS.DEPTH[DET.WS.DEPTH$WEEK >= 48,]
DET.WS.C.2$WEEK <- DET.WS.C.2$WEEK - 53
DET.WS.C <- rbind (DET.WS.C.2, DET.WS.DEPTH, DET.WS.C.1)
dmodel.13 <- gamm (DEPTH ~ s(WEEK) + DNN, random = list (ECOCEAN = ~ 1), data = DET.WS.C)

# Evaluate weights and correlations
residuals <- resid(dmodel.7$gam, type = "pearson")
predictions <- predict (dmodel.7$gam)
qplot(predictions, residuals)
qplot(DET.WS$DATETIME[!is.na(DET.WS$DEPTH)], residuals)

# No significant effect of Sex or Size when accounting for random factor
dmodel.8 <- gamm (DEPTH ~ s(WEEK) + DAY.NIGHT.C, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH)

# Hour as a smoother
PRED.GAM <- expand.grid (WEEK = 1:53, HOUR = seq (0, 24, by = 2))
PRED.GAM$R10 <- predict(dmodel.10$gam, PRED.GAM, se.fit = TRUE)$fit
PRED.GAM$R10SE <- predict(dmodel.10$gam, PRED.GAM, se.fit = TRUE)$se
ggplot(PRED.GAM, aes (x = WEEK, y = R10)) + geom_line (aes (colour = as.factor(floor(HOUR)))) + geom_ribbon (aes (ymin = R10-R10SE, ymax = R10+R10SE, fill = as.factor(floor(HOUR))), alpha = 0.2)
# Hour as a factor/numeric cheat
PRED.GAM <- expand.grid (DATE = seq(as.POSIXct("2012-01-01 00:00:00"), as.POSIXct("2012-12-31 23:59:50"), by ="week"), DNN = c(0, 0.5, 1))
PRED.GAM$WEEK <- week (PRED.GAM$DATE)
PRED.GAM$R13 <- predict(dmodel.13$gam, PRED.GAM, se.fit = TRUE)$fit
PRED.GAM$R13SE <- predict(dmodel.13$gam, PRED.GAM, se.fit = TRUE)$se
pdf ("../Figures/Depth.pdf", width = 3, height = 3/sqrt(2), pointsize = 1)
ggplot() + geom_line (data = PRED.GAM[PRED.GAM$DNN == 0.5, ], aes (x = DATE, y = R13)) + 
  geom_line (data = PRED.GAM[PRED.GAM$DNN == 0.5,], aes (x = DATE, y = R13+R13SE), linetype = 2) +
  geom_line (data = PRED.GAM[PRED.GAM$DNN == 0.5,], aes (x = DATE, y = R13-R13SE), linetype = 2) +
  #geom_line (data = PRED.GAM[PRED.GAM$DNN != 0.5, ], aes (x = WEEK, y = R13, colour = as.factor(DNN))) + 
  theme_classic() + scale_y_continuous(trans = "reverse") + ylab ("Swimming depth (m)") + xlab ("Month") +  theme(text = element_text(size = 10))
dev.off()


# Latex Tables ------------------------------------------------------------

library(xtable)
xtable(as.data.frame(summary(dmodel.13$gam)$s.table)[, -2])
print(xtable (summary(dmodel.13$gam)$p.table))
PRED.LME <- expand.grid (WEEK = 1:53, DAY.NIGHT.C = c ("DAY", "NIGHT"), ECOCEAN = levels (DET.WS.DEPTH$ECOCEAN))

# Models for the anova
dmodel.A <- gamm (DEPTH ~ s(WEEK) + DNN + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
dmodel.B <- gamm (DEPTH ~ s(WEEK) + DNN + SEX, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
dmodel.C <- gamm (DEPTH ~ s(WEEK) + DNN + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
dmodel.D <- gamm (DEPTH ~ s(WEEK) + SEX + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
dmodel.E <- gamm (DEPTH ~ s(WEEK) + DNN, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
dmodel.F <- gamm (DEPTH ~ s(WEEK) + SEX, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
dmodel.G <- gamm (DEPTH ~ s(WEEK) + SIZE, random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")
dmodel.H <- gamm (DEPTH ~ s(WEEK), random = list (ECOCEAN = ~ 1), data = DET.WS.DEPTH, method = "ML")

predictions <- predict (dmodel.7$gam)
qplot(predictions, residuals)
qplot(DET.WS$DATETIME[!is.na(DET.WS$DEPTH)], residuals)

week (DET.WS.DEPTH$DATETIME)
