# Anaylysis of range test stuff. I want to evaluate the two range tests simultaniously in the same model but using location as a random factor.

#setwd ("~/Fernando/Dropbox/Research/Whale Shark Mafia Island/Residency paper/MAP - Analysis/")
## LOAD LIBRARIES AND PREPROCESSED DATA ----------------------------------

load ("../Processed Data/AllDetections.RData")

library (sp)
library (reshape2)
library (plyr)
library (lme4)
library (boot)
library (foreach)
library (lubridate)
library (ggplot2)
library (nlme)
library (grid)

# SELECT RANGE TEST DETECTIONS --------------------------------------------

# Read tag list
RT.TAGS <- read.csv ("../Raw Data/RT_Tags.csv")
# For this paper we'll only keep V16 range test tags
DET.RT <- MAFIA.DETECTIONS[MAFIA.DETECTIONS$TRANSMITTERID %in% RT.TAGS$TRANSMITTERID[RT.TAGS$TYPE == "V16"], ]

## ASSIGN DETECTIONS TO RANGE TEST STATION ---------------------------------

RT.EVENTS <- read.csv ("../Raw Data/RTTagsEvents.csv")  # Read range test event table
# Convert to POSIXct 
RT.EVENTS$DATE <- with (RT.EVENTS, as.POSIXct (DATE, format = "%m/%d/%y %H:%M", tz = "Africa/Dar_es_Salaam"))
RT.EVENTS <- RT.EVENTS[order (RT.EVENTS$DATE), ] # Sort by date
DET.RT$RT.STATION <- NA  # Initialise RT.STATION column
for (i in 1:nrow (RT.EVENTS)){  # Cycle trough events
  if (RT.EVENTS$EVENT[i] == "DEP"){
    # If it's a deployment change the Range Test station for all future detections
    DET.RT$RT.STATION[DET.RT$DATETIME > RT.EVENTS$DATE[i] & DET.RT$TRANSMITTERID %in% RT.EVENTS$TRANSMITTERID[i]] <- as.character (RT.EVENTS$RT.STATION[i])
  } else if (RT.EVENTS$EVENT[i] == "RET") {
    # It it's a retrieval change delete the Range Test station for all future detections
    DET.RT$RT.STATION[DET.RT$DATETIME > RT.EVENTS$DATE[i] & DET.RT$TRANSMITTERID %in% RT.EVENTS$TRANSMITTERID[i]] <- NA
  }
}
DET.RT$RT.STATION <- as.factor (DET.RT$RT.STATION)  # Convert to factor 
DET.RT <- DET.RT[! is.na (DET.RT$RT.STATION), ]  # Remove NA

## CALCULATE DISTANCES ------------------------------------------------------

STATIONS <- read.csv ("../Raw Data/Stations_20130205.csv")

# Read Range test stations table
RT.STATIONS <- read.csv ("../Raw Data/RangeTestStations.csv")
# Calculate distance matrix between the range test stations and the receiver stations
dist.matrix <- spDists (as.matrix (STATIONS[,3:2]),  # Location of array stations as a matrix, first column longitude, second col. latitude
                        as.matrix (RT.STATIONS[,3:2]), longlat = TRUE)
# Put names on the matrix
colnames (dist.matrix) <- RT.STATIONS$STATION
# Convert to data frame
dist.matrix <- as.data.frame (dist.matrix)
dist.matrix$STATIONNAME <- STATIONS$STATION
dist.matrix <- melt (dist.matrix, id.vars = "STATIONNAME", variable.name= "RT.STATION", value.name= "DIST")
dist.matrix$STATION.COMB <- with (dist.matrix, paste (RT.STATION, STATIONNAME, sep = "-"))
# Include distances in the detection data frame
DET.RT <- merge (DET.RT, dist.matrix)

## GROUP IN TWO HOURS INTERVALS ----------------------------------------------

# Cut in two hour guild
DET.RT$TWOH <- cut (DET.RT$DATETIME, "2 hour")
# Join together both stations
DET.RT$STATION.COMB <- as.factor (with (DET.RT, paste (RT.STATION, STATIONNAME, sep = "-")))
# Calculate 2h detection probabilities. The delay is 8 minutes
DET.RT.A <- DET.RT[grepl("A", DET.RT$RT.STATION), ]
DET.RT.A$TWOH <- factor (DET.RT.A$TWOH)
DET.RT.A$STATION.COMB <- factor (DET.RT.A$STATION.COMB)
DET.RT.B <- DET.RT[grepl("B", DET.RT$RT.STATION), ]
DET.RT.B$TWOH <- factor (DET.RT.B$TWOH)
DET.RT.B$STATION.COMB <- factor (DET.RT.B$STATION.COMB)

RT.A <- ddply (DET.RT.A, .(TWOH, STATION.COMB), function (x) c (DP = nrow(x)/15), .drop = FALSE)
RT.B <- ddply (DET.RT.B, .(TWOH, STATION.COMB), function (x) c (DP = nrow(x)/15), .drop = FALSE)
RT <- rbind(RT.A, RT.B)
# Include distance stuff
RT$STATION.COMB <- as.character (RT$STATION.COMB)
RT <- merge (RT, dist.matrix[, 3:4])
RT$DIST <- round(RT$DIST * 1000)
# Include the corresponding depth 
for (i in 1:nrow(RT)){
  RT$DEPTH[i] <- RT.STATIONS$DEPTH[match (strsplit(RT$STATION.COMB[i], "-")[[1]][1], RT.STATIONS$STATION)]
}
RT$DEPTH.F <- as.factor (RT$DEPTH)
## BUILD MODEL ---------------------------------------------------------------

# Load wind data. We have to rename the columns manually because we have to
# delete the first three rows Field,Timestamp (dd/mm/yyy hh:mm),Barometer
# (hPa),Temp (deg C),Rel Humidity (%),Wind (km/h),Gust (km/h),Wind Dir
# (deg),Wind Gust Dir (deg),Rain (mm),Solar Rad
# (MJ/m2/h),Evaporation/Transpiration (bare) (mm)
WIND <- read.csv ("../Raw Data/mafia-wx-records-28-11-2012-05-02-2013.csv", skip = 3, header = FALSE)
names (WIND) <- c ("FIELD", "DATE", "PRESSURE", "TEMP", "REL.HUM", "WIND.SPEED", "GUST.SPEED", "WIND.DIR", "GUST.DIR", "RAIN", "SOLAR.RAD", "EVAP.TRANSP")
WIND$DATE <- as.POSIXct (WIND$DATE, format = "%d/%m/%Y %H:%M", tz = "Africa/Dar_es_Salaam")
# Put in two hour intervals
WIND$TWOH <- cut (WIND$DATE, "2 hour")
# Calculate mean values
WIND <- ddply (WIND, "TWOH", summarize, WIND = mean (WIND.SPEED), RAIN = mean (RAIN))

# Calculate time of the day
RT$TIME <- hour (RT$TWOH)
RT$DAY.NIGHT <- "DAY"
RT$DAY.NIGHT[RT$TIME %in% c(0:4, 18:22)] <- "NIGHT" 
# As numeric
RT$DAY.NIGHT.NUM[RT$DAY.NIGHT == "DAY"] <- 1
RT$DAY.NIGHT.NUM[RT$DAY.NIGHT == "NIGHT"] <- 0

# Merge weather and DT data
RT <- merge (RT, WIND)
RT$TIME <- as.numeric (RT$TWOH)
RT$DPLOGIT <- logit(RT$DP)
RT$DPLOGIT[RT$DPLOGIT == Inf] <- 10
RT$DPLOGIT[RT$DPLOGIT == -Inf] <- -10


# # Model Selection
# # Chose the random structure: The best is with autocorrelation and depth as a random factor as well
# MODEL.2 <- gls (DPLOGIT ~ DIST  + WIND + RAIN + DAY.NIGHT + DEPTH, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH))
# MODEL.3 <- gls (DPLOGIT ~ DIST  + WIND + RAIN + DAY.NIGHT + DEPTH, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB))
# MODEL.4 <- gls (DPLOGIT ~ DIST  + WIND + RAIN + DAY.NIGHT + DEPTH, data = RT , weights = varIdent (~ DEPTH))
# MODEL.5 <- gls (DPLOGIT ~ DIST  + WIND + RAIN + DAY.NIGHT + DEPTH, data = RT)
# anova (MODEL.2, MODEL.3, MODEL.4, MODEL.5)
# 
MODEL.6 <- gls (DPLOGIT ~ DIST  + WIND + RAIN + DAY.NIGHT , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.7 <- gls (DPLOGIT ~ DIST  + WIND + RAIN + DAY.NIGHT, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.8 <- gls (DPLOGIT ~ DIST  + WIND + RAIN , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.9 <- gls (DPLOGIT ~ DIST  + WIND + DAY.NIGHT , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.10 <- gls (DPLOGIT ~ DIST  + RAIN + DAY.NIGHT , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")

MODEL.11 <- gls (DPLOGIT ~ DIST  + WIND + RAIN, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.12 <- gls (DPLOGIT ~ DIST  + WIND + DAY.NIGHT, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.13 <- gls (DPLOGIT ~ DIST  + RAIN + DAY.NIGHT, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.14 <- gls (DPLOGIT ~ DIST  + WIND , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.15 <- gls (DPLOGIT ~ DIST  + RAIN , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.16 <- gls (DPLOGIT ~ DIST  + DAY.NIGHT , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")

MODEL.17 <- gls (DPLOGIT ~ DIST  + WIND, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.18 <- gls (DPLOGIT ~ DIST  + RAIN, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.19 <- gls (DPLOGIT ~ DIST  + DAY.NIGHT, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
MODEL.20 <- gls (DPLOGIT ~ DIST  , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
# 
A <-anova (MODEL.6, MODEL.7, MODEL.8, MODEL.9, MODEL.10, MODEL.11, MODEL.12, MODEL.13, MODEL.14, MODEL.15, MODEL.16, MODEL.17, MODEL.18, MODEL.19, MODEL.20)
# 
# MODEL.11 <- gls (DPLOGIT ~ DIST  + WIND + RAIN + DAY.NIGHT + DEPTH, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
# MODEL.11 <- gls (DPLOGIT ~ DIST  + WIND + RAIN + DAY.NIGHT + DEPTH, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
# 
# 
# MODEL.5 <- gls (DPLOGIT ~ DIST  + WIND + RAIN + DAY.NIGHT, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), method = "ML")
# MODEL.6 <- gls (DPLOGIT ~ DIST  + WIND  + DAY.NIGHT, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), method = "ML")
# MODEL.7 <- gls (DPLOGIT ~ DIST   + RAIN + DAY.NIGHT, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), method = "ML")
# MODEL.8 <- gls (DPLOGIT ~ DIST  + WIND + RAIN , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), method = "ML")
# MODEL.9 <- gls (DPLOGIT ~ DIST + DAY.NIGHT, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), method = "ML")
# MODEL.10 <- gls (DPLOGIT ~ DIST  + WIND , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), method = "ML")
# MODEL.11 <- gls (DPLOGIT ~ DIST  + RAIN, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), method = "ML")
# MODEL.12 <- gls (DPLOGIT ~ DIST, data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), method = "ML")

# MODEL.3 <- gls (DPLOGIT ~ DIST  + WIND + RAIN + DAY.NIGHT + DEPTH, data = RT)

# MODEL.2 <- glm (cbind (DP*15, 15-DP*15) ~ DIST  + WIND + RAIN + DAY.NIGHT , data = RT, family = "binomial")
# MODEL.3 <- glm (cbind (DP*15, 15-DP*15) ~ DIST  + WIND + RAIN + DAY.NIGHT + DEPTH, data = RT, family = "binomial")
# 
# MODEL.4 <- glm (cbind (DP*15, 15-DP*15) ~ DIST + WIND + DAY.NIGHT + DEPTH, data = RT, family = "binomial")
# MODEL.5 <- glm (cbind (DP*15, 15-DP*15) ~ DIST + RAIN + DAY.NIGHT + DEPTH, data = RT, family = "binomial")
# MODEL.6 <- glm (cbind (DP*15, 15-DP*15) ~ DIST  + DAY.NIGHT + DEPTH, data = RT, family = "binomial")
# MODEL.7 <- glm (cbind (DP*15, 15-DP*15) ~ DIST  + DAY.NIGHT + DEPTH, data = RT, family = "quasibinomial")
# # Model.7 is the chosen one

# PREDICTION --------------------------------------------------------------

# PREDICT <- expand.grid (DIST = 1:1500, DAY.NIGHT.NUM = c (0, 0.5, 1), DEPTH = c (15, 19.5, 24))
PREDICT <- expand.grid (DIST = 1:800)
# Bootstrap function
# prediction <- function (data, indices, pred, pb){  
#   d <- data[indices, ]
#   MODEL <- gls (DPLOGIT ~ DIST  , data = RT, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
#   p <- inv.logit(predict (MODEL, pred))
#   return (p)
# }
# Bootstrap results
library (doMC)
registerDoMC (cores = 30)
data <- RT
pred <- PREDICT
prediction <- foreach (i = 1:nrow(data), .combine = rbind, .inorder = FALSE) %dopar% {
  d <- data[-i, ]
  MODEL <- gls (DPLOGIT ~ DIST  , data = d, correlation = corAR1(form = ~ TIME | STATION.COMB), weights = varIdent (form = ~ 1| DEPTH), method = "ML")
  p <- inv.logit(predict (MODEL, pred))
  return (p)
}
# pred <- boot (data = RT, statistic = prediction, R = 100, stype = "i", pred = PREDICT, parallel = "multicore")

# Calculate mean
PREDICT$DP.mean <- apply (prediction, 2, mean)

low <- function (dist){
  inv.logit (4.478372-0.27447073 + (-0.013095-0.00037585) * dist)
}
high <- function (dist){
  inv.logit (4.478372+0.27447073 + (-0.013095+0.00037585) * dist)
}
PREDICT$DP.low <- low (PREDICT$DIST)
PREDICT$DP.high <- high (PREDICT$DIST)
# Calculate confidence intervals

# # Only mean values using bca, all others using normal sd
# library (doMC)
# registerDoMC (cores = 24)
# conf <- foreach (i=1:nrow (PREDICT), .combine = rbind) %dopar%  {
#   if (PREDICT$DEPTH[i] == 19.5 & PREDICT$DAY.NIGHT.NUM[i] == 0.5){
#     return (boot.ci(pred, type = "bca", index = i)$bca[4:5])
#   } else return (c(NA, NA))
# }

# save (conf, file = "conf.RData")
# load ("conf.RData")

# Fill with values
# PREDICT$DP.low <- conf[ ,1]
# PREDICT$DP.high <- conf[ ,2]
# # For the others just calculate sd
# PREDICT$DP.low[is.na(PREDICT$DP.low)] <- apply (pred$t[, is.na(PREDICT$DP.low)], 2, mean) - 1.95 * apply (pred$t[, is.na(PREDICT$DP.low)], 2, sd)
# PREDICT$DP.high[is.na(PREDICT$DP.high)] <- apply (pred$t[, is.na(PREDICT$DP.high)], 2, mean) + 1.95 * apply (pred$t[, is.na(PREDICT$DP.high)], 2, sd)

ggplot () + 
  geom_hline (aes (yintercept = c(0.5)), linetype = 2, colour = "#555555", size = 0.3) + 
  geom_vline (aes (xintercept = c( 342)), linetype = 2, colour = "#555555", size = 0.3) +
  geom_line (data = PREDICT, aes (y = DP.mean, x = DIST)) +
  geom_ribbon (data = PREDICT, 
               aes(x = DIST, ymax = DP.high, ymin = DP.low), 
               alpha = 0.15) + 
  ylim(0,1) + 
  xlab ("Distance (m)") + ylab ("Detection probability") +
  theme_linedraw () +
  theme (legend.position = c (0,1), 
         text = element_text (family = "serif"), 
         legend.justification = c (0,1), 
         legend.title.align = 0.5, 
         legend.text.align = 0.5, 
         legend.box.just = "left", 
         legend.title = element_text(face = "italic", size = 8),
         axis.title = element_text(size = 8),
         axis.text = element_text(size = 8), 
         text = element_text (),
         legend.text = element_text(size = 8), 
         panel.grid = element_line (colour = "gray"), 
         panel.grid.major = element_line (colour = "gray"), 
         panel.grid.minor = element_line (colour = "gray"), 
         plot.margin = unit (c (1,2,0,0), "pt"))
ggsave ("RT.pdf", width = 9, height = 9/sqrt(2), units = "cm", pointsize = 8)

# pdf ("RangeTest.pdf", width = 3, height = 3/sqrt(2), family = "Times", pointsize = 1)
# ggplot () + geom_line (data = PREDICT[PREDICT$DAY.NIGHT.NUM == 0.5 & PREDICT$DEPTH == 19.5, ], aes (y = DP.mean, x = DIST)) +
#   geom_ribbon (data = PREDICT[PREDICT$DAY.NIGHT.NUM == 0.5 & PREDICT$DEPTH == 19.5, ], aes(x = DIST, ymax = DP.high, ymin = DP.low), alpha = 0.2) + 
#   #geom_line(data = PREDICT[PREDICT$DAY.NIGHT.NUM != 0.5 & PREDICT$DEPTH == 19.5, ], aes (x = DIST, y = DP.mean, linetype = as.factor(DAY.NIGHT.NUM))) + 
#   ylim(0,1) + theme_classic () + theme(legend.position = "none", text = element_text(size = 10)) + xlab ("Distance (m)") + ylab ("Detection probability") +
#   geom_hline (aes (yintercept = c(0.5)), linetype = 2)
# # ggplot(range.test, aes(x = DIST, y = DP)) + geom_point() + geom_smooth(aes (colour = as.factor(DEPTH)), method = "glm", family = "binomial")
# dev.off()

pdf ("RangeTest.pdf", width = 3, height = 3/sqrt(2), family = "Times", pointsize = 6)
ggplot() + geom_line(data = PREDICT[PREDICT$DAY.NIGHT.NUM != 0.5 & PREDICT$DEPTH == 19.5, ], aes (x = DIST, y = DP.mean, linetype = as.factor(DAY.NIGHT.NUM))) + 
  ylim(0,1) + theme_classic () + theme(legend.position = "none", text = element_text(size = 10)) + xlab ("Distance (m)") + ylab ("Detection probability") 
dev.off()


# Latex tables ------------------------------------------------------------

library (xtable)
xtable (summary(MODEL.6)$tTable) %>% print (type = "html")
A <-anova (MODEL.6, MODEL.7, MODEL.8, MODEL.9, MODEL.10, MODEL.11, MODEL.12, MODEL.13, MODEL.14, MODEL.15, MODEL.16, MODEL.17, MODEL.18, MODEL.19, MODEL.20)
call <- unlist(lapply(strsplit(as.character(A$call), ","), function(x) x[1]))
call <- gsub("gls(model = DPLOGIT ~ ", "",  as.character(call), fixed = TRUE)
call <- gsub("DIST", "$d$",  as.character(call), fixed = TRUE)
call <- gsub("WIND", "$w$",  as.character(call), fixed = TRUE)
call <- gsub("RAIN", "$r$",  as.character(call), fixed = TRUE)
call <- gsub("DAY.NIGHT", "$t$",  as.character(call), fixed = TRUE)
call <- gsub("$ + $", " + ",  as.character(call), fixed = TRUE)
A$call <- call
A$dAIC <- A$AIC - min(A$AIC)
row.names(A) <- NULL
xtable(A[, c(1,3,4, 10)]) %>% print (type = "html")
