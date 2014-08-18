
load (file ="../Processed Data/AllDetections.RData")
load (file = "../Processed Data/WSDetections.RData")
load (file = "../Processed Data/EncountersWeek.RData")

library (lubridate)
library (plyr)
library (ggplot2)

## DAILY INDIVIDUALS  -------------------------------------------

# DAILY PRESENCE ABSENCE
ACO.WS.DAILY.PRES <- with (DET.WS, as.data.frame (table (TRANSMITTERID, cut (DATETIME, "day"))))
names(ACO.WS.DAILY.PRES) <- c ("TRANSMITTERID", "DATE", "DETECTION.COUNT")
ACO.WS.DAILY.PRES <- ACO.WS.DAILY.PRES[ACO.WS.DAILY.PRES$DETECTION.COUNT != 0, c(1, 2)]
ACO.WS.DAILY.PRES$DATE <- as.Date (ACO.WS.DAILY.PRES$DATE)
# Calculate periods
ACO.WS.DAILY.PRES$PERIOD.DAY <- ACO.WS.DAILY.PRES$DATE - min (ACO.WS.DAILY.PRES$DATE) + 1
ACO.WS.DAILY.PRES$PERIOD.SEC <- (ACO.WS.DAILY.PRES$DATE - min (ACO.WS.DAILY.PRES$DATE)) * 24 * 3600 + 1

# Convert transmitter ID to characters to enable marging
ACO.WS.DAILY.PRES$TRANSMITTERID <- as.character (ACO.WS.DAILY.PRES$TRANSMITTERID)
WS.TAGS$TRANSMITTERID <- as.character (WS.TAGS$TRANSMITTERID)
# Merge dataframes and order
ACO.WS.DAILY.PRES <- merge (ACO.WS.DAILY.PRES, WS.TAGS[, c (3, 6)], by = "TRANSMITTERID")
ACO.WS.DAILY.PRES <- ACO.WS.DAILY.PRES[order (ACO.WS.DAILY.PRES$DATE), ]

## NORMAL DETECTIONS
ACO.WS.PRES <- DET.WS[, c (1, 2)]
ACO.WS.PRES$PERIOD.SEC <- round (as.numeric (ACO.WS.PRES$DATETIME - min (ACO.WS.PRES$DATETIME)) + 1)
ACO.WS.PRES <- merge (ACO.WS.PRES, WS.TAGS[, c (3, 6)], by = "TRANSMITTERID")
ACO.WS.PRES <- ACO.WS.PRES[order (ACO.WS.PRES$DATETIME), ]

# WEEKLY INDIVIUDUALS IDENTIFIED ------------------------------------------

weeks <- sort (as.Date (cut (ARRAY.EVENTS$DATETIME, 'week')))
ARRAY.WEEK <- expand.grid (DATE.WEEKS = seq (min(weeks), max(weeks), 'week'), STATIONNAME = levels (ARRAY.EVENTS$STATIONNAME), RECORDING = FALSE)

# Assign station and location
for (i in length (ARRAY.EVENTS [,1]):1){  # For each event
  # If is a retrieval change the recording status to true for the past
  if (ARRAY.EVENTS$EVENT[i] == "RET"){  
    message ("    Analyzing retrieval of ", ARRAY.EVENTS$RECEIVERID[i], " on ", 
             floor_date (ARRAY.EVENTS$DATETIME[i], "day"), " at   Station ", ARRAY.EVENTS$STATIONNAME[i])
    replace.index <- (as.character (ARRAY.EVENTS$STATIONNAME[i]) == as.character(ARRAY.WEEK$STATIONNAME)) & (ARRAY.WEEK$DATE.WEEKS <= as.Date(ARRAY.EVENTS$DATETIME[i]))
    # Insert Record
    ARRAY.WEEK$RECORDING [replace.index] <- TRUE
  }
  # If is a deployment change the recording status to false for the past
  else {  
    message ("    Analyzing deployment  of ", ARRAY.EVENTS$RECEIVERID[i], " on ", 
             floor_date (ARRAY.EVENTS$DATETIME[i], "day"), " from Station ", ARRAY.EVENTS$STA[i])
    replace.index <- (as.character (ARRAY.EVENTS$STATIONNAME[i]) == as.character(ARRAY.WEEK$STATIONNAME)) & (ARRAY.WEEK$DATE.WEEKS < as.Date(ARRAY.EVENTS$DATETIME[i]))
    # Insert Record
    ARRAY.WEEK$RECORDING [replace.index] <- FALSE
  }
}

ARRAY.WEEK.TABLE <- with (ARRAY.WEEK[ARRAY.WEEK$RECORDING == TRUE, ], table (DATE.WEEKS, STATIONNAME))
DETECTIONS.WEEK <- data.frame (DATE.WEEKS = seq (min (as.Date (rownames(ARRAY.WEEK.TABLE))), max (as.Date (rownames(ARRAY.WEEK.TABLE))), 'week'), 
                               STATIONS.RECORDING = apply(ARRAY.WEEK.TABLE, 1, paste, collapse = ""), 
                               NSTATIONS = apply(ARRAY.WEEK.TABLE, 1, sum))
# Merge array info with individuals info 
DETECTIONS.WEEK <- merge (DETECTIONS.WEEK, data.frame (DATE.WEEKS = as.Date (names (with (ACO.WS.PRES, tapply (ECOCEAN, cut(DATETIME, 'week'), unique)))[-1]),
                                                       N.INDIVIDUALS = sapply (with (ACO.WS.PRES, tapply (ECOCEAN, cut(DATETIME, 'week'), unique)), length)[-1]))
# Remove first last week because of inconsistensy
DETECTIONS.WEEK <- head (DETECTIONS.WEEK, n = -1)
DETECTIONS.WEEK <- tail(DETECTIONS.WEEK, n = -1)

# Include number of sharks tagged
DETECTIONS.WEEK$SHARKS.TAGGED <- NA
for (i in 1:length (DETECTIONS.WEEK[, 1])){
  DETECTIONS.WEEK$SHARKS.TAGGED[i] <- sum (DETECTIONS.WEEK$DATE.WEEKS[i] > as.Date (cut(WS.TAGS$DATE, 'week')))
}

# Calculate raw number of detections
AUX <- as.data.frame (table(cut(DET.WS$DATETIME, 'week')))
names (AUX) <- c ("DATE.WEEKS", "N.DETECTIONS")
AUX$DATE.WEEKS <- as.Date (AUX$DATE.WEEKS)
DETECTIONS.WEEK <- merge (DETECTIONS.WEEK, AUX)

# # Manually remove sharks known to have lost the tag
# TZ-030  2013-10-14
# TZ-068	2013-12-21
# TZ-040	2014-01-07
# TZ-032	2014-01-09

TAG.LOST <- data.frame (ECOCEAN = c ("X","TZ-030", "TZ-068", "TZ-040", "TZ-032"), 
                        DATE = as.Date (c ("2012-10-26","2013-10-14", "2013-12-21", "2014-01-07", "2014-01-09")))

# Cycle over lost tags and reduce number of tagged sharks over it
for (i in 1:nrow(TAG.LOST)){
  DETECTIONS.WEEK$SHARKS.TAGGED[DETECTIONS.WEEK$DATE.WEEKS >= TAG.LOST$DATE[i]] <- DETECTIONS.WEEK$SHARKS.TAGGED[DETECTIONS.WEEK$DATE.WEEKS >= TAG.LOST$DATE[i]] - 1
}

# Calculate when sharks "enter" and leave the array
DET.WS <- merge (DET.WS, WS.TAGS)
FIRST.LAST.DET <- ddply(DET.WS, "ECOCEAN", summarize, LAST.DET = max(DATETIME), FIRST.DET = min (DATETIME))
FIRST.LAST.DET$INTERVAL <- interval (FIRST.LAST.DET$FIRST.DET, FIRST.LAST.DET$LAST.DET)
# Cycle trough each shark
DETECTIONS.WEEK$SHARKS.DET <- apply (daply (FIRST.LAST.DET, "ECOCEAN", function (x, y) y$DATE.WEEKS %within% x$INTERVAL[1], y = DETECTIONS.WEEK), 2, sum)

# Acoustic exposure
DETECTIONS.WEEK$EXP <- with (DETECTIONS.WEEK, NSTATIONS * SHARKS.TAGGED)
DETECTIONS.WEEK$EXP.D <- with (DETECTIONS.WEEK, NSTATIONS * SHARKS.DET)
DETECTIONS.WEEK$INDbySTA <- with (DETECTIONS.WEEK, N.INDIVIDUALS/NSTATIONS)
DETECTIONS.WEEK$INDbyTAG <- with (DETECTIONS.WEEK, N.INDIVIDUALS/SHARKS.TAGGED)
DETECTIONS.WEEK$INDbyEXP <- with (DETECTIONS.WEEK, N.INDIVIDUALS/EXP)
DETECTIONS.WEEK$INDbyEXP.D <- with (DETECTIONS.WEEK, N.INDIVIDUALS/EXP.D)

save (ACO.WS.DAILY.PRES, DETECTIONS.WEEK, file = "../Processed Data/AcousticPresence.RData")
# 
# pdf ("../Figures/Individuals.pdf", width = 3*2, height = 3/sqrt(2), pointsize = 1)
# ggplot() + 
#   #geom_line (aes (y = scale (INDbySTA, center = F, scale = max (INDbySTA))), colour = 2) + 
#   #geom_line (aes (y = scale (INDbyTAG, center = F, scale = max (INDbyTAG))), colour = 3) + 
#   geom_point (data = DETECTIONS.WEEK, aes (x = DATE.WEEKS ,y = scale (INDbyEXP, center = F, scale = max (INDbyEXP))), colour = 'black') +
#   # geom_point (data = DETECTIONS.WEEK, aes (x = DATE.WEEKS ,y = scale (INDbyEXP.D, center = F, scale = max (INDbyEXP.D))), colour = 'grey30') +
#   geom_line (data = DETECTIONS.WEEK, aes (x = DATE.WEEKS ,y = scale (INDbyEXP, center = F, scale = max (INDbyEXP))), colour = 'black', linetype = 1) +
#   #geom_line (data = DETECTIONS.WEEK, aes (x = DATE.WEEKS ,y = scale (INDbyEXP.D, center = F, scale = max (INDbyEXP.D))), colour = 'grey30') +
# #  geom_line(aes(y = scale (N.INDIVIDUALS, center = F, scale = max (N.INDIVIDUALS))*max (N.INDIVIDUALS))) +
#   geom_line (data = ENCOUNTERS.WEEK, aes(x = DATE.WEEK, y = scale (LTSUE, center = F, scale = max (LTSUE))), colour = 'grey50', linetype = 1) + 
#   geom_point (data = ENCOUNTERS.WEEK, aes (x = DATE.WEEK, y = scale (TSUE, center = F, scale = max (LTSUE))), colour = 'grey50') +
#   theme_classic() + ylab ("Scaled number of individuals") + xlab ("Time") +  theme(text = element_text(size = 10))
# dev.off()

