## KAUST - WHALE SHARK MAFIA STUDY Clean the data. Correct for clock drift.
# Remove first two days of data for each shark and remove false detections.

## LOAD LIBRARIES ----------------------------------------------------------

library (VTrack)
library (ggplot2)
library (lubridate)

## READ DETCTION DATA ---------------------------------------------------------

# Read CSV detection file and process it with VTrack
MAFIA.DETECTIONS <- ReadInputData (read.csv ("../Raw Data/AllMafiaDetections20140715.csv"), iHoursToAdd = 3)
MAFIA.DETECTIONS$DATETIME <- as.POSIXct (as.POSIXlt (MAFIA.DETECTIONS$DATETIME, tz="Africa/Dar_es_Salaam"))

## CORRECT CLOCK DRIFT --------------------------------------------------------

# Read CSV events file 
RECEIVER.EVENTS <- read.csv ("../Raw Data/AllMafiaEvents20140715.csv")
names (RECEIVER.EVENTS) <- c ("DATETIME", "RECEIVERID", "DESC", "DATA", "UNITS")
RECEIVER.EVENTS$DATETIME <- as.POSIXct (RECEIVER.EVENTS$DATETIME, tz ="UTC")
RECEIVER.EVENTS$DATETIME <- as.POSIXct (as.POSIXlt (RECEIVER.EVENTS$DATETIME, tz="Africa/Dar_es_Salaam"))
PC.TIMES <- RECEIVER.EVENTS [RECEIVER.EVENTS$DESC == "PC Time", c (1, 2, 4)]
# We manually checked that all PC times are in the time-zone GMT+3 so we can go ahead and convert the PC times to POSIXct class
PC.TIMES$DATA <- as.POSIXct (substr (as.character (PC.TIMES$DATA), 1, 19), tz="Africa/Dar_es_Salaam")

# # Plots to check for time zone and computer time mistakes
# ggplot(PC.TIMES) + geom_line (aes (x = DATA, y = as.vector (DATA - DATETIME)/60)) + facet_grid (. ~ RECEIVERID , scale = "fixed")
# U <- PC.TIMES[PC.TIMES$RECEIVERID == "VR2W-104845", ]
# V <- data.frame (DATETIME = seq (min(U$DATETIME), max(U$DATETIME), "week"), DIFF = approx (as.numeric(U$DATA), as.vector (U$DATA - U$DATETIME), seq (min(U$DATETIME), max(U$DATETIME), "week")))
# ggplot() + geom_line (data = U, aes (x = DATA, y = as.vector (DATA - DATETIME)/60)) + geom_point (data = V, aes (x = DATETIME, y = DIFF.y/60))

# Error for Data Upload on 2013-01-27 09:21:07 (receiver time) for VR2W-104847
PC.TIMES$DATA[PC.TIMES$DATETIME == as.POSIXct ("2013-01-27 09:21:07")] <- PC.TIMES$DATA[PC.TIMES$DATETIME == as.POSIXct ("2013-01-27 09:21:07")] + 3600*3
# Error for Data Upload on 2013-01-22 11:37:46 (receiver time) for VR2W-104848
PC.TIMES$DATA[PC.TIMES$DATETIME == as.POSIXct ("2013-01-22 11:37:46")] <- PC.TIMES$DATA[PC.TIMES$DATETIME == as.POSIXct ("2013-01-22 11:37:46")] + 3600*3
# Error for Data Upload on 2013-01-27 08:12:34 (receiver time) for VR2W-109044
PC.TIMES$DATA[PC.TIMES$DATETIME == as.POSIXct ("2013-01-27 08:12:34")] <- PC.TIMES$DATA[PC.TIMES$DATETIME == as.POSIXct ("2013-01-27 08:12:34")] + 3600*3
# Error for Data Upload on 2013-01-27 07:42:54 (receiver time) for VR2W-113484
PC.TIMES$DATA[PC.TIMES$DATETIME == as.POSIXct ("2013-01-27 07:42:54")] <- PC.TIMES$DATA[PC.TIMES$DATETIME == as.POSIXct ("2013-01-27 07:42:54")]  + 3600*3

# Correct time drift 
receiverIDs <- levels (MAFIA.DETECTIONS$RECEIVERID)
pb <- txtProgressBar(max=length (receiverIDs), style = 3)
for (i in 1:length (receiverIDs)){
  setTxtProgressBar (pb, i)
  receiver.PC.TIMES <- PC.TIMES[PC.TIMES$RECEIVERID == receiverIDs[i], ]
  drift <- approx (receiver.PC.TIMES$DATA, receiver.PC.TIMES$DATA - receiver.PC.TIMES$DATETIME, MAFIA.DETECTIONS[MAFIA.DETECTIONS$RECEIVERID == receiverIDs[i], ]$DATETIME)$y
  MAFIA.DETECTIONS[MAFIA.DETECTIONS$RECEIVERID == receiverIDs[i], ]$DATETIME <- MAFIA.DETECTIONS[MAFIA.DETECTIONS$RECEIVERID == receiverIDs[i], ]$DATETIME + drift
}
close (pb)
rm (drift, i, receiverIDs, receiver.PC.TIMES)

## ASSIGN DETECTIONS TO STATIONS ---------------------------------------------

# Read and organize events (retrievals/deployments) data
ARRAY.EVENTS <- read.csv ("../Raw Data/ArrayEvents_20140716.csv")
ARRAY.EVENTS <- ARRAY.EVENTS[(ARRAY.EVENTS$EVENT == "DEP") | (ARRAY.EVENTS$EVENT == "RET"), ]
R.EVE <- data.frame(DATETIME = as.POSIXct (ARRAY.EVENTS$DATE, tz="Africa/Dar_es_Salaam")) 
R.EVE$STATIONNAME <- factor(ARRAY.EVENTS$STATION)
R.EVE$EVENT <- ARRAY.EVENTS$EVENT
R.EVE$RECEIVERID <- ARRAY.EVENTS$REC
ARRAY.EVENTS <- R.EVE

# Read stations file and assign detections to stations
STATIONS <- read.csv ("../Raw Data/Stations_20130205.csv")
# Assign station and location
pb <- txtProgressBar(max=length (ARRAY.EVENTS [,1]), style = 3)
for (i in 1:length (ARRAY.EVENTS [,1])){  # For each event
  # If is a deployment change the station for the future
  if (ARRAY.EVENTS$EVENT[i] == "DEP"){  
    # message ("    Analyzing deployment of ", ARRAY.EVENTS$RECEIVERID[i], " on ", 
    #         floor_date (ARRAY.EVENTS$DATETIME[i], "day"), " at   Station ", ARRAY.EVENTS$STATIONNAME[i])
    replace.index <- (as.character (ARRAY.EVENTS$RECEIVERID[i]) == as.character(MAFIA.DETECTIONS$RECEIVERID)) & (MAFIA.DETECTIONS$DATETIME >= ARRAY.EVENTS$DATETIME[i]) 
    # Include station
    MAFIA.DETECTIONS$STATIONNAME [replace.index] <- as.character (ARRAY.EVENTS$STATIONNAME[i])
  }
  # If is a retrieval delete data for the future
  else {  
    # message ("    Analyzing retrieval  of ", ARRAY.EVENTS$RECEIVERID[i], " on ", 
    #        floor_date (ARRAY.EVENTS$DATETIME[i], "day"), " from Station ", ARRAY.EVENTS$STA[i])
    replace.index <- (as.character (ARRAY.EVENTS$RECEIVERID[i]) == as.character(MAFIA.DETECTIONS$RECEIVERID)) & 
      (MAFIA.DETECTIONS$DATETIME >= ARRAY.EVENTS$DATETIME[i]) 
    MAFIA.DETECTIONS$STATIONNAME[replace.index] <- NA
  }
  setTxtProgressBar (pb, i)
}
close (pb)
# Delete detections outside valid intervals
MAFIA.DETECTIONS <- MAFIA.DETECTIONS[MAFIA.DETECTIONS$STATIONNAME != 'Unknown' & !is.na (MAFIA.DETECTIONS$STATIONNAME), ]

save (MAFIA.DETECTIONS, ARRAY.EVENTS, file ="../Processed Data/AllDetections.RData")
rm (R.EVE, i, replace.index)

## SELECT WHALE SHARK TAGS, REMOVE FALSE DETECTIONS ---------------------------

# Read file with Whale Shark Tag lists
WS.TAGS <- read.csv ("../Raw Data/WSTags_20140515.csv")
WS.TAGS$DATE <- as.POSIXct (WS.TAGS$DATE, format="%d/%m/%Y", tz = "Africa/Dar_es_Salaam")
WS.TAGS$NAME <- WS.TAGS$COMMENT <- WS.TAGS$SHARK <- NULL

# Select only whale shark detections 
DET.WS <- MAFIA.DETECTIONS[!is.na (match (MAFIA.DETECTIONS$TRANSMITTERID, WS.TAGS$TRANSMITTERID)), ]

# Remove detections before 48 hours after tagging date
for (i in 1: nrow(WS.TAGS)){
  next2.days <- WS.TAGS$DATE + 60 * 60 * 24 * 2  # Add two days
  replace.index <- as.character (DET.WS$TRANSMITTERID) == as.character (WS.TAGS$TRANSMITTERID[i])
  # Delete rows that are in the tagging day
  DET.WS <- subset (DET.WS, ! (replace.index & (DET.WS$DATETIME < next2.days[i]))) 
}

# Remove unused tags from the factor list
DET.WS$TRANSMITTERID <- factor (DET.WS$TRANSMITTERID)
rm (i, next2.days, replace.index)

save (DET.WS, WS.TAGS, STATIONS, file = "../Processed Data/WSDetections.RData")

