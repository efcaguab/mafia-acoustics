# Anaylysis of range test stuff. I want to evaluate the two range tests simultaniously in the same model but using location as a random factor.


## LOAD LIBRARIES AND PREPROCESSED DATA ----------------------------------

load ("../Processed Data/AllDetections.RData")
library (sp)
library (reshape2)
library (plyr)
library (lme4)

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

# Include distances in the detection data frame
DET.RT <- merge (DET.RT, dist.matrix)

## GROUP IN TWO HOURS INTERVALS ----------------------------------------------

# Cut in two hour guild
DET.RT$TWOH <- cut (DET.RT$DATETIME, "2 hour")
# Join together both stations
DET.RT$STATION.COMB <- as.factor (with (DET.RT, paste (RT.STATION, STATIONNAME, sep = "-")))
# Calculate 2h detection probabilities. The delay is 8 minutes
RT <- ddply (DET.RT, .(TWOH, STATION.COMB), function (x, y) {
  c (DIST = round (x$DIST[1] * 1000),  # Distance in meters
     DP = nrow(x)/15,  # Detection probability
     DEPTH = y$DEPTH[match (x$RT.STATION[1], y$STATION)])
  }, y = RT.STATIONS)

## BUILD MODEL ---------------------------------------------------------------

# Calculate time of the day
RT$TIME <- hour (RT$TWOH)
RT$DAY.NIGTH <- "DAY"
RT$DAY.NIGTH[RT$TIME %in% c(0:4, 18:22)] <- "NIGHT" 

# Depth is a factor not a numeric value
RT$DEPTH <- as.factor (RT$DEPTH)

MODEL <- glm (cbind (DP*15, 15-DP*15) ~ DIST + DEPTH, data = RT, family = "binomial")

# ggplot(range.test, aes(x = DIST, y = DP)) + geom_point() + geom_smooth(aes (colour = as.factor(DEPTH)), method = "glm", family = "binomial")
