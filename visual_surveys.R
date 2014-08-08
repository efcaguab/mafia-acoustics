## KAUST - WHALE SHARK MAFIA STUDY Visual surveys analysis 
# Correction of search effort and 

# COMPARE ALL DATA WITH SOCPROG DATA --------------------------------------
encounters <- vector ("list", 2)
encounters[[1]] <- read.csv ("../Raw Data/Encounters_AllDetails_20140721.csv")
encounters[[2]] <- read.csv ("../Raw Data/Ecounters_socprog_20140721.csv")

encounters[[1]]$Notes <- factor (encounters[[1]]$Notes)

# We keep only TZ sharks
encounters[[1]] <- encounters[[1]][grepl ("TZ", encounters[[1]]$Marked.Individual), ]
encounters[[2]] <- encounters[[2]][grepl ("TZ", encounters[[2]]$ID), ]
# The socprog sheet has more sharks, I don't know why but I don't care either so I'll stick to that one
ENCOUNTERS <- encounters[[2]]

# Process input data ------------------------------------------------------

# Process the sheet
names (ENCOUNTERS)[1] <- "DATE"
ENCOUNTERS <- ENCOUNTERS[, c(1,6)]
ENCOUNTERS$DATE <- as.POSIXct (as.character(ENCOUNTERS$DATE), format = "%m/%d/%y %H:%M", tz="Africa/Dar_es_Salaam")
ENCOUNTERS$DATE <- as.Date (ENCOUNTERS$DATE)
ENCOUNTERS <- ENCOUNTERS[order (ENCOUNTERS$DATE), ]
ENCOUNTERS$PERIOD.DAY <- as.numeric (ENCOUNTERS$DATE - min(ENCOUNTERS$DATE)) + 1
ENCOUNTERS$ECOCEAN <- gsub('^(.{2})(.*)$', '\\1-\\2', ENCOUNTERS$ID)  
# Do one only for the tagged sharks 
ENCOUNTERS.TAGGED <- ENCOUNTERS[!is.na (match (ENCOUNTERS$ECOCEAN, WS.TAGS$ECOCEAN)), ]

# Process wind data
WIND <- read.csv ("../Raw Data/SatelliteBlendedWindSpeeds_20140517.csv")
WIND$DATE <- as.Date (WIND$DATE, format = '%d/%m/%Y')
WIND$WIND.SPEED[WIND$WIND.SPEED == -9999] <- NA
WIND.WEEK <- data.frame (DATE.WEEK = as.Date(unique(cut (WIND$DATE, 'week'))), WIND.SPEED = as.vector (tapply (WIND$WIND.SPEED, cut (WIND$DATE, 'week'), mean, na.rm = TRUE)))

# Process survey data
SURVEYS <- read.csv ("../Raw Data/Surveys_20140721.csv")
SURVEYS$START.TIME <- as.POSIXct (with (SURVEYS, paste (DATE, START.TIME)), format = "%d-%b-%y %H:%M",  tz = "Africa/Dar_es_Salaam")
SURVEYS$END.TIME <- as.POSIXct (with (SURVEYS, paste (DATE, END.TIME)), format = "%d-%b-%y %H:%M",  tz = "Africa/Dar_es_Salaam") 
SURVEYS$DATE <- as.Date (SURVEYS$DATE, format ='%d-%b-%y')
SURVEYS$DURATION <- with (SURVEYS, END.TIME-START.TIME)
SURVEYS$DURATION[is.na (SURVEYS$DURATION)] <- round (mean (SURVEYS$DURATION, na.rm = TRUE))
# Aggregate by week
SURVEYS.WEEK <- data.frame(DATE.WEEK = as.Date (names(with (SURVEYS, tapply (DURATION, cut(DATE, 'week'), mean, na.rm = TRUE)))),
                           S.HOURS = as.numeric (with (SURVEYS, tapply (DURATION, cut(DATE, 'week'), sum, na.rm = TRUE)))/60,
                           M.HOURS = as.numeric (with (SURVEYS, tapply (DURATION, cut(DATE, 'week'), mean, na.rm = TRUE)))/60,
                           OBSERVERS = as.numeric (with (SURVEYS, tapply (as.numeric (SNORKELERS), cut(DATE, 'week'), mean, na.rm = TRUE))),
                           N.SURVEYS = as.vector (table (cut (SURVEYS$DATE, 'week'))))

start.date <- "2012-09-10 00:00"
end.date <- "2014-11-10 00:00:00"

# Put surveys and wind data together
SURVEYS.WEEK <- merge(SURVEYS.WEEK, WIND.WEEK, all.x = TRUE)
SURVEYS.WEEK <- subset (SURVEYS.WEEK, DATE.WEEK >= start.date & DATE.WEEK <= end.date)
# Artificially add weeks in off season
SURVEYS.WEEK <- merge (data.frame (DATE.WEEK = as.Date (unique (cut (SURVEYS.WEEK$DATE, 'week')))), SURVEYS.WEEK,  all.x = TRUE)
SURVEYS.WEEK$S.HOURS[is.na (SURVEYS.WEEK$S.HOURS)] <- 0
SURVEYS.WEEK$OBSERVERS[is.na (SURVEYS.WEEK$OBSERVERS)] <- 0
SURVEYS.WEEK$WIND.SPEED[is.na (SURVEYS.WEEK$WIND.SPEED)] <- 6.8

# Number of tagged individuals sighted per week
ENCOUNTERS.WEEK <- data.frame (DATE.WEEK = as.Date (names (sapply (with (subset (ENCOUNTERS.TAGGED, DATE > start.date & DATE < end.date), tapply(ECOCEAN, cut(DATE, 'week'), unique)), length))),  
                              N.INDIVIDUALS.TAGGED = as.vector (sapply (with (subset (ENCOUNTERS.TAGGED, DATE > start.date & DATE < end.date), tapply(ECOCEAN, cut(DATE, 'week'), unique)), length)),
                              N.INDIVIDUALS = as.vector (sapply (with (subset (ENCOUNTERS, DATE > start.date & DATE < end.date), tapply(ECOCEAN, cut(DATE, 'week'), unique)), length)), 
                              N.ENCOUNTERS.TAGGED = as.vector(with (subset (ENCOUNTERS.TAGGED, DATE > start.date & DATE < end.date), tapply(ECOCEAN, cut(DATE, 'week'), length))),
                              N.ENCOUNTERS = as.vector(with (subset (ENCOUNTERS, DATE > start.date & DATE < end.date), tapply(ECOCEAN, cut(DATE, 'week'), length))))

# Include search effort data
ENCOUNTERS.WEEK <- merge(SURVEYS.WEEK, ENCOUNTERS.WEEK, by = "DATE.WEEK", all.x = TRUE)
ENCOUNTERS.WEEK$N.INDIVIDUALS[is.na (ENCOUNTERS.WEEK$N.INDIVIDUALS)] <- 0
ENCOUNTERS.WEEK$N.INDIVIDUALS.TAGGED[is.na (ENCOUNTERS.WEEK$N.INDIVIDUALS.TAGGED)] <- 0
ENCOUNTERS.WEEK$N.ENCOUNTERS.TAGGED[is.na (ENCOUNTERS.WEEK$N.ENCOUNTERS.TAGGED)] <- 0
ENCOUNTERS.WEEK$N.ENCOUNTERS[is.na (ENCOUNTERS.WEEK$N.ENCOUNTERS)] <- 0
ENCOUNTERS.WEEK$S.HOURS.SCALED <- scale (ENCOUNTERS.WEEK$S.HOURS)
ENCOUNTERS.WEEK$M.HOURS[is.na (ENCOUNTERS.WEEK$M.HOURS)] <- 0
ENCOUNTERS.WEEK$M.HOURS.SCALED <- scale (ENCOUNTERS.WEEK$M.HOURS)
ENCOUNTERS.WEEK$OBSERVERS.SCALED <- scale (ENCOUNTERS.WEEK$OBSERVERS)
ENCOUNTERS.WEEK$WIND.SPEED.SCALED <- scale (ENCOUNTERS.WEEK$WIND.SPEED)
ENCOUNTERS.WEEK$SUE <- with (ENCOUNTERS.WEEK, N.INDIVIDUALS/S.HOURS)
ENCOUNTERS.WEEK$TSUE <- with (ENCOUNTERS.WEEK, N.INDIVIDUALS.TAGGED/S.HOURS)
ENCOUNTERS.WEEK$EUE <- with (ENCOUNTERS.WEEK, N.ENCOUNTERS/S.HOURS)
ENCOUNTERS.WEEK$TEUE <- with (ENCOUNTERS.WEEK, N.ENCOUNTERS.TAGGED/S.HOURS)
ENCOUNTERS.WEEK$LSUE <- ENCOUNTERS.WEEK$SUE
ENCOUNTERS.WEEK$LTSUE <- ENCOUNTERS.WEEK$TSUE
ENCOUNTERS.WEEK$LEUE <- ENCOUNTERS.WEEK$EUE
ENCOUNTERS.WEEK$LTEUE <- ENCOUNTERS.WEEK$TEUE
ENCOUNTERS.WEEK$LSUE[is.na (ENCOUNTERS.WEEK$LSUE)] <- 0
ENCOUNTERS.WEEK$LTSUE[is.na (ENCOUNTERS.WEEK$LTSUE)] <- 0
ENCOUNTERS.WEEK$LEUE[is.na (ENCOUNTERS.WEEK$LEUE)] <- 0
ENCOUNTERS.WEEK$LTEUE[is.na (ENCOUNTERS.WEEK$LTEUE)] <- 0
# Include number of encounters

pdf ("../Figures/Surveys.pdf", width = 3, height = 3/sqrt(2), pointsize = 1)
ggplot(ENCOUNTERS.WEEK, aes (x = DATE.WEEK)) + 
  geom_line (aes(y = LTSUE * max (N.INDIVIDUALS.TAGGED))) + geom_point (aes (y = TSUE* max (N.INDIVIDUALS.TAGGED))) +
#  geom_line (aes (y = LTEUE), colour = 'gray') + geom_point (aes(y = TEUE), colour = 'gray') + 
  theme_classic() + ylab ("Encounters per hour") + xlab ("Time") +  theme(text = element_text(size = 10))
dev.off()

save (ENCOUNTERS.WEEK, file = "../Processed Data/EncountersWeek.RData")
