# Yearly model

load (file = "../Processed Data/AcousticPresence.RData")
load (file = "../Processed Data/WSDetections.RData")
load (file = "../Processed Data/EncountersWeek.RData")
load (file ="../Processed Data/AllDetections.RData")

library (plyr)
library (dplyr)
library (magrittr)
library (lubridate)
library (foreach) 
library (ggplot2)
library (mgcv)
library (doMC)
registerDoMC (cores = 14)

select <- dplyr::select
# Function to calculate presence absence data ----------------------------------------------------------------

pres.abs.rep <- function (start.date, end.date, sightings, dates){
  # Create a data frame with the detections
  sight <- data.frame (id = sightings, date = dates) %>%
    filter (date >= start.date, date <= end.date) %>%
    mutate (id = factor (id))
  
  # For each shark we'll start with the first detection only
  individuals <- levels (sight$id)
  presence.absence <- foreach (i=1:length (individuals), .combine = rbind, .inorder=FALSE) %dopar% {
    dates.tagged <- unique(sight$date)[unique(sight$date)>=sight$date[match (individuals[i], sight$id)]]
    dates.present <- sight$date[sight$id == individuals[i]]
    sight.shark <- foreach (j = 1:length (dates.present), .combine = rbind) %do% {
      dates.tagged.since <- dates.tagged[dates.tagged >= dates.present[j]]
      sight.shark.lag <- data.frame (date = dates.tagged.since, present = dates.tagged.since %in% dates.present, id = individuals[i]) %>%
        mutate (lag = as.numeric (date) - as.numeric (min (date)))
      return (sight.shark.lag)
    }
    return (sight.shark)
  } %>%
    mutate (day = yday (date), week = week (date), month = month(date))
  return (presence.absence)
}

# Acoustic presence absense --------------------------------------------------------------------

start.date = min (ACO.WS.DAILY.PRES$DATE)
end.date = max (ACO.WS.DAILY.PRES$DATE)
sightings = ACO.WS.DAILY.PRES$ECOCEAN
dates = ACO.WS.DAILY.PRES$DATE
PADet <- pres.abs.rep (start.date, end.date, sightings, dates) %>% 
  arrange (date)

# We will keep them separate for now because there are different explanatory variables 
# PADet will be explained by shark characteristics
# Remove first tagging of TZ003
WS.TAGS <- WS.TAGS[-4, ]
names (WS.TAGS) <- c ("number", "date.tagged", "id", "sex", "size", "transmitterID", "batch")
PADet <- merge (PADet, WS.TAGS[, 3:7], by =  "id")
# Add number of sharks tagged and Remove number of sharks known to have shed tags
TAG.LOST <- data.frame (id = c ("X","TZ-030", "TZ-068", "TZ-040", "TZ-032"), 
                        date = as.Date (c ("2012-10-26","2013-10-14", "2013-12-21", "2014-01-07", "2014-01-09")))
PADet <- ddply (PADet, "date", function (x, WS.TAGS, TAG.LOST){
  x$sharks.tagged <- max (which (x$date[1] >= as.Date(WS.TAGS$date))) 
  if (length (which (x$date[1] >= as.Date(TAG.LOST$date))) != 0) {
    x$sharks.tagged <- x$sharks.tagged - max (which (x$date[1] >= as.Date(TAG.LOST$date)))
  }
  return (x)
}, WS.TAGS = WS.TAGS, TAG.LOST = TAG.LOST)

# Calculate number of receivers working
ARRAY.EVENTS <- arrange (ARRAY.EVENTS, DATETIME)
times.in <- ddply (ARRAY.EVENTS, "STATIONNAME", function (x){
  deployed <- filter (x, EVENT == "DEP", lead (EVENT) == "RET")
  retrieved <- filter (x, EVENT == "RET", lag (EVENT) == "DEP")
  if (nrow (deployed) > 0) {
    y <- data.frame (station = first (x$STATIONNAME), 
                     date.in = deployed$DATETIME, 
                     date.out = retrieved$DATETIME,
                     rec.in = deployed$RECEIVERID,
                     rec.out = retrieved$RECEIVERID)
    return (y)
  } else return (NULL)
})
# Calculate array configuration and number of receivers working
PADet <- ddply (PADet, "date", function (x, times.in){
  stations.listening <- filter (times.in, x$date[1] >= as.Date (date.in), x$date[1] <= as.Date (date.out)) %>% 
    select (station) %>%
    unique()
  x$configuration <- do.call (paste, as.list(stations.listening$station))
  x$nStations <- nrow (stations.listening)
  return (x)
}, times.in = times.in)# %>%
#mutate (lag = as.numeric (date) - min (as.numeric (date)))

# # Calculate lags
# PADet <- ddply (PADet, "id", function (x){
#   x$lag <- as.numeric (x$date) - as.numeric (min(x$date)) + 1
#   return (x)
# })

# Delete data for sharks that were known to loose their tag
PADet <- ddply (PADet, "id", function (x, TAG.LOST){
  for (i in 1:length (levels (TAG.LOST$id))){
    if (as.character (TAG.LOST$id[i]) == as.character (first (x$id))){
      out <- filter (x, date < TAG.LOST$date[i])
      print (first (x$id))
    } else {out <- x}
  }
  return (out)
}, TAG.LOST = TAG.LOST)

# Lag Log
PADet <- dplyr::mutate (PADet, lagl = log (lag + 1))
# 
# # Cycle
# PADet.Cycle.1 <- filter (PADet, day <= 60)
# PADet.Cycle.1$week <- PADet.Cycle.1$week + 53
# PADet.Cycle.1$day <- PADet.Cycle.1$day + 366
# PADet.Cycle.1$month <- PADet.Cycle.1$month + 12
# 
# PADet.Cycle.2 <- filter (PADet, day >= 366 - 60)
# PADet.Cycle.2$week <- PADet.Cycle.2$week - 53
# PADet.Cycle.2$day <- PADet.Cycle.2$day - 366
# PADet.Cycle.2$month <- PADet.Cycle.2$month - 12
# 
# PADet.Cycle <- rbind (PADet.Cycle.2, PADet, PADet.Cycle.1)
# 
# 
# 
# # Create full model
# formulas <- vector ("list", 0)
# formulas[[1]] <- formula (present ~ s (day, bs = "cc") + s(lag, bs = "cr") + sex + size + nStations) 
# formulas[[2]] <- formula (present ~ s (day, bs = "cc") + s(lag, bs = "cr") + size + sharks.tagged + nStations) 
# formulas[[3]] <- formula (present ~ s (day, bs = "cc") + s(lag, bs = "cr") + size + nStations) 

args <- vector ("list", 0)
args[[1]] <- list (formula = formula (present ~ s (day, bs = "cc") + s (lag, bs = "cr") + sex + size + nStations), data = PADet, family = "binomial", gamma = 1.4)
args[[2]] <- list (formula = formula (present ~ s (day, bs = "cc") + s (lagl, bs = "cr") + sex + size  + nStations), data = PADet, family = "binomial", gamma = 1.4)
args[[3]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + nStations + s(lagl, bs = "cr")), data = PADet, family = "binomial", gamma = 1.4, random=list(id=~1))
args[[4]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + nStations + s(lagl, bs = "cr")), data = PADet, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id))
args[[5]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + nStations + s(lagl, bs = "cr")), data = PADet, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
args[[6]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + nStations + s(lagl, bs = "cr")), data = PADet, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id))
args[[7]] <- list (formula = formula (present ~ s (day, bs = "cc") + size + nStations + s(lagl, bs = "cr")), data = PADet, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id))
args[[8]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + nStations + s(lagl, bs = "cr")), data = PADet, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id))
args[[9]] <- list (formula = formula (present ~ s (day, bs = "cc") + nStations + s(lagl, bs = "cr")), data = PADet, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id))

registerDoMC (cores = 14)
md01 <- foreach (i=1:length (args)) %dopar% {
  do.call (gamm, args[[i]])
}
i <- 1
md00 <- foreach (i=1:2) %dopar% {
  do.call (gamm, args[[i]])
}

# All terms are significant for md04

md <- md04
r.md <- residuals (md$gam, type = "pearson")
f.md <- predict (md$gam, type = "response")
qplot (x = f.md, y = r.md, geom = c("smooth", "point")) 
qplot (PADet$sex, r.md, geom = "boxplot" )
qplot (as.factor(PADet$size), r.md,  geom = c("boxplot") )
qplot (PADet$sharks.tagged, y = r.md, geom = c("smooth", "point")) 
qplot (as.factor(PADet$nStations), y = r.md, geom = c("boxplot")) 

qplot (PADet$day, y = f.md, geom = c("point", "smooth")) 

plot(md$gam, all.terms = T)

pred.df <- expand.grid (date = 1:366, nStations = 13, lag = 1)
mutate (pred.df, predict (md$gam, pred.df, type = response)
        predict (md$gam, pred.df, type = "response")
        
        # Visual presence absense --------------------------------------------
        
        # Only tagged sharks & Create fake encounters for when there was surveys but not sharks
        ENCOUNTERS <- filter (ENCOUNTERS, ENCOUNTERS$ECOCEAN %in% WS.TAGS$ECOCEAN) %>%
          rbind (data.frame (DATE = SURVEYS$DATE, ID = NA,  PERIOD.DAY = NA, ECOCEAN = "XXX"))
        start.date = min (ACO.WS.DAILY.PRES$DATE)
        end.date = max (ACO.WS.DAILY.PRES$DATE)
        sightings = ENCOUNTERS$ECOCEAN
        dates = ENCOUNTERS$DATE
        PAEnc <- pres.abs (start.date, end.date, sightings, dates) %>%
          filter (id != "XXX")
        
        ggplot () + geom_smooth(data = PADet, aes (x = date, y = as.numeric (present)), method = "glm", family = "binomial", formula = y ~ ns(x,8)) + 
          geom_smooth(data = PAEnc, aes (x = date, y = as.numeric (present)), method = "glm", family = "binomial", formula = y ~ ns(x,8)) 
        
        ggplot () + geom_smooth(data = PADet, aes (x = day, y = as.numeric (present)), method = "glm", family = "binomial", formula = y ~ ns(x,8)) + 
          geom_smooth (data = PAEnc, aes (x = day, y = as.numeric (present)), method = "glm", family = "binomial", formula = y ~ ns(x,8))
        
        ggplot (PADet, aes (x = day, y = as.numeric (present))) + geom_smooth(method = "glm", family = "binomial", formula = y ~ ns(x,8))
        ggplot (filter (PADet, id == "TZ-995"), aes (x = week, y = as.numeric (present))) + geom_smooth(method = "glm", family = "binomial", formula = y ~ ns(x,8), aes (colour = id), level = 0)
        ggplot (PADet, aes (x = month, y = as.numeric (present))) + geom_smooth()
        
        