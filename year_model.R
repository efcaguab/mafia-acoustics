
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
library (parallel)
library (reshape2)
library (plot3D)
library (mail)
library (boot)
registerDoMC (cores = 14)
cl <- makeCluster (14)      
select <- dplyr::select

# Functions to calculate presence absence data ----------------------------------------------------------------

pres.abs <- function (start.date, end.date, sightings, dates){
  # Create a data frame with the detections
  sight <- data.frame (id = sightings, date = dates) %>%
    filter (date >= start.date, date <= end.date) %>%
    mutate (id = factor (id))
  
  # For each shark we'll start with the first detection only
  individuals <- levels (sight$id)
  presence.absence <- foreach (i=1:length (individuals), .combine = rbind) %do% {
    dates.tagged <- unique(sight$date)[unique(sight$date)>sight$date[match (individuals[i], sight$id)]]
    dates.present <- sight$date[sight$id == individuals[i]]
    sight.shark <- data.frame (date = dates.tagged, 
                               present = dates.tagged %in% dates.present, 
                               id = individuals[i]) %>%
      tbl_df () %>% 
      mutate (lag = as.numeric (date) - as.numeric (min(date)) + 1)
  } %>%
    mutate (day = yday (date), week = week (date), month = month(date))
  return (presence.absence)
}


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

pres.abs.lag <- function (start.date, end.date, sightings, dates){
  # Create a data frame with the detections
  sight <- data.frame (id = sightings, date = dates) %>%
    filter (date >= start.date, date <= end.date) %>%
    mutate (id = factor (id)) %>% 
    arrange (date)
  
  # For each shark we'll start with the first detection only
  individuals <- levels (sight$id)
  presence.absence <- foreach (i=1:length (individuals), .combine = rbind) %dopar% {
    dates.present <- sight$date[sight$id == individuals[i]] %>%
      as.numeric ()
    dates.tagged <- unique(sight$date)[unique(sight$date)>sight$date[match (individuals[i], sight$id)]]
    dates.comb <- as.data.frame (t (combn (dates.tagged, 2))) %>%
      tbl_df()
    names (dates.comb) <- c ("date.1", "date.2")
    dates.comb <- mutate (dates.comb, lag = date.2 - date.1, 
                          present = (date.1 %in% dates.present) & (date.2 %in% dates.present), 
                          date = as.Date(date.1, origin = "1970-01-01"), 
                          id = individuals[i]) %>%
      select (-date.1, -date.2)
    return (dates.comb)
  } %>%
    mutate (day = yday (date), week = week (date), month = month(date))
  return (presence.absence)
}


# Functions to populate data frame ----------------------------------------

add.shark.info <- function (PADet, WS.TAGS) {
  # Remove first tagging of TZ003
  WS.TAGS <- WS.TAGS[-4, ]
  names (WS.TAGS) <- c ("number", "date.tagged", "id", "sex", "size", "transmitterID", "batch")
  PADet <- merge (PADet, WS.TAGS[, 3:7], by =  "id")
  return (PADet)
}

acoustic.exposure <- function (PADet, WS.TAGS, TAG.LOST, ARRAY.EVENTS){
  # Remove first tagging of TZ003
  WS.TAGS <- WS.TAGS[-4, ]
  names (WS.TAGS) <- c ("number", "date.tagged", "id", "sex", "size", "transmitterID", "batch")
  
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
  }, times.in = times.in)
  
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
  
  return (PADet)
}

add.survey.effort <- function (PAEnc, SURVEYS) {
  # In the surveys file there are duplicated record for 2013-01-30
  SURVEYS <- filter (SURVEYS, ! duplicated (DATE))
  SURVEYS <- mutate (SURVEYS, SNORKELERS = as.numeric ( as.character (SURVEYS$SNORKELERS)), 
                     DURATION = abs (as.numeric (DURATION)))
  
  surveys.complete <- data.frame (DATE = unique (filter (PAEnc, ! (date %in% SURVEYS$DATE))$date), 
                                  START.TIME = NA, END.TIME = NA, 
                                  SNORKELERS = mean (SURVEYS$SNORKELERS, na.rm = T), 
                                  DURATION = as.vector (mean (SURVEYS$DURATION))) %>%
    rbind (SURVEYS)
  
  names (surveys.complete)[c (1, 5)] <- c("date", "duration") 
  PAEnc <- merge (PAEnc, select (surveys.complete, date, duration))
  
  return (PAEnc)
}

# Acoustic presence absense --------------------------------------------------------------------

start.date = min (ACO.WS.DAILY.PRES$DATE)
end.date = max (ACO.WS.DAILY.PRES$DATE)
sightings = ACO.WS.DAILY.PRES$ECOCEAN
dates = ACO.WS.DAILY.PRES$DATE

TAG.LOST <- data.frame (id = c ("X","TZ-030", "TZ-068", "TZ-040", "TZ-032"), 
                        date = as.Date (c ("2012-10-26","2013-10-14", "2013-12-21", "2014-01-07", "2014-01-09")))

PADet.one <- pres.abs (start.date, end.date, sightings, dates) %>% 
  arrange (date) %>%
  add.shark.info (WS.TAGS) %>%
  acoustic.exposure (WS.TAGS, TAG.LOST, ARRAY.EVENTS) %>%
  mutate (lagl = log (lag + 1))

PADet.comb <- pres.abs.lag (start.date, end.date, sightings, dates) %>% 
  arrange (date) %>%
  add.shark.info (WS.TAGS) %>%
  acoustic.exposure (WS.TAGS, TAG.LOST, ARRAY.EVENTS) %>%
  mutate (lagl = log (lag + 1))

# Reduce to weekly averages to work easier


# Acoustic models ---------------------------------------------------------

# For now just using the single combinatory

args <- vector ("list", 0)
args[[1]] <- list (formula = formula (present ~ s (day, bs = "cc") + s (lag, bs = "cr") + sex + size + nStations), data = PADet.one, family = "binomial", gamma = 1.4)
args[[2]] <- list (formula = formula (present ~ s (day, bs = "cc") + s (lagl, bs = "cr") + sex + size  + nStations), data = PADet.one, family = "binomial", gamma = 1.4)
args[[3]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + nStations + s(lag, bs = "cr")), data = PADet.one, family = "binomial", gamma = 1.4, random=list(id=~1))
args[[4]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + nStations + s(lag, bs = "cr")), data = PADet.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id))
args[[5]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + nStations + s(lag, bs = "cr")), data = PADet.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
args[[6]] <- list (formula = formula (present ~ s (day, bs = "cc") + size + nStations + s(lag, bs = "cr")), data = PADet.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id))
args[[7]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + nStations + s(lag, bs = "cr")), data = PADet.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id))
args[[8]] <- list (formula = formula (present ~ s (day, bs = "cc") + nStations + s(lag, bs = "cr")), data = PADet.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id))
args[[9]] <- list (formula = formula (present ~ s (day, bs = "cc") + size + nStations + s(lag, bs = "cr")), data = PADet.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
args[[10]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + nStations + s(lag, bs = "cr")), data = PADet.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
args[[11]] <- list (formula = formula (present ~ s (day, bs = "cc") + nStations + s(lag, bs = "cr")), data = PADet.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))

registerDoMC (cores = 14)
md01 <- foreach (i=1:length (args)) %dopar% {
  do.call (gamm, args[[i]])
}

# All terms are significant for md04
# args[[1]] <- list (formula = formula (present ~ s (day, bs = "cc") + nStations + lagl), data = PADet.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
# md02 <- do.call (gamm, args[[1]])

# save (md01, file = "../../../../../MAFIA_md01_one.RData")
# load (file = "../../../../../MAFIA_md01_one.RData")

md <- md01[[11]]
r.md <- residuals (md$gam, type = "pearson")
f.md <- predict (md$gam, type = "response")
qplot (x = f.md, y = r.md, geom = c("smooth", "point")) 
qplot (PADet.one$sex, r.md, geom = "boxplot" )
qplot (as.factor(PADet.one$size), r.md,  geom = c("boxplot") )
qplot (PADet$sharks.tagged, y = r.md, geom = c("smooth", "point")) 
qplot (as.factor(PADet.one$nStations), y = r.md, geom = c("boxplot")) 

qplot (PADet.one$day, y = f.md, geom = c("point", "smooth")) 

plot(md$gam, all.terms = T)


pred.df <- expand.grid (day = 1:366, nStations = 13, lag = log (c (1:366)))
pred.df <- mutate (pred.df, prediction = predict(md$gam, pred.df, type= "response", se.fit= T)$fit, 
          se = predict(md$gam, pred.df, type= "response", se.fit= T)$se.fit) %>% tbl_df()

ggplot(filter (pred.df, lagl %in% log (c (1, 8, 30,  90,  366))), aes (x = as.Date("2013-01-01") + day)) + geom_line (aes(y = prediction, colour = as.factor (lagl))) + geom_ribbon (aes (ymax = prediction + se, ymin = prediction - se, fill = as.factor(lagl)), alpha = 0.2) + theme_classic()

# pred.a <- acast (pred.df, day ~ lagl, value.var= "prediction")
# persp3D (pred.a, x = 1:366, y = 1:366, xlab = "Time of the year", ylab = "Lag", zlab = "Sighting probability", contour =  list(col = "black", side = c("z")), theta = -135-90, phi = 35, alpha = 0.5)

# i.date <- as.Date ("2012-10-01")
# pred.ot <- data.frame (day = yday(seq (i.date, i.date+ 365, by  ='day')), 
#                        nStations = 13, lagl = log (c (1:(366*3))), lag = c (1:(366*3)))
# pred.ot <- mutate (pred.ot, prediction = predict(md$gam, pred.ot, type= "response", se.fit= T)$fit, 
#                    se = predict(md$gam, pred.ot, type= "response", se.fit= T)$se.fit) %>% tbl_df()
# 
# ggplot (pred.ot, aes (x = i.date + lag, y = prediction)) + geom_line() + 
#   geom_ribbon (aes (ymax = prediction + se*1.92, ymin = prediction - se*1.92), alpha = 0.2, linetype = 2) + theme_classic()

# Visual presence absense --------------------------------------------

# Only tagged sharks & Create fake encounters for when there was surveys but not sharks
ENCOUNTERS <- filter (ENCOUNTERS, ENCOUNTERS$ECOCEAN %in% WS.TAGS$ECOCEAN) %>%
  rbind (data.frame (DATE = SURVEYS$DATE, ID = NA,  PERIOD.DAY = NA, ECOCEAN = "XXX"))
# All tagged sharks & Create fake encounters for when there was surveys but not sharks
# ENCOUNTERS <- filter (ENCOUNTERS, T) %>%
#   rbind (data.frame (DATE = SURVEYS$DATE, ID = NA,  PERIOD.DAY = NA, ECOCEAN = "XXX"))
start.date <- min (ACO.WS.DAILY.PRES$DATE)
end.date <- max (ACO.WS.DAILY.PRES$DATE)
sightings <- ENCOUNTERS$ECOCEAN
dates <- ENCOUNTERS$DATE

PAEnc.one <- pres.abs (start.date, end.date, sightings, dates) %>%
  filter (id != "XXX") %>%
  add.shark.info (WS.TAGS) %>%
  mutate (lagl = log (lag + 1)) %>% 
  add.survey.effort (SURVEYS) %>%
  tbl_df()

PAEnc.comb <- pres.abs.lag (start.date, end.date, sightings, dates) %>% 
  filter (id != "XXX") %>%
  add.shark.info (WS.TAGS) %>%
  mutate (lagl = log (lag + 1)) %>%
  add.survey.effort (SURVEYS) %>%
  tbl_df() 
 
# Find out which dates are not in the surveys table


# args <- vector ("list", 0)
# args[[1]] <- list (formula = formula (present ~ s (day, bs = "cc") + s (lag, bs = "cr") + sex + size + duration), data = PAEnc.one, family = "binomial", gamma = 1.4)
# args[[2]] <- list (formula = formula (present ~ s (day, bs = "cc") + s (lag, bs = "cr") + sex + size + duration), data = PAEnc.one, family = "binomial", gamma = 1.4)
# args[[3]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + duration + s(lag, bs = "cr")), data = PAEnc.one, family = "binomial", gamma = 1.4, random=list(id=~1))
# args[[4]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + duration + s(lag, bs = "cr")), data = PAEnc.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id))
# args[[5]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + duration + s(lag, bs = "cr")), data = PAEnc.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
# args[[6]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + duration + s(lag, bs = "cr")), data = PAEnc.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
# args[[7]] <- list (formula = formula (present ~ s (day, bs = "cc") + size + duration + s(lag, bs = "cr")), data = PAEnc.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
# args[[8]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + size + s(lag, bs = "cr")), data = PAEnc.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
# args[[9]] <- list (formula = formula (present ~ s (day, bs = "cc") + sex + s(lag, bs = "cr")), data = PAEnc.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
# args[[10]] <- list (formula = formula (present ~ s (day, bs = "cc") + duration + s(lag, bs = "cr")), data = PAEnc.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
# args[[11]] <- list (formula = formula (present ~ s (day, bs = "cc") + size + s(lag, bs = "cr")), data = PAEnc.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
# args[[12]] <- list (formula = formula (present ~ s (day, bs = "cc") + s(lag, bs = "cr")), data = PAEnc.one, family = "binomial", gamma = 1.4, correlation = corAR1(form = ~ lag | id), random=list(id=~1))
# 
# registerDoMC (cores = 14)
# me01 <- foreach (i=1:length (args)) %dopar% {
#   do.call (gamm, args[[i]])
# }

# save (me01, file = "../../../../../MAFIA_me01_one.RData")
load (file = "../../../../../MAFIA_me01_one.RData")
me <- me01[[12]]
pred.df <- expand.grid (day = 1:366, lag = c (1:366)) 
pred.df <- mutate (pred.df, prediction = predict(me$gam, pred.df, type= "response", se.fit= T)$fit, 
                   se = predict(me$gam, pred.df, type= "response", se.fit= T)$se.fit) %>% tbl_df()

ggplot(filter (pred.df, lag %in% c (1, 8, 30,  90,  366)), aes (x = as.Date("2013-01-01") + day)) + geom_line (aes(y = prediction, colour = as.factor (lag))) + geom_ribbon (aes (ymax = prediction + se, ymin = prediction - se, fill = as.factor(lag)), alpha = 0.2) + theme_classic()

i.date <- as.Date ("2012-10-01")
pred.ot.e <- data.frame (day = yday(seq (i.date, i.date+ 365, by  ='day')), 
                       lagl = log (c (1:(366*2))), lag = c (1:(366*2)))
pred.ot.e <- mutate (pred.ot.e, prediction = predict(me$gam, pred.ot.e, type= "response", se.fit= T)$fit, 
                   se = predict(me$gam, pred.ot.e, type= "response", se.fit= T)$se.fit ) %>% tbl_df()

pred.ot.d <- data.frame (day = yday(seq (i.date, i.date+ 365, by  ='day')), 
                       nStations = 13, lagl = log (c (1:(366*2))), lag = c (1:(366*2)))
pred.ot.d <- mutate (pred.ot.d, prediction = predict(md$gam, pred.ot.d, type= "response", se.fit= T)$fit, 
                   se = predict(md$gam, pred.ot.d, type= "response", se.fit= T)$se.fit) %>% tbl_df()

ggplot () + geom_line(data = pred.ot.e, aes (x = i.date + lag, y = prediction)) + 
  geom_ribbon (data = pred.ot.e, aes (x = i.date + lag, ymax = prediction + se, ymin = prediction - se), alpha = 0.2, linetype = 2) +
  geom_line(data = pred.ot.d, aes (x = i.date + lag, y = prediction)) + 
  geom_ribbon (data = pred.ot.d, aes (x = i.date + lag, ymax = prediction + se, ymin = prediction - se), alpha = 0.2, linetype = 2) + theme_classic() + coord_cartesian (ylim = c(-0,1))

## BY WEEK DAMN SHIT -----------------------------------------------------------------------
 
 


# ggplot () + geom_smooth(data = PADet, aes (x = date, y = as.numeric (present)), method = "glm", family = "binomial", formula = y ~ ns(x,8)) + 
#   geom_smooth(data = PAEnc, aes (x = date, y = as.numeric (present)), method = "glm", family = "binomial", formula = y ~ ns(x,8)) 
# 
# ggplot () + geom_smooth(data = PADet, aes (x = day, y = as.numeric (present)), method = "glm", family = "binomial", formula = y ~ ns(x,8)) + 
#   geom_smooth (data = PAEnc, aes (x = day, y = as.numeric (present)), method = "glm", family = "binomial", formula = y ~ ns(x,8))
# 
# ggplot (PADet, aes (x = day, y = as.numeric (present))) + geom_smooth(method = "glm", family = "binomial", formula = y ~ ns(x,8))
# ggplot (filter (PADet, id == "TZ-995"), aes (x = week, y = as.numeric (present))) + geom_smooth(method = "glm", family = "binomial", formula = y ~ ns(x,8), aes (colour = id), level = 0)
# ggplot (PADet, aes (x = month, y = as.numeric (present))) + geom_smooth()
# 

        