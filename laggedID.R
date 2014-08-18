# Homemade lagged identification rates (conditional and absolute probability)

load (file = "../Processed Data/AcousticPresence.RData")
load (file = "../Processed Data/WSDetections.RData")
load (file = "../Processed Data/EncountersWeek.RData")

library (dplyr)
library (magrittr)
library (foreach)
library (doMC)
library(splines)
library(MASS)
library (mgcv)

registerDoMC (cores = 10)

start.date = min (ACO.WS.DAILY.PRES$DATE)
end.date = max (ACO.WS.DAILY.PRES$DATE)
sightings = ENCOUNTERS$ECOCEAN
dates = ENCOUNTERS$DATE

# Lagged Identification Rates Functions
LaggedIdentificationConditionalBinomial <- function (start.date, end.date, sightings, dates){
  # Create data frame and select only dates after the chosen initial date
  sight <- data.frame (id = sightings, date = dates) %>%
    filter (date >= start.date) %>%
    mutate (id = factor (id))
  individuals <- levels (sight$id)
  LIR <- foreach (i=1:length (individuals), .combine = rbind, .inorder = FALSE) %dopar% {
    # Find out dates in which the individual was identified (baseline)
    id.dates <- sight$date[which(sight$id == individuals[i])]
    # Cycle through each date
    prob.id <- foreach (j=1:length (id.dates), .combine = rbind, .inorder = FALSE) %do%{
      # Select only dates in the future, see if the individual resighted is the chosen one or not and calculate the lag
      sight.date <- filter (sight, date >= id.dates[j]) %>%
        mutate (present = id == individuals[i], lag = as.numeric (date - min (date))) %>%
        group_by (lag)
      # For each lag count the number of successes and failures, and add the ID
      prob <- summarise (sight.date, s = sum(present), f = n() - s) %>%
        mutate (id = individuals[i])
      return (prob)
    }
  return (prob.id)
  }
  return (LIR)
}

LaggedIdentificationsAbsolute <- function (start.date, end.date, sightings, dates){
  sightings <- sightings[dates >= start.date]
  dates <- dates[dates >= start.date]
  FRAME <- expand.grid (DATE = seq (start.date, end.date, 'day'), 
                        ECOCEAN = unique (sightings), PRESENT = FALSE)
  
  # Fill data frame to see if it was present or not
  for (i in 1:length (levels (FRAME$ECOCEAN))){
    FRAME$PRESENT[!is.na (match( FRAME$DATE, dates[sightings == levels(FRAME$ECOCEAN)[i]])) & FRAME$ECOCEAN == levels (FRAME$ECOCEAN)[i]] <- TRUE
  }
  
#  LAGGED.ID  <- data.frame (ECOCEAN = NA, LAG = NA, PRESENT = NA)
  LAGGED.ID <- foreach (i=1:length (levels (FRAME$ECOCEAN)), .combine = rbind) %dopar% {
    FRAME.SHARK <- FRAME[FRAME$ECOCEAN == levels (FRAME$ECOCEAN)[i], ]
    
    # Delete previous NO-OBSERVATIONS
    if (match(TRUE, FRAME.SHARK$PRESENT) != 1){FRAME.SHARK <- tail (FRAME.SHARK, n = - match(TRUE, FRAME.SHARK$PRESENT) + 1)}
    
    days.present <- which (FRAME.SHARK$PRESENT)
   #  FRAME.SHARK.LAG <- data.frame (ECOCEAN = NA, LAG = NA, PRESENT = NA)
    FRAME.SHARK.LAG <- foreach (j=1:sum(FRAME.SHARK$PRESENT), .combine = rbind) %do% {
      #message ("i =", i, "   j = ", j)
      #message ("Length FRAME SHARK = ", length(FRAME.SHARK$PRESENT), "    days.present[j] = ", days.present[j])
      FRAME.SHARK.AUX <- data.frame (ECOCEAN = levels (FRAME$ECOCEAN)[i], 
                                     LAG = seq(1, length(FRAME.SHARK$PRESENT) - days.present[j] + 1),
                                     PRESENT = FALSE)
      if (j == 1){FRAME.SHARK.AUX$PRESENT[days.present-days.present[j]+1] <- TRUE}
      else {FRAME.SHARK.AUX$PRESENT[tail(days.present-days.present[j]+1, n = -j +1)] <- TRUE}
      return (FRAME.SHARK.AUX)
    }
   return (FRAME.SHARK.LAG)
  }
  
  LAGGED.ID <- LAGGED.ID[!is.na (LAGGED.ID$LAG) & LAGGED.ID$LAG != 1, ]
  LAGGED.ID$LAG <- LAGGED.ID$LAG -1
  
  return (LAGGED.ID)
  
}



detLaggedBin <- LaggedIdentificationConditionalBinomial (start.date = min (ACO.WS.DAILY.PRES$DATE),
                                                         end.date = max (ACO.WS.DAILY.PRES$DATE),
                                                         sightings = ACO.WS.DAILY.PRES$ECOCEAN,
                                                         dates = ACO.WS.DAILY.PRES$DATE) %>%
  mutate (p = s / (s+f))

detLaggedAbs <- LaggedIdentificationsAbsolute (start.date = min (ACO.WS.DAILY.PRES$DATE),
                                               end.date = max (ACO.WS.DAILY.PRES$DATE),
                                               sightings = ACO.WS.DAILY.PRES$ECOCEAN,
                                               dates = ACO.WS.DAILY.PRES$DATE)
detLaggedAbs <- mutate (detLaggedAbs, p = as.numeric(PRESENT))
names (detLaggedAbs) <- c ("id", "lag", "present", "p")

visLaggedBin <- LaggedIdentificationConditionalBinomial(start.date = min (ACO.WS.DAILY.PRES$DATE),
                                                        end.date = max (ACO.WS.DAILY.PRES$DATE),
                                                        sightings = ENCOUNTERS$ECOCEAN, 
                                                        dates = ENCOUNTERS$DATE) %>%
  mutate (p = s / (s+f))

visLaggedAbs <- LaggedIdentificationsAbsolute (start.date = min (ACO.WS.DAILY.PRES$DATE),
                                               end.date = max (ACO.WS.DAILY.PRES$DATE),
                                               sightings = ENCOUNTERS$ECOCEAN, 
                                               dates = ENCOUNTERS$DATE) %>%
  mutate (p = as.numeric(PRESENT))
names (visLaggedAbs) <- c ("id", "lag", "present", "p")

MDA1 <- gam (present ~ s(lag), family = binomial, data = detLaggedAbs )

# I came to the conclussion that the lagged identification rate is not really a good idea

ggplot(detLaggedBin, aes(x = lag, y = p)) + geom_smooth (method = "glm", family = "binomial", formula = y ~ ns(x,7), aes (colour = id), level = 0)
ggplot(detLaggedAbs, aes(x = lag, y = p)) + geom_smooth (method = "glm", family = "binomial", formula = y ~ ns(x,7), aes (colour = id), level = 0)

ggplot() + geom_smooth (data = detLaggedBin, aes(x = lag, y = p), method = "glm", family = "binomial", formula = y ~ ns(x,7), colour = "black") +
  geom_smooth (data = visLaggedBin, aes(x = lag, y = p), method = "glm", family = "binomial", formula = y ~ ns(x,7)) 

ggplot() + geom_smooth (data = detLaggedAbs, aes(x = lag, y = p), method = "glm", family = "binomial", formula = y ~ ns(x,7), colour = "black") + 
  geom_smooth (data = visLaggedAbs, aes(x = lag, y = p), method = "glm", family = "binomial", formula = y ~ ns(x,7))
