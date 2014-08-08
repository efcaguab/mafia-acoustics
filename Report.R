## Graphs for the report

load (file = "../Processed Data/WSDetections.RData")
library (dplyr)
library (ggplot2)
library(splines)
library(MASS)

# Aggregate in daily bins
DET.WS <- mutate (DET.WS, DAY = cut (DATETIME, 'day'))
DET.WS <- group_by (DET.WS, DAY)
DET.WS.perday <- summarise (DET.WS, count = n())
DET.WS.perday$DAY <- as.Date (DET.WS.perday$DAY)

ggplot(DET.WS.perday, aes(x = DAY, y = count)) + geom_line() + geom_smooth(method = "glm", family="poisson", formula = y ~ ns(x, 6)) + theme_minimal() + xlab ("Date") + ylab ("Daily detections")

# Depth 

DET.WS$DEPTH <- (DET.WS$SENSOR1 - 2) *3
ggplot (DET.WS, aes(x = DATETIME, y = DEPTH)) + geom_jitter(alpha = 0.05, size = 1) + geom_smooth(method = "glm", family="poisson", formula = y ~ ns(x, 10), size = 2, colour = "black")  + theme_minimal () + xlab ("Date") + ylab ("Depth")
