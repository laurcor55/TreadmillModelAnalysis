
parameterSweepPath <- getwd()
setwd("../")
source("figurePlotter.r")

knuc <- c(1000, 5000)
khyd <- 0.4
fragment <- '0.0001'
anneal <- c(1, 5)

velocitySweepPlot(parameterSweepPath, knuc, khyd, fragment, anneal)