
parameterSweepPath <- getwd()
setwd("../")
source("figurePlotter.r")

knuc <- c(5000, 10000)
khyd <- 0.5
fragment <- '0.0001'
anneal <- c(0.1, 0.5)
mciz <- c(0, 2)

gtpTurnoverSweepPlot(parameterSweepPath, knuc, khyd, fragment, anneal, mciz)
