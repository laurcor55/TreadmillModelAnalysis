
parameterSweepPath <- getwd()
setwd("../")
source("figurePlotter.r")

knuc <- 1000
khyd <- 0
fragment <- c('0.0001', '1e-05')
anneal <- c(0.1, 1, 10)

pfMixingSweepPlot(parameterSweepPath, knuc, khyd, fragment, anneal)
