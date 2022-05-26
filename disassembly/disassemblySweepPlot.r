
parameterSweepPath <- getwd()
setwd("../")
source("figurePlotter.r")

knuc <- c(500)
khyd <- 0
fragment <- c('0.0001')
anneal <- c(1, 10)
preassemble <- c(20, 120)

disassemblySweepPlot(parameterSweepPath, knuc, khyd, fragment, anneal, preassemble)
