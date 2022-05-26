
parameterSweepPath <- getwd()
setwd("../")
source("figurePlotter.r")

kswitch <- 50000
gtpase <- 0.3
ktongtp <- c(0.1, 1)
bottomon <- 10

topoff <- 6.5 #c(6.5, 10)
bottomoffcap <- c(0.3)
bottomoncap <- c(0.5, 1, 2, 5)
bottomoffgdp <- 6.5
mciz <- 2

gtpTurnoverSweepBottomCapPlot(parameterSweepPath, topoff, gtpase, kswitch, bottomon, mciz, ktongtp, bottomoffcap, bottomoncap, bottomoffgdp)
