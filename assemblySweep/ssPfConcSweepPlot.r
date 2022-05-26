
setwd("../")
source("figurePlotter.r")

kswitch <- c(5000)
gtpase <- c(0.3)
topoff <- 7.5
bottomon <- 10
bottomoff <- 3
gdpexchange <- 0.5
caponpf <- 10
startFtsz <- c(2:8)

mciz <- 2

parameterSweepPath <- '/home/lauren/Documents/research/treadmilling-model/figures/assemblySweep'
ssPfConcSweepPlot(parameterSweepPath, experimentalFileName, startFtsz, kswitch, gdpexchange, gtpase, topoff, caponpf, bottomon, bottomoff, mciz)
