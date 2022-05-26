
setwd("../")
source("figurePlotter.r")

kswitch <- c(30000)
gtpase <- c(0.4)
topoff <- 6.5
bottomon <- c(10)
bottomoff <- 5
gdpexchange <- c(1)
caponpf <- c(2)

#experimentalFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005a/chen2005_fig3.csv'
#startFtsz <- c(3.12, 2.42, 1.56)
#maxTime <- 15
#mciz <- c(0)


#experimentalFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005b/hmk/fig6c_conc.csv'
#startFtsz <- c(1.1, 2)

#experimentalFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005b/hek/fig6a_conc.csv'
#startFtsz <- c(1.6, 2.6, 3)

#startFtsz <- c(3, 2.2, 1.7)
#experimentalFileName <- 'experimental/chen2005b/mek/dataTransformed.csv'
#mciz <- 0
#maxTime <- 25


experimentalFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/bisson2015/figS5a_conc.csv'
startFtsz <- c(4, 5, 6)
mciz <- c(0)
maxTime <- 30

#experimentalFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/bisson2015/figS5b_conc.csv'
#startFtsz <- c(5, 6, 7)
#mciz <- 1
#maxTime <- 30

parameterSweepPath <- '/home/lauren/Documents/research/treadmilling-model/figures/assemblySweep'
assemblySweepPlot(parameterSweepPath, experimentalFileName, startFtsz, kswitch, gdpexchange, gtpase, topoff, caponpf, bottomon, bottomoff, mciz, maxTime)
