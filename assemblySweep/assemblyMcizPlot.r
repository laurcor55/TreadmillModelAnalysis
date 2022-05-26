
setwd("../")
source("figurePlotter.r")

kswitch <- c(15000)
gtpase <- c(0.5)
topoff <- 7.5
gdpexchange <- c(0.5)
caponpf <- c(1)
maxTime <- 30
bottomon <- 10
bottomoff <- 1

experimentalFolder <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/bisson2015/figS5'
parameterSweepPath <- '/home/lauren/Documents/research/treadmilling-model/figures/assemblySweep'

startFtsz <- c(4, 5, 6)
mciz <- c(0, 1)
label <- c('a', 'b', 'c', 'd')

rows <- length(caponpf)
cols <- length(mciz)

resultsFile <- paste0(parameterSweepPath, '/mcizAssembly.png')
png(resultsFile, width=200*cols, height=200*rows)
margin_top <- 0.25/rows
margin_left <- 0.25/cols

division <- 1/length(mciz)
plot.new()

scale <- 1.5
domain <- c(0, 30)
range <- c(0, 8)

maxConc <- max(startFtsz)+1 + max(mciz)
for (ii in c(1:length(mciz))){
  for (jj in c(1:length(caponpf))){
    margins <- c(division*(ii-1), division*ii, 0, 1)
    
    startingFtsz = startFtsz + mciz[ii]

    sweepSubplots(margin_top, margin_left, jj, ii, rows, cols, domain, range)

    experimentalFileName <- paste0(experimentalFolder, label[ii], '_conc.csv')
    resultsExp <- assemblyExtractExperimental(experimentalFileName)
    assemblyPlotExperimental(experimentalFileName, startingFtsz)
    modelFileName = paste0(parameterSweepPath, '/results/kswitch', kswitch,'_gtpase', gtpase, '_topoff', topoff,'_caponpf', caponpf[jj], '_gdpexchange', gdpexchange, '_bottomon', bottomon, '_bottomoff', bottomoff, '_mciz', mciz[ii], '.csv') #
    assemblyPlot(modelFileName, startingFtsz, FALSE)
  }
}
sweep2dLabels(margin_top, margin_left, caponpf, mciz, "kb+Cap", "MciZ", "Time (s)", expression(Free~Monomer~(mu*M)))  


dev.off()