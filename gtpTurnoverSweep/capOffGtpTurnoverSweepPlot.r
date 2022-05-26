
parameterSweepPath <- getwd()
setwd("../")
source("figurePlotter.r")


capOffSweepPlot <- function(parameterSweepPath, topoff, gtpase, kswitch, gdpexchange, bottomoff, mciz, ktongtp, bottomcapoff){
  resultsFile <- paste0(parameterSweepPath, '/capOffGtpTurnoverSweep.png')
  png(resultsFile, width=250, height=250)
  
  plot.new()
  plot(0,0, xlim=c(0,8), ylim=c(0,30))
  for (jj in c(1:length(bottomoffcap))){
    color <- colorPalette(jj, 1)
    file_name = paste0(parameterSweepPath, '/results/kswitch', kswitch,'_gtpase', gtpase, '_topoff', topoff, '_bottomoff', bottomoff, '_gdpexchange', gdpexchange, '_ktongtp', ktongtp, '_bottomoffcap', bottomoffcap[jj], '_mciz', mciz,'.csv')
    gtpTurnoverPlot(file_name, criticalConc, gtpTurnover, color, FALSE, FALSE)
  }
  dev.off()
}

kswitch <- c(10000)
gtpase <- c(0.3)
ktongtp <- c(1)
topoff <- 7.5
gdpexchange <- 0.5
bottomoff <- 1
mciz <- c(2)
criticalConc <- 0
gtpTurnover <- 0
bottomoffcap <- c(1, 0.1, 0.01, 0.001)


capOffSweepPlot(parameterSweepPath, topoff, gtpase, kswitch, gdpexchange, bottomoff, mciz, ktongtp, bottomcapoff)
