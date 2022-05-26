library(beeswarm)
library(dplyr)




parameterSweepPath <- getwd()
setwd("../")
source("figurePlotter.r")



plotResult <- function(file_name, experimentalFileName, label){
  df <- read.csv(file_name, header=FALSE, sep=',', dec='.', stringsAsFactors=FALSE)
  df <- data.matrix(df)
  maxX <- max(df)
  maxY <- max(df[2,])*1.1
  scale <- 1.5
  plot(0, 0, col="white", xlim=c(0,maxX), ylim=c(0,maxY), xlab="Time (s)", ylab="Mixed PF", cex.axis=scale, cex.lab=scale)
  time <- df[1,]
  for (ii in c(2:nrow(df))){
    lines(time, df[ii,], col=colorPalette(ii-1, 1))
  }
  dfExp <- read.csv(experimentalFileName, header=FALSE, sep=',', dec='.', stringsAsFactors=FALSE)
  dfExp <- data.matrix(dfExp)
  time <- dfExp[,1]
  fluor <- dfExp[,2]
  fluor <- -1*fluor + 1
  lines(time, fluor, col='black', lwd=2, lty=2)
  text(maxX*0.05, maxY*0.95, label, cex=scale)
}

outputFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/pfMixing/pfMixing.png'
png(outputFileName, width=400, height=800)
par(fig=c(0, 1, 0.5, 1), mar=c(5, 5, 2, 2))
file_name <- '/home/lauren/Documents/research/treadmilling-model/figures/pfMixing/results/kswitch10000_gtpase0.2_topoff7.5_bottomoff1_gdpexchange2_bottomon10.csv'
experimentalFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005a/chen2005_fig4a.csv'
plotResult(file_name, experimentalFileName, 'a)')
par(fig=c(0, 1, 0, 0.5), new=TRUE)
file_name <- '/home/lauren/Documents/research/treadmilling-model/figures/pfMixing/results/kswitch5000_gtpase0_topoff7.5_bottomoff1_gdpexchange2_bottomon10_topongtp0.5.csv'
experimentalFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005a/chen2005_fig4b.csv'
plotResult(file_name, experimentalFileName, 'b)')
dev.off()


