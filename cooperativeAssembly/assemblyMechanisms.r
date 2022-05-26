library(matrixStats)
library(vioplot)
library(scales)

assemblyGdpexchangeSweep <- function(){
  fileStart <- '/home/lauren/Documents/research/treadmilling-model/figures/assemblySweep/results/gdpexchange'
  fileEnd <- '.csv'
  gdpexchange <- c(0.5, 1, 2)
  fullColors <- c()
  plot(0, 0, xlim=c(0,15),  ylim=c(0, 4), col="white", ylab='', xlab='')
  for (ii in c(1:length(gdpexchange))){
    fileName <- paste0(fileStart, gdpexchange[ii], fileEnd)
    results <- read.csv(fileName, header=FALSE, sep=',', dec='.', stringsAsFactors=FALSE)
    results <- as.matrix(results)
    time <- results[,1]
    monomerConc <- results[,2]
    color <- colorPalette(ii, 1)
    fullColors <- c(fullColors, color)
    lines(time, monomerConc, col=color, lty=1, lwd=2)
  }
  fileNameExpected <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005a/chen2005_fig3.csv'
  assemblyPlotExperimental(fileNameExpected, 3.12)
  label <- c()
  for (ii in c(1:length(gdpexchange))){
    label <- c(label,  as.expression(bquote(.(gdpexchange[ii])~s^-1)))
  }
  grid()
  legend('topright', legend=label, title=expression('k'['GDP exchange']), col=fullColors, lty=1, lwd=2, cex=0.75)
  mtext(expression('Monomeric FtsZ ('*mu*'M)'), side=2, line=2)
  mtext('Time (s)', side=1, line=2)
}


assemblyKswitchSweep <- function(){
  fileStart <- '/home/lauren/Documents/research/treadmilling-model/figures/assemblySweep/results/knuc'
  fileEnd <- '.csv'
  knuc <- c(2500, 5000, 10000)
  fullColors <- c()
  plot(0, 0, xlim=c(0,15),  ylim=c(0, 4), col="white", ylab='', xlab='')
  for (ii in c(1:length(knuc))){
    fileName <- paste0(fileStart, knuc[ii], fileEnd)
    results <- read.csv(fileName, header=FALSE, sep=',', dec='.', stringsAsFactors=FALSE)
    results <- as.matrix(results)
    time <- results[,1]
    monomerConc <- results[,2]
    color <- colorPalette(ii, 1)
    fullColors <- c(fullColors, color)
    lines(time, monomerConc, col=color, lty=1, lwd=2)
  }
  fileNameExpected <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005a/chen2005_fig3.csv'
  assemblyPlotExperimental(fileNameExpected, 3.12)
  label <- knuc
  grid()
  legend('topright', legend=label, title=expression('K'['nuc']), col=fullColors, lty=1, lwd=3, cex=0.75)
  mtext(expression('Monomeric FtsZ ('*mu*'M)'), side=2, line=2)
  mtext('Time (s)', side=1, line=2)

}




subplot(c(0, 0.25, 0, 1))
assemblyKswitchSweep()
addLabel('a')

subplot(c(0.25, 0.5, 0, 1))
assemblyGdpexchangeSweep()
addLabel('b')




source('/home/lauren/Documents/research/treadmilling-model/figures/figurePlotter.r')


subplot(c(0.5, 0.75, 0, 1))
parameterSweepPath <- '/home/lauren/Documents/research/treadmilling-model/figures/assemblySweep'
modelFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/assemblySweep/results/threeRuns.csv'
assemblyMultiPlot(parameterSweepPath, modelFileName)
addLabel('c')




subplot(c(0.75, 1, 0, 1))
startingFtsz <- c(1.56, 2.42, 3.12)
fileName <- '/home/lauren/Documents/research/treadmilling-model/figures/assemblySweep/results/ecF268C.csv'
assemblyPlot(fileName, startingFtsz, TRUE)
fileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005a/chen2005_fig3.csv'
assemblyPlotExperimental(fileName, startingFtsz)
addLabel('d')



