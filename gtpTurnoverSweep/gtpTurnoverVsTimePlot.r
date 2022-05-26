
setwd("../")
source("figurePlotter.r")
library(pracma)


gtpTurnoverVsTimePlot <- function(scale, ratio, gtpase, bottomoffgdp){
  parameterSweepPath <- '/home/lauren/Documents/research/treadmilling-model/figures/gtpTurnoverSweep'
  palette <- c()
  
  for (ii in c(1:length(bottomoffgdp))){
    color <- colorPalette(ii, 1)
    palette <- c(palette, color)
    file_name = paste0(parameterSweepPath, '/results/ratio', ratio,'_gtpase', gtpase, '.csv')
    print(file_name)
    if (file.exists(file_name)){

      results <- gtpTurnoverExtract(file_name)
      hydDiff <- results$hydDiff
      ind <- match(NA, hydDiff)
      hydDiff[ind] <- 0
      ftsz <- results$ftsz
      time <- results$time
      if (length(hydDiff)>length(time)){
        hydDiff <- colMeans(hydDiff)
      }
      lines(time, hydDiff, col=color, lwd=scale)
      print(length(time))
      print(length(hydDiff))
    }
  }
  return(palette)
    
}

plotTurnoverSweepVsTime <- function(ratio, bottomoffgdp, gtpase, scale){
  

  points <- length(ratio)
  for (ii in c(1:points)){
    gtpTurnoverVsTimePlot(scale, ratio[ii], gtpase, bottomoffgdp)
  }
}

scale <- 2
cooperativeAssemblyPath <- '/home/lauren/Documents/research/treadmilling-model/figures/gtpTurnoverSweep'
resultsFile <- paste0(cooperativeAssemblyPath, '/gtpTurnoverVsTime.png')
plot.new()

ratio <- c(1, 2, 4, 8)
bottomoffgdp <-  c()
gtpase <- c(0.2, 0.6)

domain <- c(0, 240)
range <- c(0, 30)
if (length(ratio)<2){
  singlePlot <- TRUE
} else{
  singlePlot <- FALSE
}

if (singlePlot){
  scale <- 2
  cols <- 1
  rows <- 1
  plot.new()
  png(resultsFile, width=600, height=600)
  par(mar=c(5, 5, 1, 1))
  plot(0, 0, col='white', xlim=domain, ylim=range, xlab="", ylab="", cex.axis=scale)
} else{
  rows <- length(gtpase)
  cols <- length(ratio)
  png(resultsFile, width=200*cols+200, height=200*rows+100)
  margin_top <- 0.25/rows
  margin_left <- 0.35/cols
}


for (ii in c(1:cols)){
  for (jj in c(1:rows)){

    if (!singlePlot){
      sweepSubplots(margin_top, margin_left, jj, ii, rows, cols, domain, range)
    }
    plotTurnoverSweepVsTime(ratio[ii],  bottomoffgdp, gtpase[jj],scale)
  }
}
palette <- generateColorPalette(length(bottomoffgdp), 1)
label <- addUnits(bottomoffgdp, 2)

if (singlePlot){
  mtext('Time (s)', side=1, line=3, cex=scale)
  grid()
  mtext(expression('GTP Turnover ('*mu*'M GTP/min)'), side=2, line=3, cex=scale)
  
  legend("topleft", legend=label, col=palette, title=expression('k'['b GDP']^'off'), lty=1, cex=scale, lwd=2)
} else{
  ratiolabel <- addUnits(ratio, 3)

  sweep2dLabels(margin_top, margin_left, gtpase, ratio, expression('k'['b GTP']^'off'), 'ratio', "Time (s)", expression('GTP Turnover ('*mu*'M GTP/min)'))  
  margins <- c(0, 0.2, 0, 1)
  par(fig=margins, new=TRUE)
  legend("left", legend=label, col=palette, lty=1, lwd=3, title=expression('k'['b GDP']^'off'), cex=1.5, bty='n')

}
dev.off()

