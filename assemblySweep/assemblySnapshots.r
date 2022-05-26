source("figurePlotter.r")

plotSnapshot <- function(time){
  fileName = paste0('assemblySweep/results/snapshot_', time, '.csv')
  results <- read.csv(fileName, header=FALSE, sep=',')
  results <- as.matrix(results)
  resultsExpand <- matrix(0, nrow=nrow(results), ncol=(ncol(results)*3))
  indExpand <- 2
  for (jj in c(1:ncol(results))){
    resultsExpand[,indExpand] = results[,jj]
    indExpand <- indExpand + 1
    resultsExpand[,indExpand] = results[,jj]
    indExpand <- indExpand + 2
  }
  results <- t(apply(resultsExpand, 2, rev))
  map <- c(colorPalette(9, 0.5), colorPalette(3, 1), colorPalette(5, 1))
  xlabels <- c(0:ncol(results))
  ylabels <- c((0:nrow(results))/3)
  scale <- 1.5
  image(ylabels, xlabels, results, xlim=c(1,30), ylim=c(1, 300), col = map, xlab="", ylab="", cex.axis=scale)
  legend("bottomright", paste0(time, ' s'), bty="n",  cex=2)
}

time <- c('5','15', '20')
margins <- c(3, 3, 1, 1)

for (ii in c(1:3)){
  options(warn=2)
  area <- c(0.33*(ii-1), 0.33*ii, 0, 1)
  x <- try(mtext('bottomleft', ''), silent=TRUE)
  if(class(x)=="try-error"){
    par(fig=area, mar=margins, new=FALSE)
  } else{
    par(fig=area, mar=margins, new=TRUE)
  }

  plotSnapshot(time[ii])
}
