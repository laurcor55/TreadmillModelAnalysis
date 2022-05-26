setwd("../")
source("figurePlotter.r")


pcAssemblyPlot <- function(file_name, color, shape, printError, standAlone){
  if (file.exists(file_name)){
    if (shape==2){
      offset <- 0.1
    } else {
      offset <- 0
    }
    results <- gtpTurnoverExtract(file_name)
    hydDiff <- results$hydDiff
    ftsz <- results$ftsz
    hydDiffMean <- rowMeans(hydDiff[,(ncol(hydDiff)-20):(ncol(hydDiff))])
    if (standAlone){
      scale <- 2
      maxY <- max(hydDiffMean)*1.1
     # maxY <- 30
      plot(ftsz, hydDiffMean+offset, xlim=c(0,max(ftsz)), ylim=c(0, maxY), pch=shape, cex=scale, col=color, ylab='', xlab='', cex.axis=scale)
      grid()
    } else{
      scale <- 2
      points(ftsz, hydDiffMean, col=color, cex=scale, pch=shape)
    }
    ind <- (hydDiffMean>1)
    if (sum(ind)>1){
      model <- lm(hydDiffMean[ind]~ftsz[ind])
      coef <- model$coefficients
      ftsz <- c(0, ftsz)
      fit <- ftsz*coef[2] + coef[1]
      lines(ftsz, fit+offset, col=color, lwd=5)
      print(paste0('GTP Turnover: ', coef[2]))
      print(paste0('Critical Concentration: ', -coef[1]/coef[2]))
    }
    
    if (printError){
      #legend("topright", legend=c(paste0('L = ', s)), box.lty=0)
    }
    if (standAlone){
      mtext(expression('Total FtsZ ('*mu*'M)'), side=1, line=3, cex=scale)
      mtext(expression('PF FtsZ ('*mu*'M)'), side=2, line=3, cex=scale)
    }
   
  }
}


margins <- c(6, 6, 1, 1)
png('pc/pc.png', width=1600, height=450)
plot.new()
par(fig=c(0, 0.5, 0.9, 1), mar=c(0, 6, 0, 1), new=TRUE)

plot(0,0, xlab="",ylab="", col="white", axes=FALSE)
box()
label <- c("no PC190723", "Binds to PF", "Binds to PF and monomer")
wid <- c(0, 0.45, 0.45)
legend("left", legend=label, hor=TRUE, pch=c(1, 2, 0),bty="n", col=generateColorPalette(3,1), cex=1.5, text.width=wid, lwd=2)

par(fig=c(0, 0.25, 0, 0.9), mar=margins, new=TRUE)
file_name <- 'pc/results/assemblyPc_knuc15000_sameKinetics.csv'
pcAssemblyPlot(file_name, colorPalette(2, 1), 2, FALSE, TRUE)
file_name <- 'pc/results/assemblyPc_noPc.csv'
pcAssemblyPlot(file_name, colorPalette(1, 1), 1, FALSE, FALSE)
file_name <- 'pc/results/assemblyPc_stableT.csv'
pcAssemblyPlot(file_name, colorPalette(3, 1), 0, FALSE, FALSE)
legend("topleft", 'a)', bty="n", inset=c(-0.1, 0), cex=2)

par(fig=c(0.25, 0.5, 0, 0.9), mar=margins, new=TRUE)
file_name <- 'pc/results/gtpTurnoverPc_noPc.csv'
gtpTurnoverPlot(file_name, colorPalette(1, 1), 1, FALSE, TRUE)
file_name <- 'pc/results/gtpTurnoverPc_knuc15000_sameKinetics.csv'
gtpTurnoverPlot(file_name, colorPalette(2, 1), 2, FALSE, FALSE)
file_name <- 'pc/results/gtpTurnoverPc_stableT.csv'
gtpTurnoverPlot(file_name, colorPalette(3, 1), 0, FALSE, FALSE)
legend("topleft", 'b)', bty="n", inset=c(-0.1, 0), cex=2)

par(fig=c(0.5, 0.75, 0, 0.9), mar=margins, new=TRUE)
velocityPlot('pc/results/velocityPc.csv', c(-5,5), c(0,30), TRUE)
legend("topleft", 'c)', bty="n", inset=c(-0.1, 0), cex=2)


par(fig=c(0.75,1, 0, 0.9), mar=margins, new=TRUE)
pfLengthPlot('pc/results/pfLengthPc_stableT.csv', 3)
legend("topleft", 'd)', bty="n", inset=c(-0.1, 0), cex=2)


dev.off()
