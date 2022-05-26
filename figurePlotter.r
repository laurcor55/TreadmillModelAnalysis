library(beeswarm)
library(dplyr)
library(png)

diagramPlot <- function(file_name){
  im <- readPNG(file_name)
  plot.new()
  rasterImage(im, 0, 0, 1, 0.8)
}


gtpTurnoverVsVelocityPlot <- function(file_name){
  df <- read.csv(file_name, header=TRUE, sep=',')
  scale <- 1.25

  #beeswarm(vel~hyd, data=df, pch=16, cex=0.75, xlab="", ylab="", cex.axis=scale*1.25, method='swarm')
  #bxplot(vel~hyd, data=df, outline=FALSE,  add=TRUE)
  #plot(vel~hyd, data=df, xlab="", ylab="",  cex.axis=scale*1.25)
  velocitymean <- aggregate(df[, 2], list(df$hyd), function(x) mean = mean(x))

  plot(velocitymean[,1], velocitymean[,2], xlab="", ylab="", ylim=c(-1, 8))
  velocitysd <- aggregate(df[, 2], list(df$hyd), function(x) sd =sd(x))
  arrows(x0=velocitysd[,1], y0=velocitymean[,2] - velocitysd[,2], x1=velocitysd[,1], y1=velocitymean[,2] + velocitysd[,2], code=3, angle=90, length=0.05, lwd=scale)
  grid()
  abline(h=0, col=colorPalette(9, 0.8), lty=5, lwd=1)
  abline(v=0, col=colorPalette(9, 0.8), lty=5, lwd=1)

  mtext(expression('k'['hyd']~'('*s^-1*')'), side=1, line=2)
  mtext('PF Velocity (subunits/s)', side=2, line=2)
}

pfMixingPlot <- function(file_name, experimentalFileName, tau, magnitude, standAlone){
  if (file.exists(file_name)){
    df <- read.csv(file_name, header=FALSE, sep=',', dec='.', stringsAsFactors=FALSE)
    df <- data.matrix(df)
    maxX <- max(df)
    maxY <- max(df[2,])*1.1
    if (standAlone){
      plot(0, 0, col="white", xlim=c(0,maxX), ylim=c(0,maxY), ylab='', xlab='', main='')
    }
    time <- df[1,]
    if (ncol(df)>2){
      mixing <- df[2:nrow(df),]
      mixing <-  t(apply(mixing, 2, rev))
      boundary <- assemblyBoundary(mixing)
      color <- colorPalette(1, 0.5)
      polygon(c(time,rev(time)),c(boundary$upper,rev(boundary$lower)),col=color, border=color)
    } else{
      mixing <- df[2,]
      lines(time, mixing, col=color)
    }
    
    dfExp <- read.csv(experimentalFileName, header=FALSE, sep=',', dec='.', stringsAsFactors=FALSE)
    dfExp <- data.matrix(dfExp)
    time <- dfExp[,1]
    fluor <- dfExp[,2]
    fluor <- -1*fluor + 1
    fit <- exp(-time/tau)
    fit <- -magnitude*fit+magnitude
    lines(time, fluor, col='black', lwd=2, lty=2)
    if (standAlone){
      mtext('Time (s)', side=1, line=2)
      mtext('Proportion of Mixed PF', side=2, line=2)
      grid()
    }
    
  }
  
}



pfMixingSweepPlot <- function(parameterSweepPath, knuc, khyd, fragment, anneal){
  resultsFile <- paste0(parameterSweepPath, '/pfMixingSweep.png')
  cols <- length(anneal)
  rows <- length(fragment)
  png(resultsFile, width=250*cols+100, height=200*rows+100)
  
  margin_top <- 0.25/rows
  margin_left <- 0.2/cols
  domain <- c(0, 550)
  range <- c(0, 1)
  experimentalFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005a/chen2005_fig4b.csv'
  
  plot.new()
  for (ii in c(1:rows)){
    for (jj in c(1:cols)){
      sweepSubplots(margin_top, margin_left, ii, jj, rows, cols, domain, range)
      fileName <- paste0(parameterSweepPath, '/results/knuc', knuc,'_khyd', khyd, '_fragment', fragment[ii], '_anneal', anneal[jj], '.csv')
      pfMixingPlot(fileName, experimentalFileName, 127, 0.4, FALSE)
    }
  }
  #expression('k'['hyd'])
  sweep2dLabels(margin_top, margin_left, fragment, anneal, 'fragment', 'anneal', expression('Time (s)'), expression('Monomeric FtsZ ('*mu*'M)'))  
  dev.off()
}

pfMixingPlotMulti <- function(fileNames){
  maxX <- 550
  maxY <- 1.2
  scale <- 2
  plot(0, 0, col="white", xlim=c(0,maxX), ylim=c(0,maxY), ylab='', xlab='', cex.axis=scale*1.25, main='')
  
  for (ii in c(1:2)){
    df <- read.csv(fileNames[ii], header=FALSE, sep=',', dec='.', stringsAsFactors=FALSE)
    df <- data.matrix(df)
    time <- df[1,]
    if (ncol(df)>2){
      mixing <- df[2:nrow(df),]
      mixing <-  t(apply(mixing, 2, rev))
      boundary <- assemblyBoundary(mixing)
      color <- colorPalette(ii, 0.5)
      polygon(c(time,rev(time)),c(boundary$upper,rev(boundary$lower)),col=color, border=color)
    } else{
      mixing <- df[2,]
      lines(time, mixing, col=color)
    }
  }
  label <- addUnits(c(20, 240), 4)
  legend("topright", title='Pre-Assembly Time', legend=label, col=generateColorPalette(2, 1), lty=1, lwd=3, cex=1.5)
   
  mtext('Time (s)', side=1, line=3, cex=scale*1.25)
  mtext('Proportion of Mixed PF', side=2, line=3, cex=scale*1.25)
  grid()
}


disassemblyHalftime <- function(file_name, startingFtsz){
  results <- assemblyExtract(file_name)
  monomerConc <- results$monomerConc
  time <- results$time
  startingFtszModel <- monomerConc[1,]
  ftszDiff <- abs(startingFtszModel-startingFtsz)
  runChoice <- which(ftszDiff < 0.05)
  singleConc <- monomerConc[2:nrow(monomerConc),runChoice]
  signalDiff <- abs(singleConc[1,]-singleConc[nrow(singleConc),])
  runChoice <- which(signalDiff > 0.05)
  
  singleConc <- singleConc[,runChoice]
  singleConc <- rowMeans(singleConc)
  finalConc <- mean(singleConc[length(singleConc)])
  startingConc <- mean(singleConc[1])
  middleConc <- (finalConc - startingConc)/2
  ind <- which.min(abs(middleConc - singleConc))
  print(paste('Half-time:', time[ind]))
}

disassemblySweepPlot <- function(parameterSweepPath, knuc, khyd, fragment, anneal, preassembleTime){
  resultsFile <- paste0(parameterSweepPath, '/disassemblySweep.png')
  cols <- length(anneal)
  rows <- length(fragment)
  png(resultsFile, width=250*cols+100, height=200*rows+100)
  
  margin_top <- 0.25/rows
  margin_left <- 0.2/cols
  domain <- c(0, 120)
  range <- c(0, 8)
  
  plot.new()
  for (ii in c(1:rows)){
    for (jj in c(1:cols)){
      sweepSubplots(margin_top, margin_left, ii, jj, rows, cols, domain, range)
      fileName1 <- paste0(parameterSweepPath, '/results/knuc', knuc,'_khyd', khyd, '_fragment', fragment[ii], '_anneal', anneal[jj], '_preassemble', preassembleTime[1], '.csv')
      fileName2 <- paste0(parameterSweepPath, '/results/knuc', knuc,'_khyd', khyd, '_fragment', fragment[ii], '_anneal', anneal[jj], '_preassemble', preassembleTime[2], '.csv')
      fileNames <- c(fileName1, fileName2)
      disassemblyPlot(fileNames, preassembleTime, c(12, 90), FALSE)
    }
  }
  #expression('k'['hyd'])
  sweep2dLabels(margin_top, margin_left, fragment, anneal, 'fragment', 'anneal', expression('Time (s)'), expression('Monomeric FtsZ ('*mu*'M)'))  
  dev.off()
}


disassemblyPlot <- function(fileNames, times, tau, standAlone){
  xmax <- max(read.csv(fileNames[1]))
  if (standAlone){
    plot(0, 0, xlim=c(0, xmax), ylim=c(0, 7), col="white", ylab='', xlab='')#, cex.axis=scale*1.25)
    grid()
    abline(h=0, col=colorPalette(9, 0.8), lty=5, lwd=1.5)
    abline(v=0, col=colorPalette(9, 0.8), lty=5, lwd=1.5)
  }
  for (ii in c(1:length(fileNames))){
    fileName <- fileNames[ii]
    assemblyPlot(fileName, c(numeric(ii-1), 6), FALSE)
    disassemblyHalftime(fileName, 6)
  }
  alpha <- 6

  label <- c()
  for (ii in c(1:length(tau))){
    if (tau[ii]>0){
      x <- seq(0, xmax, length.out=100)
      fit <- -(alpha*exp(-x/tau[ii]))+alpha
      lines(x, fit, col='black', lty=2, lwd=3)
      label <- c(label, addUnits(times[ii], 4))
    }
  }
  if (length(times)>1){
    legend("bottomright", title='Pre-Assembly', legend=label, col=generateColorPalette(2, 1), lty=1, lwd=5, inset=c(0, 0.1), cex=0.75)
  }
  if (standAlone){
    mtext('Time (s)', side=1, line=2)#, cex=scale*1.25)
    mtext(expression('Monomeric FtsZ ('*mu*'M)'), side=2, line=2)#, cex=scale*1.25)
  }
}


##################################################################
################## SS PF CONC ##########################


ssPfConcSweepPlot <- function(parameterSweepPath, experimentalFileName, startingFtsz, kswitch, gdpexchange, gtpase, topoff, caponpf, bottomon, bottomoff, mciz){
  resultsFile <- paste0(parameterSweepPath, '/parameterSweepResult.png')
  rows <- length(gtpase)
  cols <- length(kswitch)
  png(resultsFile, width=250*cols+100, height=200*rows+100)

  maxConc <- max(startingFtsz)+1
  margin_top <- 0.25/rows
  margin_left <- 0.25/cols
  domain <- c(0, maxConc)
  range <- c(0, maxConc)
  plot.new()
  for (ii in c(1:cols)){
    for (jj in c(1:rows)){
      file_name = paste0(parameterSweepPath, '/results/kswitch', kswitch[ii],'_gtpase', gtpase[jj], '_topoff', topoff,'_caponpf', caponpf, '_gdpexchange', gdpexchange, '_bottomon', bottomon, '_bottomoff', bottomoff,  '_mciz', mciz, '_array', '.csv') #
      sweepSubplots(margin_top, margin_left, jj, ii, rows, cols, domain, range)
      ssPfConcPlot(file_name, 'black', 1, FALSE)
    }
  }
  sweep2dLabels(margin_top, margin_left, gtpase, kswitch, "khyd", "Knuc", expression('Total FtsZ ('*mu*'M)'), expression('PF FtsZ ('*mu*'M)'))  
  dev.off()
}

ssPfConcPlot <- function(file_name, color, shape, standAlone){
  results <- assemblyExtract(file_name)
  monomerConc <- results$monomerConc
  time <- results$time
  finalMonomerConc <- colMeans(monomerConc[round(length(time)*2/3):length(time),])
  startingFtszModel <- monomerConc[1,]
  pfFtsz <- startingFtszModel - finalMonomerConc
  if (standAlone){
    plot(startingFtszModel, pfFtsz, ylab='', xlab='', pch=shape, xlim=c(0,max(startingFtszModel)), ylim=c(0,max(startingFtszModel)), col=color)
    grid()
    mtext(expression('Total FtsZ ('*mu*'M)'), side=1, line=2)
    mtext(expression('PF FtsZ ('*mu*'M)'), side=2, line=2)
  } else{
    points(startingFtszModel, pfFtsz, col=color,  pch=shape)
  }
  criticalConcentration(startingFtszModel, pfFtsz, color)
  abline(h=0, col=colorPalette(9, 0.8), lty=5, lwd=1)
  abline(v=0, col=colorPalette(9, 0.8), lty=5, lwd=1)
  
}

ssPfConcLoss <- function(file_name, expectedCc){
  results <- assemblyExtract(file_name)
  monomerConc <- results$monomerConc
  time <- results$time
  finalMonomerConc <- colMeans(monomerConc[(length(time)-30):length(time),])
  startingFtszModel <- monomerConc[1,]
  pfFtsz <- startingFtszModel - finalMonomerConc
  ind <- (pfFtsz>0.2)
  if (sum(ind)>1){
    model <- lm(pfFtsz[ind]~startingFtszModel[ind])
    coef <- model$coefficients
    actualCc <- -coef[1]/coef[2]
    
  }
  loss <- abs(expectedCc - actualCc)
  return(loss)
}

##################################################################
################## ASSEMBLY ##########################

assemblyExtract <- function(fileName){
  results <- matrix(0, nrow=4, ncol=4)
  if (file.exists(fileName)){
    results <- read.csv(fileName, header=FALSE, sep=',')
    results <- as.matrix(results)
  } 
  time <- results[2:nrow(results),1]
  monomerConc <- results[, 2:ncol(results)]
  ind <- which(monomerConc==NA)
  #monomerConc[ind] <- 0
  results <- list('time'=time, 'monomerConc'=monomerConc)
  return(results)
}

assemblyPlot <- function(file_name, startingFtsz, standAlone){
  results <- assemblyExtract(file_name)
  monomerConc <- results$monomerConc
  time <- results$time
  startingFtszModel <- monomerConc[1,]
  if (standAlone){
    plot(0, 0, xlim=c(0,max(time)), ylim=c(0, max(startingFtsz)*1.3), col="white", ylab='', xlab='')#, cex.axis=scaleAll*0.75)
    grid()
    abline(h=0, col=colorPalette(9, 0.8), lty=5, lwd=1)
    abline(v=0, col=colorPalette(9, 0.8), lty=5, lwd=1)
  }
  for (kk in c(1:length(startingFtsz))){
    ftszDiff <- abs(startingFtszModel-startingFtsz[kk])
    runChoice <- which(ftszDiff < 0.05)
    singleConc <- monomerConc[2:nrow(monomerConc),runChoice]
    signalDiff <- abs(singleConc[1,]-singleConc[nrow(singleConc),])
    runChoice <- which(signalDiff > 0.05)
    color <- colorPalette(kk, 0.6)
    if (length(runChoice)==1){
      lines(time, singleConc[, runChoice], col=color, lwd=3)
    }
    else {
      singleConc <- singleConc[,runChoice]
      boundary <- assemblyBoundary(singleConc)
      polygon(c(time,rev(time)),c(boundary$upper,rev(boundary$lower)), lwd=3,col=color, border=color)
    }
  }
  if (standAlone){
    label <- addUnits(startingFtsz, 1)
    legend("topright", title='Total FtsZ', legend=label, col=generateColorPalette(3, 0.8), lty=1, lwd=3, cex=0.75)
    mtext('Time (s)', side=1, line=2)#, cex=scaleAll*1.25)
    mtext(expression('Monomeric FtsZ ('*mu*'M)'), side=2, line=2)#, cex=scaleAll*1.25)
  }
}

assemblyMultiPlot <- function(parameterSweepPath, modelFileName){
  results <- assemblyExtract(modelFileName)
  monomerConc <- results$monomerConc
  time <- results$time

  plot(0, 0, col="white", xlim=c(0,15), ylim=c(0, 4), xlab='', ylab='')
  grid()
  startingFtszModel <- monomerConc[1,]
  startingFtsz <- unique(startingFtszModel)
  for (kk in c(1:length(startingFtsz))){
    ftszDiff <- abs(startingFtszModel-startingFtsz[kk])
    runChoice <- which(ftszDiff < 0.05)
    singleConc <- monomerConc[2:nrow(monomerConc),runChoice]
    signalDiff <- (singleConc[1,]-singleConc[nrow(singleConc),])
    runChoice <- which(signalDiff > 0.05)
    color <- colorPalette(kk, 1)
    for (ll in c(1:length(runChoice))){
      lines(time, singleConc[, runChoice[ll]], col=color, lwd=2)
    }
  }
  ftszUnique <- unique(startingFtsz)
  print(ftszUnique)
  label <- addUnits(ftszUnique, 1)
  legend("topright", title='Total FtsZ', legend=label, col=generateColorPalette(3, 1), lty=1, lwd=3, cex=0.75)
    
  mtext(expression('Monomeric FtsZ ('*mu*'M)'), side=2, line=2)
  mtext('Time (s)', side=1, line=2)
}

assemblyBoundary <- function(singleConc){
  sampleCount <- ncol(singleConc)
  dataPoints <- nrow(singleConc)*ncol(singleConc)
  rowMean <- rowMeans(singleConc)
  rowStd <- rowSds(singleConc)
  upper <- rowMean + 1.96*rowStd/sqrt(sampleCount)
  lower <- rowMean - 1.96*rowStd/sqrt(sampleCount)
  boundary <- list("upper"=upper, "lower"=lower)
  return(boundary)
}

assemblyExtractExperimental <- function(fileName){
  results <- matrix(4)
  if (file.exists(fileName)){
    results <- read.csv(fileName, header=FALSE, sep=',', dec='.', stringsAsFactors=FALSE)
  }
  time <- results[,1]
  monomerConc <- results[, 2:ncol(results)]
  results <- list('time'=time, 'monomerConc'=monomerConc)
  return(results)
}

assemblyPlotExperimental <- function(file_name, startingFtsz){
  results <- read.csv(file_name, header=FALSE, sep=',', dec='.', stringsAsFactors=FALSE)
  results <- as.matrix(results)
  time <- results[,1]
  monomerConc <- results[, 2:ncol(results)]
  for (ii in c(1:length(startingFtsz))){
    ftszDiff <- abs(monomerConc[1,]-startingFtsz[ii])
    runChoice <- which.min(ftszDiff)
    color <- colorPalette(ii, 1.7)
    lines(time, monomerConc[,runChoice], col=color, lty=2, lwd=2)
  }
}

assemblyLoss <- function(filenameExp, filenameModel, startingFtsz, totalTime){
  loss <- 100
  if (file.exists(filenameModel)){
    resultsExp <- assemblyExtractExperimental(filenameExp)
    resultsModel <- assemblyExtract(filenameModel)
    startingFtszModel <- resultsModel$monomerConc[1,]
    startingFtszExp <- resultsExp$monomerConc[1,]
    timeSamples <- seq(0.1, totalTime, length=100)
    if (sum(resultsModel$monomerConc[1:2, 2])>0.1){
      loss <- 0
      for (kk in c(1:length(startingFtsz))){
        ftszDiff <- abs(startingFtszModel-startingFtsz[kk])
        runChoice <- which(ftszDiff < 0.05)
        if (length(runChoice)>0){
          singleConc <- resultsModel$monomerConc[,runChoice]
          if (length(runChoice)>1){
            singleConc <- rowMeans(singleConc)
          }
          runChoice <- which.min(abs(startingFtszExp-startingFtsz[kk]))
          concExp <- resultsExp$monomerConc[,runChoice]
          for (ii in c(1:length(timeSamples))){
            expTimePoint <- which.min(abs(resultsExp$time - timeSamples[ii]))
            modelTimePoint <- which.min(abs(resultsModel$time - timeSamples[ii]))
            if (length(modelTimePoint)>1){
              mean <- mean(singleConc[modelTimePoint,])
            }
            else{
              mean <- singleConc[modelTimePoint]
            }
            distance <- abs(concExp[expTimePoint] - mean)
            loss <- loss + distance
          }
        }
      }
      loss <- loss/(100*length(startingFtsz))
    }
  }
  return(loss)
}


assemblySweepPlot <- function(parameterSweepPath, experimentalFileName, startingFtsz, kswitch, gdpexchange, gtpase, topoff, caponpf, bottomon, bottomoff, mciz, maxTime){
  resultsFile <- paste0(parameterSweepPath, '/parameterSweepResult.png')
  rows <- length(gtpase)
  cols <- length(kswitch)
  png(resultsFile, width=250*cols+100, height=200*rows+100)

  resultsExp <- assemblyExtractExperimental(experimentalFileName)
  maxConc <- max(startingFtsz)+1
  margin_top <- 0.25/rows
  margin_left <- 0.25/cols
  domain <- c(0, maxTime)
  range <- c(0, max(startingFtsz))
  plot.new()
  for (ii in c(1:cols)){
    for (jj in c(1:rows)){
      file_name = paste0(parameterSweepPath, '/results/kswitch', kswitch[ii],'_gtpase', gtpase[jj], '_topoff', topoff,'_caponpf', caponpf, '_gdpexchange', gdpexchange, '_bottomon', bottomon, '_bottomoff', bottomoff,  '_mciz', mciz, '.csv') #
      sweepSubplots(margin_top, margin_left, jj, ii, rows, cols, domain, range)
      assemblyPlot(file_name, startingFtsz, FALSE)
      assemblyPlotExperimental(experimentalFileName, startingFtsz)
      loss <- assemblyLoss(experimentalFileName, file_name, startingFtsz, maxTime)
      print(loss)
      text(maxTime*0.4, maxConc-0.5, paste0('L = ',round(loss, 3)))
    }
  }
  sweep2dLabels(margin_top, margin_left, gtpase, kswitch, expression('k'['hyd']), expression('K'['nuc']), "Time (s)", expression('Monomeric FtsZ ('*mu*'M)'))  
  dev.off()
}




#####################################################################
########## GTP TURNOVER #########

gtpTurnoverSweepPlot <- function(parameterSweepPath, knuc, khyd, fragment, anneal, mciz){
  resultsFile <- paste0(parameterSweepPath, '/gtpTurnoverSweep.png')
  cols <- length(anneal)
  rows <- length(knuc)
  png(resultsFile, width=250*cols+100, height=200*rows+100)
  
  margin_top <- 0.25/rows
  margin_left <- 0.2/cols
  domain <- c(0, 9)
  range <- c(0, 40)
  
  plot.new()
  for (ii in c(1:rows)){
    for (jj in c(1:cols)){
      sweepSubplots(margin_top, margin_left, ii, jj, rows, cols, domain, range)
      for (kk in c(1:length(mciz))){
        color <- colorPalette(kk, 1)
        file_name = paste0(parameterSweepPath, '/results/knuc', knuc[ii],'_khyd', khyd, '_fragment', fragment, '_anneal', anneal[jj], '_mciz', mciz[kk], '.csv')
       # gtpTurnoverVsTimePlot(file_name)
        gtpTurnoverPlot(file_name, color, 1, FALSE, FALSE)
      }
      
    }
  }
  #expression('k'['hyd'])
  sweep2dLabels(margin_top, margin_left, knuc, anneal, 'knuc', expression('k'['anneal']), expression('Total FtsZ ('*mu*'M)'), expression('GTP Turnover ('*mu*'M GTP/min)'))  
  dev.off()
}

gtpTurnoverSweepBottomCapPlot <- function(parameterSweepPath, topoff, gtpase, kswitch, bottomon, mciz, ktongtp, bottomoffcap, bottomoncap, bottomoffgdp){
  resultsFile <- paste0(parameterSweepPath, '/gtpTurnoverSweep.png')
  cols <- length(bottomoncap)
  rows <- length(bottomoffcap)
  png(resultsFile, width=250*cols+100, height=200*rows+100)
  
  margin_top <- 0.25/rows
  margin_left <- 0.2/cols
  domain <- c(0, 9)
  range <- c(0, 40)
  
  plot.new()
  for (ii in c(1:rows)){
    for (jj in c(1:cols)){
      sweepSubplots(margin_top, margin_left, ii, jj, rows, cols, domain, range)
      for (kk in c(1:length(ktongtp))){
        color <- colorPalette(kk, 1)
        file_name = paste0(parameterSweepPath, '/results/kswitch', kswitch,'_gtpase', gtpase, '_bottomoncap', bottomoncap[jj], '_ktongtp', ktongtp[kk],'_bottomoffcap', bottomoffcap[ii], '_bottomon', bottomon, '_topoff', topoff, '_bottomoffgdp', bottomoffgdp, '_mciz', mciz, '.csv')
        
        gtpTurnoverPlot(file_name, color, 1, FALSE, FALSE)
      }
      
    }
  }
  
  sweep2dLabels(margin_top, margin_left, bottomoffcap, bottomoncap, expression('k'['off']~'Cap'), expression('k'['on']~'Cap'), expression('Total FtsZ ('*mu*'M)'), expression('GTP Turnover ('*mu*'M GTP/min)'))  
  dev.off()
}


gtpTurnoverLoss <- function(file_name, idealCc, idealGtpase, maxFtsz){
  idealFtsz <- seq(0, maxFtsz, length=100)
  ind <- which.min(abs(idealFtsz-idealCc))
  idealFit <- idealFtsz*idealGtpase-idealCc*idealGtpase
  idealFit[1:ind] <- 0
  if (file.exists(file_name)){
    results <- gtpTurnoverExtract(file_name)
    hydDiff <- results$hydDiff
    ftsz <- results$ftsz
    hydDiffMean <- rowMeans(hydDiff)
    loss <- 0 
    for (mm in c(1:length(ftsz))){
      ind <- which.min(abs(idealFtsz-ftsz[mm]))
      loss <- loss + abs(hydDiffMean[mm] - idealFit[ind])
    }
    loss <- round(loss/length(ftsz), 2)
  }
  else{
    loss <- 100
  }
  return(loss)
}




gtpTurnoverExtract <- function(file_name){
  s <- 0
  results <- matrix(0, nrow=4, ncol=24)
  if ((file.exists(file_name)) & (file.size(file_name)>0)){
    resultsOriginal <- read.csv(file_name, header=FALSE, sep=',', stringsAsFactors=FALSE)
    if (nrow(resultsOriginal)==1){
      results <- matrix(0, nrow=2, ncol=ncol(resultsOriginal))
      results[1,] <- as.numeric(resultsOriginal)
      results[2,] <- as.numeric(resultsOriginal)
    } else{
      results <- as.matrix(resultsOriginal)
      ind <- which(results=='NaN')
      results[ind] <- 0
    }
  }
  time <- results[1,2:ncol(results)]
  ftsz <- results[2:nrow(results),1]
  hydDiff <- results[2:nrow(results),2:ncol(results)]
  results <- list("hydDiff"=hydDiff, "ftsz"=ftsz, "time"=time)
  return(results)
}

gtpTurnoverPlot <- function(file_name, color, shape, printError, standAlone){
  if (file.exists(file_name)){
    results <- gtpTurnoverExtract(file_name)
    hydDiff <- results$hydDiff
    ftsz <- results$ftsz
    hydDiffMean <- rowMeans(hydDiff[,round(ncol(hydDiff)*0.5):(ncol(hydDiff))])
    if (standAlone){
      maxY <- 40
      plot(ftsz, hydDiffMean, xlim=c(0,max(ftsz)), ylim=c(0, maxY), pch=shape, col=color, ylab='', xlab='')
      grid()
      

    } else{
      points(ftsz, hydDiffMean, col=color, pch=shape)
    }
    criticalConcentration(ftsz, hydDiffMean, color)
    abline(h=0, col=colorPalette(9, 0.8), lty=5, lwd=1)
    abline(v=0, col=colorPalette(9, 0.8), lty=5, lwd=1)
    if (standAlone){
      mtext(expression('Total FtsZ ('*mu*'M)'), side=1, line=2)
      mtext(expression('GTP Turnover ('*mu*'M GTP/min)'), side=2, line=2)
    }
  }
}

gtpTurnoverVsHydPlot <- function(){
  fileStart <- 'gtpTurnoverSweep/results/kswitch30003_gtpase'
  fileEnd <- '_bottomoncap0_ktongtp1_bottomoffcap0_bottomon10_topoff6.5_bottomoffgdp2.7_bottomoffgtp5_mciz0.csv'
  gtpase <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  scale <- 1.25
  plot(0, 0, xlim=c(0, 0.6), ylim=c(0, 15), col="white", main="", xlab="", ylab="", cex.axis=scale*1.25)
  for (ii in c(1:length(gtpase))){
    fileName <- paste0(fileStart, gtpase[ii], fileEnd)
    results <- gtpTurnoverExtract(fileName)
    hydDiff <- results$hydDiff
    ftsz <- results$ftsz
    means <- rowMeans(hydDiff[,round(ncol(hydDiff)*0.5):(ncol(hydDiff))])
    sdev <- sd(means)
    means <- mean(means)
    khyd <- rep(gtpase[ii], length(means))
    points(khyd, means, cex=scale)
    if (sdev >0.1){
      arrows(khyd, means-sdev, khyd, means+sdev, length=0.05, angle=90, code=3)

    }
    
  } 
  abline(h=0, col=colorPalette(9, 0.8), lty=5, lwd=1)
  abline(v=0, col=colorPalette(9, 0.8), lty=5, lwd=1)
  mtext(expression('k'['hyd']~'('*s^-1*')'), side=1, line=4, cex=scale*1.25)
  mtext(expression('GTP Turnover ('*mu*'M GTP/min)'), side=2, line=3, cex=scale*1.25)
  grid()
}

gtpTurnoverVsTimePlot <- function(file_name){
  if (file.exists(file_name)){
    results <- gtpTurnoverExtract(file_name)
    hydDiff <- results$hydDiff
    ind <- match(NA, hydDiff)
    hydDiff[ind] <- 0
    ftsz <- results$ftsz
    ftszUnique <- unique(ftsz)
    time <- results$time
    scale <- 1.25

    plot(0, 0, col="white", xlim=c(0, 300), ylim=c(0, max(hydDiff)*1.3), lwd=scale, ylab='', xlab='', cex.axis=scale*1.25, main='')

    for (ii in c(1:length(ftszUnique))){
      ind <- which(ftsz == ftszUnique[ii])
      if (length(ind)>1){
        hydDiffUnique <- colMeans(hydDiff[ind,])
      } else{
        hydDiffUnique <- hydDiff

      }
      lines(time, hydDiffUnique, col=colorPalette(ii, 1), lwd=3)
    }
    ftsz <- addUnits(ftszUnique, 1)
    legend('topright', legend=ftsz, title="FtsZ", col=generateColorPalette(2), lty=1, lwd=3, cex=1.25)
    
    grid()
    mtext(expression('GTP Turnover ('*mu*'M GTP/min)'), side=2, line=3, cex=scale*1.25)
    mtext('Time (s)', side=1, line=3, cex=scale*1.25)

  }
}


#################################################
########### PF LENGTH ###########################

pfLengthPlot <- function(file_name, colorNumber){
  df <- read.csv(file_name, header=TRUE, sep=',', dec='.', stringsAsFactors=TRUE)
  #boxplot(df$len[df$time=="Nucleation"], df$len[df$time=="Steady State"], names=c("Nucleation", "Steady State"), ylab='PF Length', cex.lab=scale, cex.axis=scale, cex.main=scale, cex.sub=2)
  df <- filter(df, time>60)
  maxX <- max(df$len)+20
  color <- colorPalette(colorNumber, 0.4)
  hist(df$len, breaks=seq(0, maxX, 10), ylim=c(0, 90), col=color, ylab='', xlab='', main='')
  #boxplot(df$len[df$ftsz=="2"], df$len[df$ftsz=="3"], names=c("2 uM", "3 uM"), ylab='PF Length (subunits)', cex.lab=scale, cex.axis=scale, cex.main=scale, cex.sub=2)
  grid()
  mtext('PF Length (subunits)', side=1, line=2)
  mtext('Frequency', side=2, line=2)

  print(paste('Length: ', mean(df$len)))
  print(paste('Standard Deviation: ',sd(df$len)))
 # vioplot(df$len[df$time=="Nucleation"], df$len[df$time=="Steady State"], names=c("Nucleation", "Steady State"))
}


pfLengthPlotMciz <- function(file_name){
  print(file_name)
  if (file.exists(file_name)){
    df <- read.csv(file_name, header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE)
    df <- filter(df, time>60)
    if (nrow(df)>1){
      values <- length(unique(df$mciz))
      color <- generateColorPalette(values, 0.5)
      scale <- 1.25
      bxplot(len~mciz, data=df, outline=FALSE, xlab="", ylab="", type = "n")
      beeswarm(len~mciz, data=df, pch=16, col=color, add=TRUE, cex=0.5)
      grid()
      mtext(expression('MciZ ('*mu*'M)'), side=1, line=2)
      mtext(expression('PF Length (subunits)'), side=2, line=2)
      print(aggregate(df[, 3], list(df$mciz), function(x) c(mean = mean(x), sd = sd(x))))
    }
  }
}


pfLengthPlotTime <- function(file_name){
  scale <- 1.25
  df <- read.csv(file_name, header=TRUE, sep=',', dec='.', stringsAsFactors=TRUE)
  color <- generateColorPalette(2, 0.5)

  bxplot(len~time, data=df, outline=FALSE, xlab="", ylab="", type = "n", cex.axis=scale*1.25)
  beeswarm(len~time, data=df, pch=16, col=color, add=TRUE, cex=0.5)
  print(t.test(len~time, data=df))
  grid()
  mtext('Pre-Assembly Time (s)', side=1, line=3, cex=scale*1.25)
  mtext('PF Length (subunits)', side=2, line=3, cex=scale*1.25)

  print(paste('Length 20 s: ', mean(df$len[df$time==20])))
  print(paste('SD: ', sd(df$len[df$time==20])))

  print(paste('Length 4 min: ', mean(df$len[df$time>100])))
  print(paste('SD: ', sd(df$len[df$time>100])))
 # vioplot(df$len[df$time=="Nucleation"], df$len[df$time=="Steady State"], names=c("Nucleation", "Steady State"))
}

###########################################################
######### VELOCITY ###########################

velocitySweepPlot <- function(parameterSweepPath, gtpase, topoff, bottomoncap, kswitch, bottomon, topongtp,bottomoffcap){
  resultsFile <- paste0(parameterSweepPath, '/velocitySweep.png')
  cols <- length(knuc)
  rows <- length(anneal)
  png(resultsFile, width=250*cols+100, height=200*rows+100)
  
  margin_top <- 0.25/rows
  margin_left <- 0.2/cols
  domain <- c(0, 9)
  range <- c(0, 40)
  
  plot.new()
  for (ii in c(1:rows)){
    for (jj in c(1:cols)){
      sweepSubplots(margin_top, margin_left, ii, jj, rows, cols, domain, range)
      file_name = paste0(parameterSweepPath, '/results/knuc', knuc[jj],'_khyd', khyd, '_fragment', fragment, '_anneal', anneal[ii],  '.csv')
      print(file_name)
      velocityPlotMciz(file_name, FALSE)

    }
  }
  #expression('k'['hyd'])
  sweep2dLabels(margin_top, margin_left, anneal, knuc, 'anneal', expression('K'['nuc']), expression('Total FtsZ ('*mu*'M)'), expression('GTP Turnover ('*mu*'M GTP/min)'))  
  dev.off()
}

velocityPlotMciz <- function(file_name, standAlone){
  if (file.exists(file_name)){
    df <- read.csv(file_name, header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE)
    if (nrow(df)>1){
      df <- filter(df, cap<1)
      values <- length(unique(df$cap))
      color <- generateColorPalette(values, 0.8)
      scale <- 1.25
      if (standAlone){
        beeswarm(bottom~cap, ylim=c(0, 20), data=df, corral='wrap', xlab="", ylab="", pch=16, cex=0.3, col=color, add=FALSE)

      } else {
        beeswarm(bottom~cap, data=df, xlab="", ylab="", cex.axis=scale*1.25, pch=16, cex=0.3, col=color, add=TRUE)

      }

      #bxplot(center~cap, data=df, outline=FALSE, xlab="", ylab="", type = "n", cex.axis=scale*1.25)
      #beeswarm(center~cap, data=df, pch=16, col=color, add=TRUE, cex=0.75)
      grid()
      mtext(expression('MciZ ('*mu*'M)'), side=1, line=2)
      mtext(expression('Velocity (subunits/s)'), side=2, line=2)
      print(aggregate(df[, 2], list(df$cap), function(x) c(mean = mean(x), sd = sd(x))))
      df <- filter(df, bottom>4)
      print(aggregate(df[, 2], list(df$cap), function(x) c(mean = mean(x), sd = sd(x))))
    }
  }
}

velocityPlot <- function(file_name, domain, range, standAlone, mciz){
  scale <- 1.25
  dfOriginal <- read.csv(file_name, header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE)
  if (nrow(dfOriginal)>1){
    df <- filter(dfOriginal, cap==mciz)
    totalPfs <- nrow(df)
    x <- seq(domain[1], domain[2], 1)

    if (standAlone & nrow(df)>1){
      
      hist(df$center, breaks=x, col=colorPalette(2, 0.4), add=FALSE, freq=TRUE,  ylab='', xlab='', cex.axis=scale*1.25, main='')
      #hist(-1*df$top, breaks=seq(domain[1], domain[2], 1), col=colorPalette(3, 0.4), add=TRUE)
      legend('topright', legend=c('Bottom', 'Top'), col=c(colorPalette(2, 0.7), colorPalette(3, 0.7)), cex=1.5, pch=15)
      print(sapply(df, mean))
      print(sapply(df, sd))
      totalPfs <- nrow(df)
      df <- filter(df, top>3)
      df <- filter(df, bottom>3)
      print(sapply(df, mean))
      print(sapply(df, sd))
      mobilePfs <- nrow(df)
      print(paste('Mobile Pfs:', mobilePfs/totalPfs))
      
      grid()
    } else{
      hist(df$center, breaks=seq(domain[1], domain[2], 1), col=colorPalette(3, 0.4), freq=FALSE, ylim=range, xlab="", ylab="", main="", cex.lab=scale, cex.lab=scale, cex.axis=scale, cex.main=scale, cex.sub=2, xaxt='n', yaxt='n', axes=FALSE)
      
      
    }
    if (FALSE){
      mciz <- unique(dfOriginal$cap)
      for (ii in c(1:length(mciz))){
        df <- filter(dfOriginal, cap==mciz[ii])
        df <- filter(df, top>3)
        df <- filter(df, bottom>3)
        print(sapply(df, mean))
      }
    }
    

    
  }
  if (standAlone){
    mtext(expression('Velocity (subunits/s)'), side=1, line=3, cex=scale*1.25)
    mtext(expression('Frequency'), side=2, line=3, cex=scale*1.25)
  }
  
}


velocityLoss <- function(file_name){
  loss <- 100
  if (file.exists(file_name)){
    pfVelocity <- read.csv(file_name, header=FALSE, sep=',', dec='.', stringsAsFactors=FALSE)
    pfVelocity <- as.numeric(pfVelocity)
    mean <- round(mean(pfVelocity), 2)
    loss <- abs(6.5 - mean)
  }
  return(loss)
}



#########################################################################################################
##########      PLOT GEOMETRY #########################

sweepSubplots <- function(margin_top, margin_left, ii, jj, rows, cols, domain, range){
  len <- dividePlotVertical(rows, margin_top)
  wid <- dividePlotHorizontal(cols, margin_left)  
  label <- cols*(ii-1)+jj+96
  label <- paste0(chr(label), ')')
  borderColor <- colorPalette(7, .5)
  margins <- c(margin_left+(jj-1)*wid, margin_left+jj*wid, 1-margin_top-(ii)*len, 1-margin_top-len*(ii-1))
  par(fig=margins, new=TRUE, mar=c(0, 0, 0, 0))
  plot(seq(domain[1], domain[2], length=100), seq(range[1], range[2]/0.93, length=100), col="white", type = "n", ylab='', xlab='', cex.lab=1.5, xaxt='n', yaxt='n', axes=FALSE)
  box(col=borderColor)
  text(domain[1]+1, range[2]*0.9, label)
  grid()
  if (jj==cols){
    par(mgp = c(0, 1.5, 1))
    axis(4, cex.axis=1.25, tck=-0.02) #, at=round(seq(range[1], range[2], length=5), digits=2))
  }
  if (ii==rows){
    par(mgp = c(0, 1.5, 1))
    axis(1, cex.axis=1.25, tck=-0.02)
  }
}

sweepSubplotsBoxplots <- function(margin_top, margin_left, ii, jj, rows, cols, domain, range){
  len <- dividePlotVertical(rows, margin_top)
  wid <- dividePlotHorizontal(cols, margin_left)  
  label <- cols*(ii-1)+jj+96
  label <- paste0(chr(label), ')')
  borderColor <- colorPalette(7, .5)
  margins <- c(margin_left+(jj-1)*wid, margin_left+jj*wid, 1-margin_top-(ii)*len, 1-margin_top-len*(ii-1))
  par(fig=margins, new=TRUE, mar=c(0, 0, 0, 0))
  c3 <- 'white'
  boxplot(c(0, 0), xlim=c(0.5, length(domain)+0.5), ylim=range, medcol=c3, whiskcol=c3, staplecol=c3, boxcol=c3, outcol=c3, type = "n", ylab='', cex.lab=1.5, yaxt='n', axes=FALSE)
  box(col=borderColor)
  grid()
  if (jj==cols){
    par(mgp = c(0, 1, 1))
    axis(4, cex.axis=0.75, tck=-0.02)
  }
  if (ii==rows){
    par(mgp = c(0, 1, 1))
    axis(1, cex.axis=0.75, tck=-0.02, at=c(1, 2, 3), labels=domain)
  }
}
dividePlotHorizontal <- function(cols, margin_left){
  wid <- (1-margin_left*2)/cols
  return(wid)
}

dividePlotVertical <- function(rows, margin_top){
  len <- (1-margin_top*2)/rows
  return(len)
}


sweep2dLabels <- function(margin_top, margin_left, parameter1, parameter2, label1, label2, labelX, labelY){
  rows <- length(parameter1)
  cols <- length(parameter2)
  len <- dividePlotVertical(rows, margin_top)
  wid <- dividePlotHorizontal(cols, margin_left)
  fillColor <- colorPalette(8, 1)
  borderColor <- colorPalette(7, .5)

  colorRow <- 1
  margins <- c(0, margin_left/2, margin_top, 1-margin_top)
  par(fig=margins, new=TRUE, mar=c(0, 0, 0, 0))
  newEmptyPlot()
  
  if (rows>1){
    color <- colorPalette(colorRow, 0.5)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = color, border=color)
    text(0.5, 0.5, label1, srt=90, adj=0.5, cex=2)  
    for (ii in c(1:rows)){
      margins <- c(margin_left/2, margin_left, 1-margin_top-(ii)*len, 1-margin_top-len*(ii-1))      
      par(fig=margins, new=TRUE, mar=c(0, 0, 0, 0))
      newEmptyPlot()
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = fillColor, border=borderColor)
      text(0.5, 0.5, parameter1[ii], srt=90, cex=1.25) 
    }
  }

  if (cols>1){
    colorRow <- 6
    margins <- c(margin_left, margin_left+cols*wid, 1-margin_top/2, 1)
    par(fig=margins, new=TRUE, mar=c(0, 0, 0, 0))
    newEmptyPlot()
    color <- colorPalette(colorRow, 0.5)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = color, border=color)
    text(0.5, 0.5, label2, adj=0.5, cex=2)  
    for (jj in c(1:cols)){
      margins <- c(margin_left+jj*wid-wid, margin_left+jj*wid, 1-margin_top, 1-margin_top/2)
      par(fig=margins, new=TRUE, mar=c(0, 0, 0, 0))
      newEmptyPlot()
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = fillColor, border=borderColor)
      text(0.5, 0.5, parameter2[jj], cex=1.25)
    }
  }
  
  par(fig=c(margin_left, 1-margin_left, margin_top, 1-margin_top))
  mtext(labelX, side=1, line=3, cex=1.5) 
  mtext(labelY, side=4, line=3, cex=1.5)

}


colorPalette <- function(colorRow, alpha){
  fileName <- '/home/lauren/Documents/research/treadmilling-model/figures/colorPalette.csv'
  palette <- read.csv(fileName, header=TRUE)
  palette$r <- palette$r/255
  palette$g <- palette$g/255
  palette$b <- palette$b/255
  if (alpha<1.01){
    color <- rgb(palette[colorRow,]$r, palette[colorRow,]$g, palette[colorRow,]$b, alpha=alpha)
  } else{
    colorScaler <- 2 - alpha 
    color <- rgb(palette[colorRow,]$r*colorScaler, palette[colorRow,]$g*colorScaler, palette[colorRow,]$b*colorScaler, alpha=1)

  }
  return(color)
}

generateColorPalette <- function(totalColors, alpha){
  color <- c()
  for (ii in c(1:totalColors)){
    color <- c(color, colorPalette(ii, alpha))
  }
  return(color)
}

addUnits <- function(inputs, units){
  label <- c()

  if (units==1){ #'uM'
    for (input in inputs){
      label <- c(label,  as.expression(bquote(.(input)~mu*'M')))
    }
  } else if (units==2){ #s^-1
    for (input in inputs){
      label <- c(label,  as.expression(bquote(.(input)~'s'^-1)))
    }
  } else if (units==3){ #s^-1 uM^-1
    for (input in inputs){
      label <- c(label,  as.expression(bquote(.(input)~'s'^-1~mu*'M'^-1)))
    }
  } else if (units==4){ #s
    for (input in inputs){
      label <- c(label,  as.expression(bquote(.(input)~'s')))
    }
  }

  return(label)
}

newEmptyPlot <- function(){
  plot(c(0,1), c(0, 1), type = "n", ylab='', xlab='', cex.lab=1.5, axes=FALSE)
}

criticalConcentration <- function(ftsz,  magnitude, color){
  span <- 4*length(which(ftsz==ftsz[1]))
  if (length(ftsz)<(span+2)){
    ind <- c(2:length(ftsz))
  } else{
    r2 <- numeric((length(ftsz)-span))
    for (ii in c(1:(length(ftsz)-span))){
      ind <- c(ii:(ii+span))
      r2[ii] <- summary(lm(magnitude[ind]~ftsz[ind]))$r.squared
    }
    ind <- c(which.max(r2):(which.max(r2)+span)) 
  }
  if (length(ind)<4){
    ind <- c(1:length(ftsz))
  }
  #ind <- which(ftsz>0.5)
  if (sum(ind)>1){
    ftszInd <- ftsz[ind]
    model <- lm(magnitude[ind]~ftszInd)
    coef <- model$coefficients
    ftszCc <- c(min(ftsz)-1, ftszInd, max(ftsz)+1)
    fit <- ftszCc*coef[2] + coef[1]
    lines(ftszCc, fit, col=color, lwd=3)
    print(paste0('Slope: ', coef[2]))
    print(paste0('Critical Concentration: ', -coef[1]/coef[2]))
  }
}

subplot <- function(area){
  margins <- c(4, 4, 1, 1)
  options(warn=2)
  x <- try(mtext('bottomleft', ''), silent=TRUE)
  if(class(x)=="try-error"){
    par(fig=area, mar=margins, new=FALSE)
  } else{
    par(fig=area, mar=margins, new=TRUE)
  }

}

addLabel <- function(label){

  legend("topleft", paste0(label, ')'), bty="n", inset=c(-0.13, 0), cex=1)

}