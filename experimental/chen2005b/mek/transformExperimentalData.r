library(stats)

transformData <- function(lightSignal, startFtsz, criticalConcentration){
  signalLength <- length(lightSignal)
  kernelLength <- 2
  kernel <- array(-10, c(1, kernelLength))
  kernel <- kernel/kernelLength
  lightSignalTransform <- matrix(lightSignal[1], 1, signalLength+2)
  lightSignalTransform[2:(signalLength+1)] <- lightSignal
  lightSignalTransform[signalLength+2] <- lightSignalTransform[signalLength+1] 
  lightSignal <- convolve(lightSignalTransform, kernel, type="open")
  lightSignal <- lightSignal[kernelLength:(signalLength+1)]
  lightSignal <- lightSignal - lightSignal[signalLength]
  dlightdconc <- (lightSignal[1] - lightSignal[signalLength])/(startFtsz-criticalConcentration)
  concMonomer <- array(startFtsz, c(signalLength, 1))
  for (ii in c(2:length(lightSignal))){
    concMonomer[ii] <- (lightSignal[ii]-lightSignal[ii-1])/dlightdconc + concMonomer[ii-1]
  }
  return(concMonomer)
}

startFtsz <- c( 1.7, 2.2, 3.0)
concNames <- c('17', '22', '30')
for (ii in c(1:length(concNames))){
  fileName <- paste0('rawData/MEKdata', concNames[ii], '.csv')
  data <- read.csv(fileName)
  ind <- which.min(abs(data[,1] - 15))
  print(ind)
  data <- data[1:ind,]
  dataTransform <- data
  dataTransform[,2] <- transformData(data[,2], startFtsz[ii], 0.3)
  if (ii==1){
    dataAll <- data[,]
    dataAllTransform <- dataTransform[,]
  } else{
    dataAll <- cbind(dataAll, data[,2])
    dataAllTransform <- cbind(dataAllTransform, dataTransform[,2])

  }
}


plot(1, type="n", xlab="", ylab="", xlim=c(0, 30), ylim=c(0, 3.5))
for (ii in c(1:length(concNames))){
  points(dataAllTransform[,1], dataAllTransform[,ii+1])
}

firstRow <- c(0, startFtsz)
dataAllTransform <- rbind(firstRow, dataAllTransform)
write.table(dataAllTransform, 'dataTransformed.csv', row.names=FALSE, col.names=FALSE, sep=",")
