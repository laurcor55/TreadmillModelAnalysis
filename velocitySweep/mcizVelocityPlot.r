library(beeswarm)
library(dplyr)
caponpf <- c(10)
bottomoffcapgdp <- c(10)

kswitch <-5000
topongtp <- c(0.3)


cols <- length(bottomoffcapgdp)
rows <- length(caponpf)

outputFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/velocitySweep/mcizVelocity.png'
png(outputFileName, width=220*cols, height=220*rows)


parameterSweepPath <- getwd()
setwd("../")
source("figurePlotter.r")

mcizVelocityPlot <- function(parameterSweepPath, topongtp, bottomoffcapgdp, caponpf, kswitch){
  file_name = paste0(parameterSweepPath,'/results/topongtp', topongtp, '_bottomoffcapgdp', bottomoffcapgdp, '_caponpf', caponpf, '_kswitch', kswitch, '.csv')
  print(file_name)
  if (file.exists(file_name)){
    df <- read.csv(file_name, header=TRUE, sep=',', dec='.', stringsAsFactors=FALSE)
    if (nrow(df)>1){
      scale <- 1.5
      #df$bottom <- abs(df$bottom)
      df <- filter(df, bottom>1)
      boxplot(center~mciz, data=df, outline=FALSE, add=TRUE, xlab="", ylab="", type = "n", xaxt='n', yaxt='n', axes=FALSE)
      print(aggregate(df[, 2:4], list(df$mciz), mean))
      #if (nrow(df)<40){
        beeswarm(bottom~mciz, data=df, cex=0.25, method="hex", add=TRUE, xlab="", ylab="", type = "n", xaxt='n', yaxt='n', axes=FALSE)
     #}
    
    }
    
  }
}



margin_top <- 0.25/rows
margin_left <- 0.25/cols

plot.new()

domain <- c(0.1, 0.2, 0.4)
range <- c(0, 20)

for (ii in c(1:rows)){
  for (jj in c(1:cols)){
    sweepSubplotsBoxplots(margin_top, margin_left, ii, jj, rows, cols, domain, range)
    mcizVelocityPlot(parameterSweepPath, topongtp, bottomoffcapgdp[jj], caponpf[ii], kswitch)

  }
}
bottomoffcapgdp <- addUnits(bottomoffcapgdp, 2)
sweep2dLabels(margin_top, margin_left, caponpf, bottomoffcapgdp, expression(k[cap]^on), expression(k[cap]^off), expression(MciZ~(mu*M)), 'Velocity (subunits/s)')

dev.off()
