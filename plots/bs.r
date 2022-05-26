print('_____________________________')
print('B. SUBTILUS')


subplot(c(0, 0.5, 0.67, 1))
startingFtsz <- c(4, 5, 6)
fileName <- 'assemblySweep/results/bs2.csv'
assemblyPlot(fileName, startingFtsz, TRUE)
fileName <- 'experimental/bisson2015/figS5a_conc.csv'
assemblyPlotExperimental(fileName, startingFtsz)
addLabel('a')


subplot(c(0.5, 1, 0.67, 1))
fileName <- 'assemblySweep/results/bsMcizArray2.csv'
ssPfConcPlot(fileName, colorPalette(2, 1), 2, TRUE)
fileName <- 'assemblySweep/results/bsArray2.csv'
ssPfConcPlot(fileName, colorPalette(1, 1), 1, FALSE)
label <- addUnits(c(0, 2), 1)
legend("topleft", legend=label, col=generateColorPalette(2, 1), lty=1, cex=0.75, title='MciZ', pch=c(1, 2), lwd=2, inset=c(0.2, 0))
addLabel('b')



subplot(c(0, 0.5, 0.33, 0.67))
fileName <- 'gtpTurnoverSweep/results/bsMciz4.csv' 
gtpTurnoverPlot(fileName, colorPalette(2, 1), 2, FALSE, TRUE)
fileName <- 'gtpTurnoverSweep/results/bs4.csv' 
gtpTurnoverPlot(fileName, colorPalette(1, 1), 1, FALSE, FALSE)
label <- addUnits(c(0, 2), 1)
legend("topleft", legend=label, col=generateColorPalette(2, 1), lty=1, cex=0.75, title='MciZ', pch=c(1, 2), lwd=2, inset=c(0.2, 0))
addLabel('c')


subplot(c(0.5, 1, 0.33, 0.67))
velocityFileName <- 'velocitySweep/results/bs3.csv'
velocityPlotMciz(velocityFileName, TRUE)
addLabel('d')

subplot(c(0, 0.5, 0, 0.33))
pfLengthPlotMciz('pfLength/results/bs.csv')
addLabel('e')

subplot(c(0.5, 1, 0, 0.33))
fileName <- 'velocitySweep/results/hydVsVelocity.csv'
gtpTurnoverVsVelocityPlot(fileName)
addLabel('f')