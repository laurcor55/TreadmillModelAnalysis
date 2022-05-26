print('_____________________________')
print('sulA')


margins <- c(3, 3, 1, 1)
par(mar=margins)


fileName <- 'gtpTurnoverSweep/results/cap3.csv' 
gtpTurnoverPlot(fileName, colorPalette(3, 1), 3, FALSE, TRUE)
fileName <- 'gtpTurnoverSweep/results/seq3.csv' 
gtpTurnoverPlot(fileName, colorPalette(2, 1), 2, FALSE, FALSE)
fileName <- 'gtpTurnoverSweep/results/ecF268C.csv'
gtpTurnoverPlot(fileName, colorPalette(1, 1), 1, FALSE, FALSE)
label <- c('No SulA', 'Sequesterer', 'Capper')
legend("topleft", legend=label, col=generateColorPalette(3, 1), lty=1, cex=0.75, pch=c(1, 2, 3), lwd=1.25, inset=c(0, 0))
