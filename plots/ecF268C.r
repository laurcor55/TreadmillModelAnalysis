print('_____________________________')
print('ecF268C')
subplot(c(0, 0.5, 0.67, 1))
fileName <- 'assemblySweep/results/ecF268C_array.csv'
ssPfConcPlot(fileName, colorPalette(1, 1), 1, TRUE)
addLabel('a')

subplot(c(0.5, 1, 0.67, 1))
fileName <- 'gtpTurnoverSweep/results/ecF268C.csv'
gtpTurnoverPlot(fileName, colorPalette(1, 1), 1, FALSE, TRUE)
fileName <- 'gtpTurnoverSweep/results/ecF268C_1.csv'
#gtpTurnoverPlot(fileName, colorPalette(1, 1), 1, FALSE, FALSE)
addLabel('b')

subplot(c(0, 0.5, 0.33, 0.67))
pfLengthPlot('pfLength/results/ecF268C.csv', 2)
addLabel('c')

subplot(c(0.5, 1, 0.33, 0.67))
fileName <- '/home/lauren/Documents/research/treadmilling-model/paper/3-methods/mixingDiagram.png'
diagramPlot(fileName)
addLabel('d')

subplot(c(0, 0.5, 0, 0.33))
fileName <- 'pfMixing/results/ecF268C.csv'
experimentalFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005a/chen2005_fig4a.csv'
pfMixingPlot(fileName, experimentalFileName, 9.3, 1, TRUE)
addLabel('e')

subplot(c(0.5,1, 0, 0.33))
fileName <- 'pfMixing/results/ecF268C_EDTA.csv'
experimentalFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005a/chen2005_fig4b.csv'
pfMixingPlot(fileName, experimentalFileName, 127, 0.4, TRUE)
addLabel('f')