print('_____________________________')
print('ecL68W')

subplot(c(0, 0.25, 0, 1))
startingFtsz <- c(2, 4, 6)
fileName <- 'assemblySweep/results/ecL68W.csv'
assemblyPlot(fileName, startingFtsz, TRUE)
fileName <- 'experimental/chen2005b/hmk/fig6c_conc.csv'
assemblyPlotExperimental(fileName, startingFtsz)
addLabel('a')


subplot(c(0.25, 0.5, 0, 1))
startingFtsz <- c(1.7, 2.2, 3)
fileName <- 'assemblySweep/results/ecL68W_EDTA.csv'
assemblyPlot(fileName, startingFtsz, TRUE)
fileName <- 'experimental/chen2005b/mek/dataTransformed.csv'
assemblyPlotExperimental(fileName, startingFtsz)
addLabel('b')

#subplot(c(0.5, 0.75, 0, 1))
#pfLengthPlotTime('pfLength/results/ecL68W.csv')
#addLabel('c')

subplot(c(0.5, 0.75, 0, 1))
fileNames <- c('disassembly/results/ecL68W_EDTA_20preassemble.csv', 'disassembly/results/ecL68W_EDTA_240preassemble.csv')

times <- c(20, 240)
disassemblyPlot(fileNames, times, c(12, 95), TRUE)
addLabel('c')

subplot(c(0.75, 1, 0, 1))
fileNames <- c('disassembly/results/ecL68W_20preassemble.csv', 'disassembly/results/ecL68W_240preassemble.csv')
times <- c(20, 240)
disassemblyPlot(fileNames, times, c(8, 8), TRUE)
addLabel('d')

#subplot(c(0.75, 1, 0, 1))
#fileName <- 'velocitySweep/results/ecL68W.csv'
#velocityPlot(fileName, c(-10, 10), c(0, 1), TRUE, 0)
#addLabel('d')

