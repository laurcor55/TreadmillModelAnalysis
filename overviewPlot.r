source("figurePlotter.r")

doublePlot <- function(sourceFileName, outputFileName, dimensions){
  subplotWidth <- 3200
  subplotWidth <- 2.5

#  tiff(paste0(outputFileName, '.tiff'), width = dimensions[2]*subplotWidth, height = dimensions[1]*subplotWidth, units = "px", res = 300)

  tiff(paste0(outputFileName, '.tiff'), width = dimensions[2]*subplotWidth, height = dimensions[1]*subplotWidth, units = "in", res = 600)
  source(sourceFileName)
  dev.off()

  png(paste0(outputFileName, '.png'), width = dimensions[2]*subplotWidth, height = dimensions[1]*subplotWidth, units = "in", res = 100)
  source(sourceFileName)
  dev.off()
}

sourceFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/assemblySweep/assemblySnapshots.r'
outputFileName <- 'plots/figure5'
dimensions <- c(1, 3)
#doublePlot(sourceFileName, outputFileName, dimensions)


sourceFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/cooperativeAssembly/assemblyMechanisms.r'
outputFileName <- 'plots/figure6'
dimensions <- c(1, 4)
doublePlot(sourceFileName, outputFileName, dimensions)

sourceFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/plots/sulA.r'
outputFileName <- 'plots/figure10'
dimensions <- c(1, 1)
doublePlot(sourceFileName, outputFileName, dimensions)



sourceFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/plots/ecF268C.r'
outputFileName <- 'plots/figure7'
dimensions <- c(3, 2)
doublePlot(sourceFileName, outputFileName, dimensions)

sourceFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/plots/ecL68W.r'
outputFileName <- 'plots/figure8'
dimensions <- c(1, 4)
doublePlot(sourceFileName, outputFileName, dimensions)


sourceFileName <- '/home/lauren/Documents/research/treadmilling-model/figures/plots/bs.r'
outputFileName <- 'plots/figure9'
dimensions <- c(3, 2)
doublePlot(sourceFileName, outputFileName, dimensions)

