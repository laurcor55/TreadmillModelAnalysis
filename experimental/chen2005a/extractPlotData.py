import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def filterByNeighbors(originalPlot, neighborhoodSize, neighborThreshold):
  filteredPlot = np.zeros((originalPlot.shape[0], originalPlot.shape[1]))
  maxValue = np.max(np.max(originalPlot))
  intensityThreshold = maxValue*0.5

  for ii in range(neighborhoodSize, originalPlot.shape[1]-neighborhoodSize):
    jj = neighborhoodSize
    while (jj<originalPlot.shape[0]-neighborhoodSize):
      jj+=1
      totalNeighbors = 0
      if (originalPlot[jj, ii] > intensityThreshold):
        neighborhood = originalPlot[jj-neighborhoodSize:jj+neighborhoodSize+1, ii-neighborhoodSize:ii+neighborhoodSize+1]
        totalNeighbors = np.sum(np.sum(neighborhood))/maxValue
      if (totalNeighbors>neighborThreshold):
        filteredPlot[jj, ii] = 1
        jj+=2*neighborhoodSize
  return filteredPlot

def lowPassFilter(originalPlot, size):
  kernel = np.ones((size, size))
  output = signal.convolve2d(originalPlot, kernel, boundary='symm', mode='same')
  return output

def findPlotValues(plot):
  fluorescence = np.zeros((6, plot.shape[1]))
  timeIndex = 0
  while (timeIndex<plot.shape[1]):
    ftszIndex = 0
    amplitude = 0
    while ((ftszIndex<6) and (amplitude<plot.shape[0])):
      if (plot[amplitude, timeIndex]>0):
        fluorescence[ftszIndex, timeIndex] = plot.shape[0] - amplitude
        ftszIndex += 1
        amplitude += 20
      else:
        fluorescence[ftszIndex, timeIndex] = 0
      amplitude += 1
    timeIndex += 1
  return fluorescence

def lightSignalToConcMonomer(lightSignal, startConcentration, criticalConcentration):
  lightSignalConv = signal.convolve(lightSignal, [0.25, 0.25, 0.25, 0.25], mode='same')
  lightSignalConv[0:3] = lightSignal[0:3]
  lightSignalConv[-4:-1] = lightSignal[-4:-1]
  fluorPerConc = (np.mean(lightSignalConv[0:5]) - np.mean(lightSignalConv[-6:-1]))/(startConcentration - criticalConcentration)
  dataPoints = lightSignalConv.shape[0]
  concMonomer = np.zeros(dataPoints)
  concMonomer[0] = startConcentration
  for ii in range(1, dataPoints):
    fluorDiff = lightSignalConv[ii] - lightSignalConv[ii-1]
    concDiff = fluorDiff/fluorPerConc
    concMonomer[ii] = concMonomer[ii-1] + concDiff
  concMonomer[-1] = criticalConcentration
  return concMonomer
  

imOriginal = plt.imread('chen2005_fig4a.png')
im = imOriginal[:, :, 0]
im = np.multiply(im, -1)
im = np.subtract(im, np.min(im))
im = np.pad(im, 10, 'edge')

len = 2
wid = 2
plt.subplot(len, wid, 1)
plt.imshow(im)

filteredPlot0 = lowPassFilter(im, 3)
plt.subplot(len, wid, 2)
plt.imshow(filteredPlot0)

filteredPlot1 = filterByNeighbors(filteredPlot0, 2, 4)
filteredPlot1 = filteredPlot1[10:-10, 10:-10]
plt.subplot(len, wid, 3)
plt.imshow(filteredPlot1)

fluorescence = findPlotValues(filteredPlot1)
space = np.linspace(0, fluorescence.shape[1], num=fluorescence.shape[1])  
time = np.linspace(0, 16, num=fluorescence.shape[1])
plt.subplot(len, wid, 4)
#plt.imshow(imOriginal)
for ii in range(6):
  plt.scatter(space, fluorescence[ii,:], alpha=0.1)
#plt.show()


startFtsz = [3.12, 2.42, 1.98, 1.56, 1.22, 0.85]
concMonomer = np.zeros((fluorescence.shape))
plt.figure()
for ii in range(concMonomer.shape[0]):
  concMonomer[ii,:] = lightSignalToConcMonomer(fluorescence[ii,:], startFtsz[ii], 0.53)
  plt.plot(time, concMonomer[ii,:])
plt.show()
result = np.vstack((concMonomer, time))
result = np.rot90(result, k=3)

np.savetxt("chen2005_fig4a.csv", result, delimiter=",", fmt='%f')
