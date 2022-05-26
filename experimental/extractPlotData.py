import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from skimage.io import imread
import math

def filterByNeighbors(originalPlot, neighborhoodSize, neighborThreshold, startFtsz):
  outputSignal1 = np.zeros((len(startFtsz), originalPlot.shape[1]))
  maxValue = np.max(np.max(originalPlot))
  intensityThreshold = 0.2
  outputImage = np.zeros((originalPlot.shape))
  skipSize = 2
  for ii in range(originalPlot.shape[1]-neighborhoodSize):
    jj = neighborhoodSize
    linesFound = 0
    while (linesFound < len(startFtsz)) and jj < originalPlot.shape[0]-neighborhoodSize:
      jj+=1
      totalNeighbors = 0
      if (originalPlot[jj, ii] > intensityThreshold):
        verticalNeighborhood = originalPlot[jj:jj+neighborhoodSize,ii]
        maxInd = np.argmax(verticalNeighborhood)
        totalNeighbors = np.sum(verticalNeighborhood>intensityThreshold)
        if (totalNeighbors>neighborThreshold):
          outputSignal1[linesFound, ii] = jj+maxInd
          outputImage[jj, ii] = 1
          linesFound += 1
          jj+=skipSize*neighborhoodSize+maxInd
  outputSignal2 = np.zeros((len(startFtsz), originalPlot.shape[1]))
  for ii in range(originalPlot.shape[1]-neighborhoodSize):
    jj = originalPlot.shape[0]
    linesFound = 0
    while (linesFound<len(startFtsz)) and jj > 0:
      jj+=-1
      totalNeighbors = 0
      if (originalPlot[jj, ii] > intensityThreshold):
        verticalNeighborhood = originalPlot[jj-neighborhoodSize:jj,ii]
        maxInd = np.argmax(verticalNeighborhood)
        totalNeighbors = np.sum(verticalNeighborhood>intensityThreshold)
        if (totalNeighbors>neighborThreshold):
          outputSignal2[linesFound, ii] = jj-neighborhoodSize + maxInd
          linesFound += 1
          jj += -1*(skipSize*neighborhoodSize+maxInd)
  outputSignal2 = np.flipud(outputSignal2)
  outputSignal = np.add(outputSignal1, outputSignal2)
  outputSignal = np.divide(outputSignal, 2)
  return outputSignal, outputImage

def thresholdSignal(originalPlot, intensityThreshold):
  filteredPlot = (originalPlot > intensityThreshold) * originalPlot
  return filteredPlot

def lowPassFilter(originalPlot, size):
  kernel = np.ones((size, size))
  kernel = np.divide(kernel, size*size)
  output = signal.convolve2d(originalPlot, kernel, boundary='symm', mode='same')
  return output


def lightSignalToConcMonomer(lightSignal, startFtsz, criticalConcentration):
  concMonomer = np.zeros((lightSignal.shape))

  for jj in range(lightSignal.shape[0]):
    fluorPerConc = (np.mean(lightSignal[jj, 0:50]) - np.mean(lightSignal[jj, -51:-1]))/(startFtsz[jj] - criticalConcentration)
    print(fluorPerConc)
    concMonomer[jj, 0:2] = startFtsz[jj]
    for ii in range(2, lightSignal.shape[1]):
      fluorDiff = lightSignal[jj, ii] - lightSignal[jj, ii-1]
      concDiff = fluorDiff/fluorPerConc
      concMonomer[jj, ii] = concMonomer[jj, ii-1] + concDiff
  return concMonomer


def preprocessImage(imageFile, neighborhoodSize):
  imOriginal = imread(imageFile, as_gray=True)

  im = imOriginal
  im = np.subtract(im, np.min(np.min(im)))
  im = np.multiply(-1, im)
  im = np.add(im, -1*np.min(np.min(im)))
  im = np.divide(im, np.max(np.max(im)))
  padSize = neighborhoodSize
  im = np.pad(im, padSize, 'edge')

  return im

def postprocessSignal(inputSignal, padSize):
  outputSignal = inputSignal
  kernelLen = 5
  kernel = np.divide(np.ones(kernelLen), kernelLen)
  for ii in range(inputSignal.shape[0]):
    outputSignal[ii,:] = np.convolve(inputSignal[ii,:], kernel, 'same')
  outputSignal = outputSignal[:, padSize:padSize-1]
  return outputSignal

def iterativeLpfThresholdPlot(im, filterRounds, lpfSize, threshold):
  rows = filterRounds + 1
  wid = 1
  plt.figure()
  plt.subplot(rows, wid, 1)
  plt.imshow(im)
  filteredPlot = im
  for ii in range(filterRounds):
    filteredPlot = lpfThreshold(filteredPlot, lpfSize, threshold)
    plt.subplot(rows, wid, ii+2)
    plt.imshow(filteredPlot)
  plt.show()
  return filteredPlot

def lpfThreshold(im, lpfSize, threshold):
  filteredPlot = im
  filteredPlot = lowPassFilter(filteredPlot, lpfSize)
  filteredPlot = thresholdSignal(filteredPlot, threshold)
  return filteredPlot

def plotResult(time, fluorescence, concMonomer, im):
  rows = 1
  wid = 2
  plt.figure()
  for ii in range(fluorescence.shape[0]):
    plt.subplot(rows, wid, 1)
    plt.plot(fluorescence[ii,:])
    plt.subplot(rows, wid, 2)
    plt.plot(time, concMonomer[ii, :])
  plt.subplot(rows, wid, 1)
  plt.imshow(im, aspect='auto')
  plt.show()

def exportResult(time, concMonomer, imageFile):
  result = np.vstack((concMonomer, time))
  result = np.rot90(result, k=3)
  fileName = imageFile[0:-4] + '_conc.csv'
  print(fileName)
  np.savetxt(fileName, result, delimiter=",", fmt='%f')

imageFile = 'chen2012/fig3a.png'
startFtsz = [6, 5, 4, 3, 2]
criticalConcentration = 1
totalTime = 13
filterRounds = 1
neighborhoodSize = 30
neighborhoodThreshold = 10

im = preprocessImage(imageFile, neighborhoodSize)


filteredPlot = iterativeLpfThresholdPlot(im, 3, 20, 0.1)

fluorescence, outputImage = filterByNeighbors(filteredPlot, neighborhoodSize, neighborhoodThreshold, startFtsz)
plt.imshow(outputImage)
plt.show()
#fluorescence = postprocessSignal(fluorescence)
fluorescence = fluorescence[:, neighborhoodSize:-neighborhoodSize]
concMonomer = lightSignalToConcMonomer(fluorescence, startFtsz, criticalConcentration)
time = np.linspace(0, totalTime, num=concMonomer.shape[1])
plotResult(time, fluorescence, concMonomer, im)
exportResult(time, concMonomer, imageFile)




