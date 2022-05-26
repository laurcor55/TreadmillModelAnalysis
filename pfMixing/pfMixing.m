clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

Kinetics = kineticsEcF268C();
outputFileName = 'results/ecF268C.csv'
pfMixingMeasure(Kinetics, 50, outputFileName, 30);

