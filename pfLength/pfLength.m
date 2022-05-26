clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');
Kinetics = kineticsEcF268C();
ftsz = 2;
mciz = [0];
outputFileName = 'results/ecF268C.csv';

pfLengthMeasure(Kinetics, ftsz, mciz, outputFileName)