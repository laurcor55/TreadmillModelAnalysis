clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

ftsz = [1.56, 2.42, 3.12];% %[1.7, 2.2, 3]
totalTime = 15;
mciz = 0;
repeats = 5;
outputFileName = 'results/ecF268C.csv';
Kinetics = kineticsEcF268C();

assemblyMeasure(Kinetics, ftsz, mciz, totalTime, repeats, outputFileName);

