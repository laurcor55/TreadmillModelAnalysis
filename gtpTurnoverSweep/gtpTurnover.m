clear all; close all;

addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

appendFile = false; 
ftsz = [1:10];
Kinetics = kineticsEcF268C();
Kinetics.kcaponpf = 0
repeats = 2;
totalTime = 60;

mciz = 3;
outputFileName = 'results/seq3.csv';
calculateHydrolysis(Kinetics, ftsz, mciz, repeats, totalTime, outputFileName, appendFile);

mciz = 2;
outputFileName = 'results/bsMciz4.csv';
%calculateHydrolysis(Kinetics, ftsz, mciz, repeats, totalTime, outputFileName, appendFile);