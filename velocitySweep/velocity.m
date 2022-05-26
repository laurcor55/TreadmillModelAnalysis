clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');


Kinetics = kineticsBs();


outputFileName = 'results/bs4.csv'

mciz = [2];
appendFile = true;
repeats = 5;
velocityMeasure(outputFileName, Kinetics, repeats, appendFile, mciz);