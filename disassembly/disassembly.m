clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

Kinetics = kineticsEcL68W();

preassembleTime = 20;
repeats = 5;
outputFileName = 'results/ecL68W_20preassemble.csv';
disassemblyMeasure(Kinetics, preassembleTime, repeats, outputFileName);

preassembleTime = 240;
outputFileName = 'results/ecL68W_240preassemble.csv';
disassemblyMeasure(Kinetics, preassembleTime, repeats, outputFileName);


