clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');
totalTime = 20;
Kinetics = kineticsBs();
time = [5, 15, 20];
Parameters = parameters();
  Parameters.totalTime = totalTime;
  Parameters.concTotalFtsZ = 3;
  
Outputs = runExperiment(Parameters, Kinetics, false);
outputFileNameStart = 'results/snapshot_';
matrixPFSnapshots = cellsTo3DMatrix(Outputs.savePfs, Outputs.saveLocations);
for ii=1:length(time)
  distanceToTime = Outputs.time - time(ii);
  [~, ind] = min(abs(distanceToTime));
  outputFileName = strcat(outputFileNameStart, num2str(time(ii)), '.csv')
  csvwrite(outputFileName, matrixPFSnapshots(:, :, ind))
end