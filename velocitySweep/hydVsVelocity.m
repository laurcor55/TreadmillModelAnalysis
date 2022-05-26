clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

fileName = 'results/hydVsVelocity.csv';
fid = fopen(fileName, 'w');
fprintf(fid, 'hyd,vel\n');

gtpase = [0:0.1:0.5];
for ii=1:length(gtpase)
  runHydrolysisComparison(gtpase(ii), fid);
end

fclose(fid);

 

function runHydrolysisComparison(gtpase, fid)
  for jj=1:5
    jj
    gtpase
    Parameters = parameters();
      Parameters.totalTime = 60;
      Parameters.concTotalFtsZ = 10;
      Parameters.concCap = 0;
      Parameters.mixPFs = 0;
    Kinetics = kineticsBs();
      Kinetics.khyd = gtpase;
    
    Outputs = runExperiment(Parameters, Kinetics, false);
    [~, velocity, ~] = velocityCalculate(Outputs, 5);
    for ii = 1:length(velocity)
      fprintf(fid, '%f,%f\n', gtpase, velocity(ii));
    end
  end
end