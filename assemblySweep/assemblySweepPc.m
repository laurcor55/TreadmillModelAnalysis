clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/classes', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

ftsz = [0.25, 0.5, 1, 1.5, 2, 3, 5, 7];

totalTime = 60;
time = 0:0.1:totalTime;
  
resultsOverall = time';
Parameters = ParametersClass;
  Parameters.totalTime = totalTime;
  Parameters.concCap = 0;
  Parameters.mixPFs = 0;
Kinetics = KineticsClass;
  Kinetics.kswitch = 15000;
  Kinetics.indGtpaseRate = 0.5;
  Kinetics.kgdpexchange = 0.5;

  Kinetics.kbongtp = 10;
  Kinetics.kboffgtp = 1;
  Kinetics.ktongtp = 1;
  Kinetics.ktoffgtp = 0.1;

%  Kinetics.kbongdp = 10;
%  Kinetics.kboffgdp = 1;
%  Kinetics.ktongdp = 1;
%  Kinetics.ktoffgdp = 0.1;

fileName = 'results/assemblyPc_noPc.csv'
                  
for mm = 1:length(ftsz)
  Parameters.concTotalFtsZ = ftsz(mm);
  result = multipleRuns(Parameters, Kinetics, time, 1);
  resultsOverall = [resultsOverall, result];
end
csvwrite(fileName, resultsOverall);

function results = multipleRuns(Parameters, Kinetics, time, repeats)
  results = zeros(length(time), repeats);
  for jj=1:repeats
    monomerConc = zeros(1, length(time));
    Outputs = runExperiment(Parameters, Kinetics, false);
    rowInd = 1;
    for kk=1:length(time)
      distanceToTime = Outputs.time - time(kk);
      [~, ind] = min(abs(distanceToTime));
      results(rowInd, jj) = Outputs.monomerConc(ind) + Outputs.capDimerConc(ind);
      rowInd = rowInd + 1;
    end
  end
end