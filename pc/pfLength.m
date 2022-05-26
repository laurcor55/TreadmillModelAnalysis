clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/classes', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');
Parameters = ParametersClass;
  Parameters.totalTime = 60;
  Parameters.concCap = 0;
  Parameters.mixPFs = 0;
  Parameters.concTotalFtsZ = 3;
Kinetics = KineticsClass;
  Kinetics.kswitch = 100;
  Kinetics.indGtpaseRate = 0.5;
  Kinetics.kgdpexchange = 0.5;

  Kinetics.kbongtp = 10;
  Kinetics.kboffgtp = 1;
  Kinetics.ktongtp = 1;
  Kinetics.ktoffgtp = 0.1;

  Kinetics.kbongdp = 10;
  Kinetics.kboffgdp = 1;
  Kinetics.ktongdp = 1;
  Kinetics.ktoffgdp = 0.1;

fileExtension = 'stableT'; %'knuc15000_sameKinetics';

fileName = strcat('results/pfLengthPc_', fileExtension, '.csv');
delete(fileName)
fid = fopen(fileName, 'w');
fprintf(fid, 'time,len\n');

for qq=1:5
  singleRunPfLength(Parameters, Kinetics, fid)
end

fclose(fid);

    

function singleRunPfLength(Parameters, Kinetics, fid)
  timePoints = [20, Parameters.totalTime];
  timeNames = {'Nucleation', 'Steady State'};

  Outputs = runExperiment(Parameters, Kinetics, false);
  for jj=1:length(timePoints)
    [~, ind] = min(abs(timePoints(jj)-Outputs.time));
    pfs = Outputs.savePfs{ind};
    for ii=1:length(pfs)
      if (length(pfs{ii})>1)
        fprintf(fid, '%s,%f\n', timeNames{jj},  length(pfs{ii}));
      end
    end
  end
end