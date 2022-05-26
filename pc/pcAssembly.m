clear all; close all;

addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/classes', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');


totalTime = 60;
ftsz = [0.25, 0.5, 1, 1.5, 2, 3, 5, 7];
iterations = 1;
Parameters = ParametersClass;
  Parameters.totalTime = totalTime;
  Parameters.concCap = 0;
  Parameters.mixPFs = 0;
Kinetics = KineticsClass;
  Kinetics.kswitch = 1000;
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

fileExtension = 'stableT';

fileName = strcat('results/gtpTurnoverPc_', fileExtension, '.csv');
delete(fileName)

fileName = strcat('results/assemblyPc_', fileExtension, '.csv');
delete(fileName)


for ii=1:length(ftsz)
  for jj=1:iterations
    Parameters.concTotalFtsZ = ftsz(ii);
    Outputs = runExperiment(Parameters, Kinetics, false);
    writeGtpTurnover(Outputs, Parameters, fileExtension);
    writeAssembly(Outputs, Parameters, fileExtension);
  end
end


function writeGtpTurnover(Outputs, Parameters, fileExtension)
  fileName = strcat('results/gtpTurnoverPc_', fileExtension, '.csv');
  timeStep = round(Parameters.totalTime./2);

  if ~(isfile(fileName))
    fid = fopen(fileName, 'wt');
    fprintf(fid, '0');
    for kk=1:Parameters.totalTime
      fprintf(fid, ', %f', kk);
    end
    fprintf(fid, '\n');
  else 
    fid = fopen(fileName, 'a');
  end
  fprintf(fid, '%f', Parameters.concTotalFtsZ);
  for kk=1:Parameters.totalTime
    [~, ind1] = min(abs(Outputs.time - kk));
    [~, ind0] = min(abs(Outputs.time - kk + timeStep));
    dt = Outputs.time(ind1) - Outputs.time(ind0);
    hyd = num2conc(Outputs.hydrolysisCount(ind1) - Outputs.hydrolysisCount(ind0))*60/dt;
    fprintf(fid, ', %f', hyd);
  end
  fprintf(fid, '\n');
  fclose(fid);
end

function writeAssembly(Outputs, Parameters, fileExtension)
  fileName = strcat('results/assemblyPc_', fileExtension, '.csv');
  timeStep = round(Parameters.totalTime./2);
  if ~(isfile(fileName))
    fid = fopen(fileName, 'wt');
    fprintf(fid, '0');
    for kk=1:Parameters.totalTime
      fprintf(fid, ', %f', kk);
    end
    fprintf(fid, '\n');
  else 
    fid = fopen(fileName, 'a');
  end
  fprintf(fid, '%f', Parameters.concTotalFtsZ);
  for kk=1:Parameters.totalTime
    [~, ind1] = min(abs(Outputs.time - kk));
    fprintf(fid, ', %f', Outputs.pfFtszConc(ind1));
  end
  fprintf(fid, '\n');
  fclose(fid);
end