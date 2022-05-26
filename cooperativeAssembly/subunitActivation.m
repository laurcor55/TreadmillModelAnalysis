clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

Parameters = parameters();
  Parameters.totalTime = 60;
  Parameters.concCap = 0;
  Parameters.mixPFs = 0;
  Parameters.concTotalFtsZ = 3.12;
Kinetics = kineticsEcF268C();

nucleationMode = 2

if nucleationMode == 1
  %ftsz = [0.1:0.1:0.5, 1, 1.5, 2, 3];
  %Kinetics.kgdpexchange = 0.1;
  %resultsOverall = concentrationSweep(Parameters, Kinetics, ftsz);
  %fileName = strcat('results_activation/gdpexchange', num2str(Kinetics.kgdpexchange),'.csv');
  %csvwrite(fileName, resultsOverall); 

  gdpexchange = [0.01, 0.1, 1, 10];
  resultsOverall = gdpexchangeSweep(Parameters, Kinetics, gdpexchange);
  fileName = 'results_activation/gdpexchangesweep.csv';
  csvwrite(fileName, resultsOverall); 

elseif nucleationMode == 2
  resultsOverall = concentrationSweep(Parameters, Kinetics, ftsz);
  fileName = strcat('results_rtequilibrium/kswitch', num2str(Kinetics.kswitch),'.csv');
  csvwrite(fileName, resultsOverall); 

  %kswitch = [100, 1000, 10000, 100000];
  %resultsOverall = kswitchSweep(Parameters, Kinetics, kswitch);
  %fileName = 'results_rtequilibrium/kswitchsweep.csv';
  %csvwrite(fileName, resultsOverall); 

elseif nucleationMode == 3
  %ftsz = [0.1:0.1:0.5, 1, 1.5, 2, 3];
  %Kinetics.kgdpexchange = 2;
  %Kinetics.kswitch = 10000;
  %resultsOverall = concentrationSweep(Parameters, Kinetics, ftsz);
  %fileName = strcat('results_both/kswitch', num2str(Kinetics.kswitch), '_gdpexchange', num2str(Kinetics.kgdpexchange),'.csv');
  %csvwrite(fileName, resultsOverall); 

  kswitch = [100, 1000, 10000, 100000];
  resultsOverall = kswitchSweep(Parameters, Kinetics, kswitch);
  fileName = 'results_both/kswitchsweep.csv';
  csvwrite(fileName, resultsOverall); 
end



function resultsOverall = gdpexchangeSweep(Parameters, Kinetics, gdpexchange)
  time = [0:0.1:Parameters.totalTime]';
  resultsOverall = time;
  for ii=1:length(gdpexchange)
    Kinetics.kgdpexchange = gdpexchange(ii);
    result = multipleRuns(Parameters, Kinetics, time, 1);
    resultsOverall = [resultsOverall, result];
  end
  resultsOverall = [0, gdpexchange; resultsOverall];
end


function resultsOverall = kswitchSweep(Parameters, Kinetics, kswitch)
  time = [0:0.1:Parameters.totalTime]';
  resultsOverall = time;
  for ii=1:length(kswitch)
    Kinetics.kswitch = kswitch(ii);
    result = multipleRuns(Parameters, Kinetics, time, 1);
    resultsOverall = [resultsOverall, result];
  end
  resultsOverall = [0, kswitch; resultsOverall];
end

function resultsOverall = concentrationSweep(Parameters, Kinetics, ftsz)
  time = [0:0.1:Parameters.totalTime]';
  resultsOverall = time;multipleRuns
  for ii=1:length(ftsz)
    Parameters.concTotalFtsZ = ftsz(ii);
    result = multipleRuns(Parameters, Kinetics, time, 1);
    resultsOverall = [resultsOverall, result];
  end
end

function results = multipleRuns(Parameters, Kinetics, time, repeats)
  results = zeros(length(time), repeats);
  for jj=1:repeats
    monomerConc = zeros(1, length(time));
    Outputs = runExperiment(Parameters, Kinetics, false);
    rowInd = 1;
    for kk=1:length(time)
      distanceToTime = Outputs.time - time(kk); 
      [~, ind] = min(abs(distanceToTime));
      results(rowInd, jj) = Outputs.monomerConc(ind);
      rowInd = rowInd + 1;
    end
  end
end