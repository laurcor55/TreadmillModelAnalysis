clear all; close all;

addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/classes', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');



ftsz = 3;
ratio = [1, 2, 4, 8];
gtpase = [0.2, 0.6];


for ii=1:length(ratio)
    for kk=1:length(gtpase)
      runHydrolysisGdp(ratio(ii), gtpase(kk), ftsz, 3);
    end
end


function runHydrolysisGdp(ratio,  gtpase, ftsz, iterations)
  totalTime = 240;
  Parameters = ParametersClass;
    Parameters.totalTime = totalTime;
    Parameters.concCap = 0;
    Parameters.mixPFs = 0;
  Kinetics = KineticsClass;
    
    Kinetics.kboffgdp = Kinetics.ktoffgdp./ratio;
    Kinetics.kbongdp = Kinetics.kboffgdp./100;

    Kinetics.ktongtp = 10./ratio;
    Kinetics.ktoffgtp = Kinetics.ktongtp./3.3;


    Kinetics.indGtpaseRate = gtpase;
    Kinetics.kswitch = 30000;
  fileName = strcat('results/ratio', num2str(ratio), '_gtpase', num2str(gtpase), '.csv')
  outputTime = 0:1:totalTime;
  %if ~(isfile(fileName))
    fid = fopen(fileName, 'wt');
    fprintf(fid, '0');
    for kk=2:length(outputTime)
      fprintf(fid, ', %f', outputTime(kk));
    end
    fprintf(fid, '\n');
 % else 
   % fid = fopen(fileName, 'a');
 % end
  
  
  for ii=1:length(ftsz)
    for jj=1:iterations
      Parameters.concTotalFtsZ = ftsz(ii);
      Outputs = runExperiment(Parameters, Kinetics, false)
      fprintf(fid, '%f', ftsz(ii));
      for kk=2:length(outputTime)
        [~, ind1] = min(abs(Outputs.time - outputTime(kk)));
        [~, ind0] = min(abs(Outputs.time - outputTime(kk-1)));
        dt = Outputs.time(ind1) - Outputs.time(ind0);
        hyd = num2conc(Outputs.hydrolysisCount(ind1) - Outputs.hydrolysisCount(ind0))*60/dt;
        fprintf(fid, ', %f', hyd);
      end
      fprintf(fid, '\n');
    end
  end
  fclose(fid)
end