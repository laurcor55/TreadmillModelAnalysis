clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/classes', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

ftsz = 6;

Parameters = ParametersClass;
  Parameters.totalTime = 60;
  Parameters.concTotalFtsZ = ftsz;
  Parameters.mixPFs = 0;
Kinetics = KineticsClass;
  Kinetics.indGtpaseRate = 0.3;
  Kinetics.kgdpexchange = 1;
  Kinetics.ktoffgdp = 6.5;
  Kinetics.kboffgtp = 1.5;
  Kinetics.kbongtp = 5;
  Kinetics.kbongdp = 0.001;
  Kinetics.kboffgdp = 6.5;

kswitch = [30000];
topongtp = [0.001, 1];
bottomoffcap = [1];
caponpf = [2];
mciz = [0, 0.5, 1];

for ii=1:length(topongtp)
  for kk=1:length(bottomoffcap)
    for ll=1:length(caponpf)
      for mm=1:length(kswitch)
        outputFilename = strcat('results/topongtp', num2str(topongtp(ii)), '_bottomoffcap', num2str(bottomoffcap(kk)), '_caponpf', num2str(caponpf(ll)), '_kswitch', num2str(kswitch(mm)), '.csv')
     %    if ~(isfile(outputFilename))
          fid = fopen(outputFilename, 'w');
          fprintf(fid, 'mciz,bottom,center,top\n');
    %    else 
    %      fid = fopen(outputFilename, 'a');
    %    end
        for oo=1:length(mciz)
          Kinetics.ktongtp = topongtp(ii);
          Kinetics.ktoffgtp = topongtp(ii)./3.3;
          Kinetics.kcaponpf = caponpf(ll);
          Kinetics.kcapoffpf = bottomoffcap(kk); 
          Kinetics.kswitch = kswitch(mm);

          Parameters.concCap = mciz(oo);
        
          
          for pp=1:2
            Outputs = runExperiment(Parameters, Kinetics, false);
            [velocityTop, velocityBottom, velocityCenter] = velocityMeasure(Outputs, 0.5);
            for nn=1:length(velocityTop)
              fprintf(fid, '%f,%f,%f,%f\n', mciz(oo), velocityBottom(nn), velocityCenter(nn), velocityTop(nn));
            end
          end
        end
        fclose(fid);
      end
    end
  end
end
