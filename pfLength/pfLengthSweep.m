clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/classes', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

Parameters = ParametersClass;
  Parameters.totalTime = 300;
  Parameters.concCap = 0;
  Parameters.mixPFs = 0;

kswitch = [30000];
gtpase = [0.44];
topoff = 6.5;
gdpexchange = [1];
bottomon = [10];
bottomoff = 5;
topongtp = [1];
ftsz = [3];
mciz = [0, 0.3];

for ii=1:length(kswitch)
  for jj=1:length(gtpase)
    for kk=1:length(topoff)
      for ll=1:length(gdpexchange)
        for nn=1:length(bottomon)
          for pp=1:length(topongtp)
            
            for oo=1:length(ftsz)
              fileName = strcat('results/kswitch', num2str(kswitch(ii)),'_gtpase', num2str(gtpase(jj)), '_topoff', num2str(topoff(kk)), '_bottomon', num2str(bottomon(nn)), '_gdpexchange', num2str(gdpexchange(ll)),'_topongtp', num2str(topongtp(pp)), '_bottomoff', num2str(bottomoff), '_ftsz', num2str(ftsz(oo)), '.csv')
              fid = fopen(fileName, 'w');
              fprintf(fid, 'time,mciz,len\n');
              for rr = 1:length(mciz)
                  Parameters.concCap = mciz(rr);
                  Parameters.mixPFs = 0;
                  Parameters.concTotalFtsZ = ftsz(oo);
                Kinetics = kinetics();
                  Kinetics.kswitch = kswitch(ii);
                  Kinetics.indGtpaseRate = gtpase(jj);
                  Kinetics.ktoffgdp = topoff(kk);
                  Kinetics.kgdpexchange = gdpexchange(ll);
                  Kinetics.kbongtp = bottomon(nn);
                  Kinetics.ktongtp = topongtp(pp);

                  Kinetics.kboffgtp = bottomoff;
                  Kinetics.ktoffgtp = topoff(kk)./10.*bottomoff;

                  Kinetics.ktongdp = topoff(kk)./1000;

                for qq=1:5
                  singleRunPfLength(Parameters, Kinetics, fid)
                end
              end
            end
            fclose(fid);

          end
        end
      end
    end
  end
end

function singleRunPfLength(Parameters, Kinetics, fid)
  timePoints = [20, Parameters.totalTime];
  timeNames = {'Nucleation', 'Steady State'};

  Outputs = runExperiment(Parameters, Kinetics, false);
  for jj=1:length(timePoints)
    [~, ind] = min(abs(timePoints(jj)-Outputs.time));
    pfs = Outputs.savePfs{ind};
    for ii=1:length(pfs)
      if (length(pfs{ii})>1)
        fprintf(fid, '%s,%f,%f\n', timeNames{jj}, Parameters.concCap, length(pfs{ii}));
      end
    end
  end
end