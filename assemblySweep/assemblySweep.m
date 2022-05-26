clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/classes', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

ftsz = [3.12, 2.42, 1.56];  %[4:6]; %  [1.7, 2.2, 3]%      %% % 

kswitch = [9000];
gtpase = [0.44];
topoff = [6.5];
gdpexchange = [1];
bottomon = [10];
bottomoff = [5];
kcaponpf = 1.4;
mciz = [0];

parameterSweep(ftsz, kswitch, gtpase, topoff, gdpexchange, kcaponpf, bottomon, bottomoff, mciz);


function parameterSweep(ftsz, kswitch, gtpase, topoff, gdpexchange, kcaponpf, bottomon, bottomoff, mciz)
  totalTime = 15;
  repeats = 5;
  time = 0:1:totalTime;
  for ii=1:length(kswitch)
    for jj=1:length(gtpase)
      for kk=1:length(topoff)
        for ll=1:length(gdpexchange)
          for nn=1:length(kcaponpf)
            for pp=1:length(bottomon)
              for qq=1:length(bottomoff)
                for oo=1:length(mciz)
                  resultsOverall = [0; time'];
                  Parameters = ParametersClass;
                    Parameters.totalTime = totalTime;
                    Parameters.concCap = mciz(oo);
                    Parameters.mixPFs = 0;
                  Kinetics = kinetics();

                    Kinetics.kswitch = kswitch(ii);
                    Kinetics.indGtpaseRate = gtpase(jj);
                    Kinetics.kgdpexchange = gdpexchange(ll);
                    
                    Kinetics.kbongtp = bottomon(pp);
                    Kinetics.kboffgtp = bottomoff(qq);
                    Kinetics.kbongdp = 0.001;
                    Kinetics.kboffgdp = bottomoff(qq);
                    Kinetics.ktoffgdp = topoff(kk);

                    Kinetics.ktongtp = 1;
                    Kinetics.ktoffgtp = Kinetics.ktongtp./10.*bottomoff(qq);

                    Kinetics.kcaponpf = kcaponpf(nn);
                    Kinetics.ktongdp = topoff(kk)./1000;

                  if length(ftsz)>4
                    addedStuff = '_array';
                  else
                    addedStuff = '';
                  end
                  outputFileName = strcat('results/kswitch', num2str(kswitch(ii)),'_gtpase', num2str(gtpase(jj)), '_topoff', num2str(topoff(kk)), '_caponpf', num2str(kcaponpf(nn)), '_gdpexchange', num2str(gdpexchange(ll)),'_bottomon', num2str(bottomon(pp)), '_bottomoff', num2str(bottomoff(qq)),'_mciz', num2str(mciz(oo)), addedStuff, '.csv')
                  assemblyMeasure(ftsz, mciz, totalTime, repeats, outputFileName);
                end
              end
            end
          end
        end
      end
    end
  end
end

