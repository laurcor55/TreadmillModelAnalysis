clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');


knuc = [1000, 10000];
khyd = [0];
fragment = [0.0001, 0.00001];
anneal = [0.1, 1, 10];


iteration = 0;

for ii=1:length(knuc)
  for jj=1:length(khyd)
    for kk=1:length(fragment)
      for mm=1:length(anneal)

        [outputFileName, Kinetics] = createFileName(knuc(ii), khyd(jj), fragment(kk), anneal(mm));
        pfMixingMeasure(Kinetics, 550, outputFileName, 60);
        
        iteration = iteration + 1
      end
    end
  end
end

function [outputFileName, Kinetics] = createFileName(knuc, khyd, fragment, anneal)

  Kinetics = kineticsEcF268C_EDTA();
    Kinetics.knuc = knuc;
    Kinetics.khyd = khyd;
    Kinetics.kfragment = fragment;
    Kinetics.kanneal = anneal;

  outputFileName = strcat('results/knuc', num2str(knuc),'_khyd', num2str(khyd),'_fragment', num2str(fragment), '_anneal', num2str(anneal), '.csv')
end