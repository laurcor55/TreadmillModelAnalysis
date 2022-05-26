clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

knuc = [1000, 5000];
khyd = [0.4, 0.6];
fragment = [0.0001];
anneal = [1, 5];


iteration = 0;

mciz = [0.1, 0.5];
appendFile = true;
repeats = 5;

for ii=1:length(knuc)
  for jj=1:length(khyd)
    for kk=1:length(fragment)
      for mm=1:length(anneal)

        [outputFileName, Kinetics] = createFileName(knuc(ii), khyd(jj), fragment(kk), anneal(mm));
        velocityMeasure(outputFileName, Kinetics, repeats, appendFile, mciz);
        iteration = iteration + 1
      end
    end
  end
end

function [outputFileName, Kinetics] = createFileName(knuc, khyd, fragment, anneal)

  Kinetics = kineticsBs();
    Kinetics.knuc = knuc;
    Kinetics.khyd = khyd;
    Kinetics.kfragment = fragment;
    Kinetics.kanneal = anneal;

  outputFileName = strcat('results/knuc', num2str(knuc),'_khyd', num2str(khyd),'_fragment', num2str(fragment), '_anneal', num2str(anneal), '.csv')
end