clear all; close all;

addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');

appendFile = false; 

knuc = [5000, 10000];
khyd = [0.3, 0.5];
fragment = [0.0001];
anneal = [0.1, 0.5];


mciz = [0, 2];
ftsz = [3:0.5:6];

repeats =  1;
totalTime = 60;

iteration = 0;

for ii=1:length(knuc)
  for jj=1:length(khyd)
    for kk=1:length(fragment)
      for mm=1:length(anneal)
        for ll=1:length(mciz)

          runHydrolysis(knuc(ii), khyd(jj), fragment(kk), anneal(mm), mciz(ll), totalTime, repeats, ftsz, appendFile);
          iteration = iteration + 1
        end
      end
    end
  end
end
