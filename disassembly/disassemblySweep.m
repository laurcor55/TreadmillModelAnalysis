
clear all; close all;

addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');


knuc = [500];
khyd = [0];
fragment = [0.00001];
anneal = [1, 10];

preassembleTime = [20, 120];

iteration = 0;

for ii=1:length(knuc)
  for jj=1:length(khyd)
    for kk=1:length(fragment)
      for mm=1:length(anneal)
        for ll=1:length(preassembleTime)

          runDisassembly(knuc(ii), khyd(jj), fragment(kk), anneal(mm), preassembleTime(ll));
          iteration = iteration + 1
        end
      end
    end
  end
end

function runDisassembly(knuc, khyd, fragment, anneal, preassembleTime)

  Kinetics = kineticsEcL68W_EDTA();
    Kinetics.knuc = knuc;
    Kinetics.khyd = khyd;
    Kinetics.kfragment = fragment;
    Kinetics.kanneal = anneal;

  repeats = 3;
  outputFileName = strcat('results/knuc', num2str(knuc),'_khyd', num2str(khyd),'_fragment', num2str(fragment), '_anneal', num2str(anneal), '_preassemble', num2str(preassembleTime), '.csv')
  disassemblyMeasure(Kinetics, preassembleTime, repeats, outputFileName);

end