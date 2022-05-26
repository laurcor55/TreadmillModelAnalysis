clear all; close all;

addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/data', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');
addpath('../velocitySweep');
addpath('../assemblySweep');
addpath('../pfMixing');
addpath('lossFunctions');
%names = {'kbongtp', 'kboffgtp', 'kbongdp', 'kboffgdp', 'ktongdp', 'ktoffgdp', 'ktongtp', 'ktoffgtp', 'capKd', 'kcaponpf', 'kcapoffpf', 'kcapoffpfgdp', 'khyd', 'kgdpexchange', 'knuc', 'kanneal'};

names = {'kbongtp', 'kboffgtp', 'kbongdp', 'kboffgdp',  'ktoffgdp', 'ktongtp', 'khyd', 'knuc', 'kanneal'};

Kinetics = kineticsBs();
%names = fieldnames(Kinetics);

kineticsValues = pullKinetics(Kinetics, names);
n = length(names);
totalIterations = 100;
A = totalIterations.*0.05;
alpha = 0.602;
gamma = 0.101;

minLoss = 100;

%pc = parcluster('local')
%parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')))

for k=1:totalIterations
  kineticValues = pullKinetics(Kinetics, names);
  a = kineticsValues.*0.1.*rand((size(kineticsValues)));
  c = a;

  ak = a;%/(A+k).^alpha;
  ck = c;%./k.^gamma;
  dk = (2*round(rand(n,1))-1);
  loss = zeros(1, 2);

  kineticsCell = perturbKinetics(Kinetics, names, ck, dk);
  for ii=1:2
    kineticsCell{ii}
    loss(ii) =  calculateLossBs(kineticsCell{ii});
  end
  gk = ((loss(2) - loss(1))./(2.*ck.*dk));
  Kinetics = updateKinetics(Kinetics, names, ak, gk)
  lossMean(k) = mean(loss)
  plot(lossMean)
  if (minLoss>lossMean(k))
    [~, ind] = min(loss);
    lowestKinetics = kineticsCell{ind};
    lowestLoss = loss;
  end
  lowestKinetics
end
