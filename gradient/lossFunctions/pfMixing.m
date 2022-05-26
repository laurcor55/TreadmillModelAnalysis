function loss = pfMixing(Kinetics)
  Parameters = ParametersClass;
    Parameters.totalTime = 30;
    Parameters.concCap = 0;
    Parameters.mixPFs = 1;
    Parameters.concTotalFtsZ = 6;
  time = 0:1:Parameters.totalTime;
  repeats = 1;
  mixedProportion = zeros(length(time), 1);
  for ii=1:repeats
    Outputs = runExperiment(Parameters, Kinetics, false);
    proportion = calculateMixing(Outputs.savePfs);
    for qq=2:length(time)
      [~, ind] = min(abs(Outputs.time - time(qq)));
      mixedProportion(qq) = mixedProportion(qq) + proportion(ind);
    end
  end
  mixedProportion = mixedProportion./repeats;
  experimentalFileName = '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005a/chen2005_fig4a.csv';
  dfExp = csvread(experimentalFileName);
  timeExperimental = dfExp(:,1);
  fluor = dfExp(:,2);
  fluor = -1*fluor + 1;
  mixedProportionExperimental = zeros(length(time), 1);
  for qq=1:length(time)
    [~, ind] = min(abs(timeExperimental - time(qq)));
    mixedProportionExperimental(qq) = fluor(ind);
  end
  loss = sum(abs(mixedProportionExperimental - mixedProportion));
end