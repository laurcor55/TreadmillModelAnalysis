

function loss = ssHyd(Kinetics)
  totalTime = 300; 

  Parameters = parameters();
  Parameters.totalTime = totalTime;
  Parameters.concCap = 0;
  Parameters.mixPFs = 0;
  ftsz = [3, 4];
  hydSsDiff = zeros(1, length(ftsz));
  hydSs = zeros(1, length(ftsz));
  hydSsExp = [4.3, 6.8]
  for ii=1:length(ftsz)
    Parameters.concTotalFtsZ = ftsz(ii);
    Outputs = runExperiment(Parameters, Kinetics, false);
    ind3 = length(Outputs.time);
    [~, ind2] = min(abs(Outputs.time - 270));
    [~, ind1] = min(abs(Outputs.time - 60));
    [~, ind0] = min(abs(Outputs.time - 30));
    dt01 = Outputs.time(ind1) - Outputs.time(ind0);
    dt23 = Outputs.time(ind3) - Outputs.time(ind2);

    hydSs01 = num2conc(Outputs.hydrolysisCount(ind1) - Outputs.hydrolysisCount(ind0))*60/dt01;
    hydSs23 = num2conc(Outputs.hydrolysisCount(ind3) - Outputs.hydrolysisCount(ind2))*60/dt23;
    hydSsDiff(ii) = abs(hydSs01 - hydSs23);
    hydSs(ii) = abs(hydSs23 - hydSsExp(ii));
  end

  loss = mean(hydSsDiff) + mean(hydSs);
  loss = loss./3;
end