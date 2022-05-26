

function loss = hydVsConc(Kinetics, mciz, slopeExpected, ccExpected)
  totalTime = 30; 

  Parameters = parameters();
  Parameters.totalTime = totalTime;
  Parameters.concCap = mciz;
  Parameters.mixPFs = 0;
  ftsz = [4, 4, 6, 6] + mciz;
  hydSs = zeros(1, length(ftsz));
  concSs = zeros(1, length(ftsz));

  for ii=1:length(ftsz)
    Parameters.concTotalFtsZ = ftsz(ii);
    Outputs = runExperiment(Parameters, Kinetics, false);

    [~, ind1] = min(abs(Outputs.time - totalTime));
    [~, ind0] = min(abs(Outputs.time - (totalTime).*3./4));
    dt = Outputs.time(ind1) - Outputs.time(ind0);
    hydSs(ii) = num2conc(Outputs.hydrolysisCount(ind1) - Outputs.hydrolysisCount(ind0))*60/dt;
    concSs(ii) = mean(Outputs.monomerConc(ind0:ind1));
  end
  fitLine = polyfit(ftsz, hydSs, 1);
  slope = fitLine(1);
  cc = -1.*fitLine(2)./fitLine(1);
  loss = abs(slopeExpected-slope) + 1.*abs(ccExpected - cc) + 1.*abs(ccExpected-mean(concSs));
end