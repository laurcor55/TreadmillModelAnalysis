function results = multipleRuns(Parameters, Kinetics, time, repeats)
  results = zeros(length(time)+1, repeats);
  for jj=1:repeats
    jj
    Outputs = runExperiment(Parameters, Kinetics, false);
    rowInd = 2;
    results(1, jj) = Parameters.concTotalFtsZ;
    for kk=1:length(time)
      distanceToTime = Outputs.time - time(kk);
      [~, ind] = min(abs(distanceToTime));
      results(rowInd, jj) = Outputs.monomerConc(ind);
      rowInd = rowInd + 1;
    end
    results(end, jj) = Outputs.monomerConc(end);
  end
  
end