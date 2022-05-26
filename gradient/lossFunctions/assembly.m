

function loss = assembly(Kinetics, ftsz, totalTime, experimentalFile)
  time = 0:totalTime;
  repeats = 2;
  Parameters = ParametersClass;
  Parameters.totalTime = totalTime;
  Parameters.concTotalFtsZ = ftsz;
  results = multipleRuns(Parameters, Kinetics, time, repeats);

  experimentalResults = csvread(experimentalFile);
  [~, columnInd] = min(abs(results(1, 2:end) - Parameters.concTotalFtsZ));
  loss = 0;
  for ii=1:length(time)
    [~, experimentalInd] = min(abs(time(ii) - experimentalResults(2:end, 1)));
    loss = loss + abs(mean(results(ii, :)) - experimentalResults(experimentalInd+1, columnInd+1));
  end
  loss = loss./length(time);
end