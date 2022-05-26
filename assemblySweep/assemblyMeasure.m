
function assemblyMeasure(Kinetics, ftsz, mciz, totalTime, repeats, outputFileName)
  
  time = 0:0.1:totalTime;
  resultsOverall = [0; time'];
  Parameters = parameters();
    Parameters.totalTime = totalTime;
    Parameters.concCap = mciz;
    Parameters.mixPFs = 0;

  for mm = 1:length(ftsz)
    Parameters.concTotalFtsZ = ftsz(mm);% + mciz;
    result = multipleRuns(Parameters, Kinetics, time, repeats);
    resultsOverall = [resultsOverall, result];
    csvwrite(outputFileName, resultsOverall);
  end

end