
function disassemblyMeasure(Kinetics, preassembleTime, repeats, outputFileName)
  
  
  Parameters = parameters();
    Parameters.totalTime = 40;
    Parameters.initialRoundTime = preassembleTime;
    Parameters.concCap = 0;
    Parameters.mixPFs = 0;
    Parameters.disassemblePfs = 1;
    Parameters.concTotalFtsZ = 6;
  time = 0:0.1:Parameters.totalTime;
  resultsOverall = [0; time'];

  result = multipleRuns(Parameters, Kinetics, time, repeats);
  resultsOverall = [resultsOverall, result];
  csvwrite(outputFileName, resultsOverall);

end