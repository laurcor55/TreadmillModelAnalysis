function pfMixingMeasure(Kinetics, totalTime, outputFileName, premixTime)
  Parameters = parameters();
    Parameters.totalTime = totalTime;
    Parameters.concCap = 0;
    Parameters.mixPFs = 1;
    Parameters.concTotalFtsZ = 6;
    Parameters.initialRoundTime = premixTime;


  fid = fopen(outputFileName, 'w');
  timeStep = Parameters.totalTime./50;
  time = 0:5:Parameters.totalTime;
  printTime(fid, time);

  for oo=1:5
    oo
    Outputs = runExperiment(Parameters, Kinetics, false);
    printResult(fid, time, Outputs);
  end
  fclose(fid);
  
end



function printTime(fid, time)
  fprintf(fid, "0");
  for oo=2:length(time)
    fprintf(fid, ",%.2f", time(oo));
  end
  fprintf(fid, "\n");
end

function printResult(fid, time, Outputs)
  pfMixingKinetics = calculateMixing(Outputs.savePfs);
  fprintf(fid, "0");
  for qq=2:length(time)
    [~, ind] = min(abs(Outputs.time - time(qq)));
    fprintf(fid, ",%.2f", pfMixingKinetics(ind));
  end
  fprintf(fid, "\n");
end