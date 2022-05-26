function pfLengthMeasure(Kinetics, ftsz, mciz, outputFileName)
  fid = fopen(outputFileName, 'w');
  fprintf(fid, 'time,mciz,len\n');

  Parameters = parameters();
      Parameters.mixPFs = 0;
      Parameters.concTotalFtsZ = ftsz;
      Parameters.totalTime = 4*60;
  repeats = 5;
  for ll=1:repeats
    for kk=1:length(mciz)
      Parameters.concCap = mciz(kk);
      
      Outputs = runExperiment(Parameters, Kinetics, false);
      
      timePoints = [20, 60, Parameters.totalTime];
      for jj=1:length(timePoints)
        [~, ind] = min(abs(timePoints(jj)-Outputs.time));
        pfs = Outputs.savePfs{ind};
        for ii=1:length(pfs)
          if (length(pfs{ii})>1)
            fprintf(fid, '%i,%f,%f\n', round(timePoints(jj)), Parameters.concCap, length(pfs{ii}));
          end
        end
      end
    end
  end
  fclose(fid);
  
end