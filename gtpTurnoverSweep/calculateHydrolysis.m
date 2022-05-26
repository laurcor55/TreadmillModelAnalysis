function calculateHydrolysis(Kinetics,ftsz, mciz, repeats, totalTime, outputFileName, appendFile)
  Parameters = parameters();
  Parameters.totalTime = totalTime;
  Parameters.concCap = mciz;
  if (appendFile) & (isfile(outputFileName))
    fid = fopen(outputFileName, 'a');
    writeFile = true;
  else
    writeFile = true;
    fid = fopen(outputFileName, 'wt');
    fprintf(fid, '0');
    for kk=5:5:Parameters.totalTime
      fprintf(fid, ', %f', kk);
    end
    fprintf(fid, '\n');
  end
  
  if (writeFile)
    
    for ii=1:length(ftsz)
      for jj=1:repeats
        Parameters.concTotalFtsZ = ftsz(ii);
        Outputs = runExperiment(Parameters, Kinetics, false)
        

        fprintf(fid, '%f', Parameters.concTotalFtsZ);
        for kk=5:5:Parameters.totalTime
          [~, ind1] = min(abs(Outputs.time - kk));
          [~, ind0] = min(abs(Outputs.time - kk + 5));
          dt = Outputs.time(ind1) - Outputs.time(ind0);
          hyd = num2conc(Outputs.hydrolysisCount(ind1) - Outputs.hydrolysisCount(ind0))*60/dt;
          fprintf(fid, ', %f', hyd);
        end
        fprintf(fid, '\n');
      end
    end
    fclose(fid)
  end
end