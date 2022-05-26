function velocityMeasure(outputFileName, Kinetics, repeats, appendFile, mciz)
  
  Parameters = parameters();
    Parameters.totalTime = 60;
    Parameters.concTotalFtsZ = 10;
    Parameters.mixPFs = 0;
  if (appendFile)
    fid = fopen(outputFileName, 'a');
  else
    fid = fopen(outputFileName, 'w'); 
    fprintf(fid, 'cap,bottom,center,top\n');             
  end
  
    
    
  for rr = 1:length(mciz)
    Parameters.concCap = mciz(rr)
    for oo=1:repeats
      Outputs = runExperiment(Parameters, Kinetics, false)
      [velocityTop, velocityBottom, velocityCenter] = velocityCalculate(Outputs, 5);
      for mm=1:length(velocityTop)
        fprintf(fid, '%f,%f,%f,%f\n', Parameters.concCap, velocityBottom(mm), velocityCenter(mm), velocityTop(mm));

      end
    end
  end
    
  fclose(fid);
end