
function [velocityTop, velocityBottom, velocityCenter] = velocityCalculate(Outputs, timeMeasure)
  [~, ssInd] = min(abs(Outputs.time - (Outputs.time(end) - timeMeasure)));
  ssPfs = Outputs.savePfs{ssInd};
  finalPfs = Outputs.savePfs{end};
  ssLocations = Outputs.saveLocations{ssInd};
  finalLocations = Outputs.saveLocations{end};
  ssTime = Outputs.time(end) - Outputs.time(ssInd);

  minPfs = min([length(ssPfs), length(finalPfs)]);
  velocityBottom = [];
  velocityTop = [];
  velocityCenter = [];
  for ii=1:minPfs 
    ssLength = length(ssPfs{ii});
    finalLength = length(finalPfs{ii});
    if (ssLength>5 && finalLength>5)
      velocityBottom = [velocityBottom, double(finalLocations(ii) - ssLocations(ii))./ssTime];
      ssTopLocation = ssLocations(ii) - ssLength;
      finalTopLocation = finalLocations(ii) - finalLength;
      velocityTop = [velocityTop, double(finalTopLocation - ssTopLocation)./ssTime];
      centerSs =  ssLocations(ii) - ssLength/2;
      centerFinal =  finalLocations(ii) - finalLength/2;
      velocityCenter = [velocityCenter, double(centerFinal - centerSs)./ssTime];
    end
  end
end