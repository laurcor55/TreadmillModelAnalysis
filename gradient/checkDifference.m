function meanLossDifference = checkDifference(iteration, lossHistory)
  iterationCheckSpan = 3;
  if (iteration>iterationCheckSpan+1)
    lossDifference = abs(lossHistory(2:end) - lossHistory(1:end-1));
    meanLossDifference = mean(lossDifference(end-iterationCheckSpan:end));
  end
end