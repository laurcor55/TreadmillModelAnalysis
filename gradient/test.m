startTime = clock();
repetitions = 100000;
parfor ii=1:repetitions
  rand([100, 100]);
end
finishTime = clock();
totalTime = finishTime - startTime