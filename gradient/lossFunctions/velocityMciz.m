

function loss = velocityMciz(Kinetics)
  Parameters = ParametersClass;
  Parameters.totalTime = 60;
  Parameters.concTotalFtsZ = 6;
  Parameters.mixPFs = 0;
  mciz = [0, 0.5];
  velocityMeanMoving = zeros(1, length(mciz));
  velocityMeanTotal = zeros(1, length(mciz));

  for ii=2:length(mciz)
    Parameters.concCap = mciz(ii);
    velocityCenter = [];
    for jj=1:3
      Outputs = runExperiment(Parameters, Kinetics, false);
      [~, ~, velocityCenterRun] = velocityCalculate(Outputs, 1);
      velocityCenter = [velocityCenter, velocityCenterRun];
    end
    velocityMeanTotal(ii) = mean(velocityCenter);
    ind = find(velocityCenter>3);
    velocityMeanMoving(ii) = mean(velocityCenter(ind));
  end
  loss = 0;
  %loss = abs(velocityMeanMoving(1)-6.5);
  loss = loss + abs(velocityMeanMoving(2)-11);
end