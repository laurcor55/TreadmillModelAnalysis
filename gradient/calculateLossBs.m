function loss = calculateLossBs(Kinetics)
  repeats = 1;
  loss = zeros(1, repeats);
  for ii=1:repeats
    hydLoss0 = hydVsConc(Kinetics, 0, 2.5, 1.2);
    hydLoss2 =  hydVsConc(Kinetics, 2, 5.5, 3.5);
    %velLoss = velocityMciz(Kinetics);
    %ssHydLoss = ssHyd(Kinetics);
    loss(ii) =  hydLoss0 + hydLoss2;%+ velLoss; %ssHydLoss; %
  end
  loss = mean(loss)./10;
end
