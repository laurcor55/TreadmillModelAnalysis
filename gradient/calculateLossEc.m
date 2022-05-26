function loss = calculateLossEc(Kinetics)
  loss = 0;
  fileName = '/home/lauren/Documents/research/treadmilling-model/figures/experimental/chen2005b/hmk/fig6c_conc.csv'
  for ii=1:1
  %  loss = loss + pfMixing(Kinetics);
    loss = loss + assembly(Kinetics, 6, 5, fileName);
    loss = loss + assembly(Kinetics, 4, 5, fileName);
  % loss = loss + hydVsConc(Kinetics, 0, 5, 0.5);
  end
  loss = loss;
end