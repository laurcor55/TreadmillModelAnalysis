function runHydrolysis(knuc, khyd, fragment, anneal, mciz, totalTime, repeats, ftsz, appendFile)
  ftsz = ftsz + mciz;
  Parameters = parameters();
    Parameters.totalTime = totalTime;
    Parameters.concCap = mciz;
    Parameters.mixPFs = 0;
  Kinetics = kineticsBs();
    Kinetics.knuc = knuc;
    Kinetics.khyd = khyd;
    Kinetics.kfragment = fragment;
    Kinetics.kanneal = anneal;
  outputFileName = strcat('results/knuc', num2str(knuc),'_khyd', num2str(khyd),'_fragment', num2str(fragment), '_anneal', num2str(anneal), '_mciz', num2str(mciz), '.csv')
  calculateHydrolysis(Kinetics, ftsz, mciz, repeats, totalTime, outputFileName, appendFile);
end