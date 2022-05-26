function kineticsCell = perturbKinetics(Kinetics, names, ck, dk)
  kinetics1 = Kinetics;
  kinetics2 = Kinetics;
  for ii=1:length(names)
    kinetics1.(names{ii}) = kinetics1.(names{ii}) - ck(ii).*dk(ii);
    kinetics2.(names{ii}) = kinetics2.(names{ii}) + ck(ii).*dk(ii);
  end
  kinetics1.ktoffgtp = kinetics1.ktongtp./kinetics1.kbongtp.*kinetics1.kboffgtp;
  kinetics2.ktoffgtp = kinetics2.ktongtp./kinetics2.kbongtp.*kinetics2.kboffgtp;
  

  kineticsCell = {kinetics1, kinetics2};
end
