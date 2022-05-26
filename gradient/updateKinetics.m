function Kinetics = updateKinetics(Kinetics, names, ak, gk)
  kinetics = pullKinetics(Kinetics, names);
  for ii=1:length(names)
    Kinetics.(names{ii}) = kinetics(ii) - ak(ii).*gk(ii);
  end
end
