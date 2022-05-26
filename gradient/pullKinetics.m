

function kinetics = pullKinetics(Kinetics, names)
  kinetics = zeros(length(names), 1);
  for ii=1:length(names)
    kinetics(ii) = Kinetics.(names{ii});
  end
end