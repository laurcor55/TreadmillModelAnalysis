clear all; close all;
addpath('../../model/lib')
addpath('../../model/lib/analysis', '../../model/lib/classes', '../../model/lib/dependencies', '../../model/lib/export', '../../model/lib/gillespie');
addpath('../');
addpath('../velocitySweep');

ftsz = 6;

Parameters = ParametersClass;
  Parameters.totalTime = 60;
  Parameters.concTotalFtsZ = ftsz;
  Parameters.mixPFs = 0;
  Parameters.concCap = 0;
Kinetics = KineticsClass;
  Kinetics.kswitch = 100;
  Kinetics.indGtpaseRate = 0.5;
  Kinetics.kgdpexchange = 0.5;

  Kinetics.kbongtp = 10;
  Kinetics.kboffgtp = 1;
  Kinetics.ktongtp = 1;
  Kinetics.ktoffgtp = 0.1;

  Kinetics.kbongdp = 10;
  Kinetics.kboffgdp = 1;
  Kinetics.ktongdp = 1;
  Kinetics.ktoffgdp = 0.1;


outputFilename = 'results/velocityPc.csv'
fid = fopen(outputFilename, 'w');
fprintf(fid, 'cap,bottom,center,top\n');


for oo=1:5
  Outputs = runExperiment(Parameters, Kinetics, false);
  [velocityTop, velocityBottom, velocityCenter] = velocityMeasure(Outputs, 5);
  for mm=1:length(velocityTop)
    fprintf(fid, '%f,%f,%f,%f\n', Parameters.concCap, velocityBottom(mm), velocityCenter(mm), velocityTop(mm));

  end
end

fclose(fid);

