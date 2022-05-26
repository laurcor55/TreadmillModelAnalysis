clear all; close all;

v = VideoWriter('PFMovie.mp4');

addpath('../model/lib')
addpath('../model/lib/analysis', '../model/lib/data', '../model/lib/dependencies', '../model/lib/export', '../model/lib/gillespie');


%% Full model can run from this script. User can alter Parameters and Kinetics here.
% Parameters Input:
Parameters = parameters();
  Parameters.totalTime = 10;
  Parameters.concTotalFtsZ = 3;

% Kinetics Input:
Kinetics = kineticsBs();

Outputs = runExperiment(Parameters, Kinetics, false);
matrixPfSnapshots = cellsTo3DMatrix(Outputs.savePfs, Outputs.saveLocations);
[x, y, z] = size(matrixPfSnapshots);
frameLen = mean(Outputs.time(2:end) - Outputs.time(1:end-1));
v.FrameRate = 1/frameLen;


fig = figure
set(fig,'units','pixels','position',[0 0 1280 960])
v = VideoWriter('pfMovie.avi');
tstep = 1./30;
t = tstep:tstep:Parameters.totalTime;
set(gca,'nextplot','replacechildren');
Colormap = colormap(protocolormap);

open(v);
matrix = zeros(x,y,length(t));
for ii = 1:length(t)
  [~, ind] = min(abs(t(ii)-Outputs.time));
  matrix(:,:,ii) = flipud(matrixPfSnapshots(:,:,ind));
  image(matrix(:,:,ii)+1);
%  xlim([0,300])
%  ylim([0,300])
  axis off
  pause(frameLen);
  frame = getframe;
  writeVideo(v,frame);
end

close(v);

t = [2, 5, 10, 30];

matrix = zeros(x,y,length(t));
for ii = 1:length(t)
  [~, ind] = min(abs(t(ii)-Outputs.time));
  matrix(:,:,ii) = matrixPfSnapshots(:,:,ind);
  subplot(2, 2, ii)
  image(matrix(:,:,ii)+1);
  axis off
end
