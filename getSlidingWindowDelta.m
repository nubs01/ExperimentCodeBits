function [winDelta,winTime] = getSlidingWindowDelta(expDir,animID,day,tet,winSize,step)
% winSize and step in seconds
dataDir = sprintf('%s%s%s_direct%s',expDir,filesep,animID,filesep);
deltaFiles = dir(sprintf('%sEEG%s%sdelta%02i-*-%02i.mat',dataDir,filesep,animID,day,tet));
nEpochs = numel(deltaFiles);
deltaFiles = sort({deltaFiles.name});
for i=1:numel(deltaFiles)
    d = load([dataDir EEG filesep deltaFiles{i}]);
    d = d.delta;
    %if ~exist('delta','var')
    %    delta = d;
    %else
    %    delta{day}{i}{tet} = d{day}{i}{tet};
    %end
    deltaTime{i} = d{day}{i}{tet}.starttime:1/d{day}{i}{tet}.samprate:d{day}{i}{tet}.endtime;
    deltaPow{i} = d{day}{i}{tet}.deta(:,3);
    winN = fix(winSize*d{day}{i}{tet}.samprate);
    stepN = fix(step*d{day}{i}{tet}.samprate);
    avgDelta = slidingAvg([deltaTime{i}' deltaPow{i}],winN,stepN);
    winDelta{i} = avgDelta(:,2);
    winTime{i} = avgDelta(:,1);
end
        

