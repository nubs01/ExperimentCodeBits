%% Setup initial Data
animID = 'RW3';
dataDir = [pwd filesep];
fileparts(pwd)
expDir = [fileparts(pwd) filesep];
ee = 1;
tt = 2;
dd = 3;
eegDir = [dataDir 'EEG' filesep];
eegFN = @(s,d,e,t) sprintf('%s%s%02i-%02i-%02i.mat',[eegDir animID],s,d,e,t);
eeg = load(eegFN('eeg',dd,ee,tt));
eeg = eeg.eeg;
eeg = eeg{dd}{ee}{tt};
eegTime = eeg.starttime:1/eeg.samprate:eeg.endtime;
eegDat = eeg.data;
pos = load(sprintf('%spos%02i.mat',animID,dd));
pos = pos.pos{dd}{ee};
posTime = pos.data(:,1);
vel = pos.data(:,5);
[vel,velTime] = cleanVelocity(vel,posTime);
newVel = interp1(velTime,vel,eegTime);
freqBands = {[0.5 4.5],[0.5 9],[0.5 20],[0.5 55]};
artifactBand = [50 400];
[b,a] = butter(10,artifactBand/(eeg.samprate/2),'bandpass');
eegArt = filtfilt(b,a,eegDat);
h = viewTrace([eegDat eegArt newVel'],eegTime,10);
artifactTimes = detectArtifacts(eegTime,eegDat,eeg.samprate);
for l=1:numel(artifactTimes)
    a = artifactTimes(l)-.2;
    b = artifactTimes(l)+.2;
    idx = find(eegTime>a & eegTime<b);
    hold on
    plot(eegTime(idx),eegDat(idx),'r')
end

%% Test Clustering of state
[ratioVec,~,epochTimes] = scoreStateLFP(eegTime,eegDat,eeg.samprate);
epochVel = zeros(size(epochTimes,1),1);
for i=1:size(epochTimes,1)
    idx = (velTime>epochTimes(i,1) & velTime<epochTimes(i,2));
    epochVel(i) = mean(vel(idx));
end
clustDat = [cell2mat(ratioVec) epochVel];
figure()
plot(ratioVec{2},ratioVec{1},'k.')
[kidx,kC] = kmeans(clustDat,3);
colors = 'cmgkbry';
figure()
for i=1:numel(kidx)
    if ~isnan(kidx(i))
        hold on
        plot(ratioVec{2}(i),ratioVec{1}(i),[colors(kidx(i)) '.'])
    end
end
title('K-Means clusters')
Tidx = clusterdata(clustDat,'linkage','ward','maxclust',3);
colors = 'cmgkbry';
figure()
for i=1:numel(kidx)
    if ~isnan(kidx(i))
        hold on
        plot(ratioVec{2}(i),ratioVec{1}(i),[colors(Tidx(i)) '.'])
    end
end
title('Clusterdata Clustering')
