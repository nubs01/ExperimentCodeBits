function [meanPow,stdPow] = getBaselinePower(eegDat,eegTime,fs,freqRange,winSize,excludeTimes)
    % [meanPow,stdPow] = getBaselinePower(eegDat,eegTime,freqRange,winSize,excludeTimes)
    % TODO: do this using Hilbert amplitude as well 
    nFreqs = 128;
    freqs = linspace(freqRange(1),freqRange(2),nFreqs);
    [S,freqs,windowTimes] = getPSD(eegDat,fs,'freqRange',freqRange,'winSize',winSize,'freqRes',mean(diff(freqs)),'noverlap',winSize/2,'psdType','stft','specgram',1); 
    windowTimes = windowTimes' - winSize/2;
    windowTimes = [windowTimes windowTimes+winSize];
    windowTimes = windowTimes + eegTime(1);
    for k=1:size(excludeTimes,1)
        t1 = excludeTimes(k,1);
        t2 = excludeTimes(k,2);
        rmvIdx = ((windowTimes(:,1)>=t1 & windowTimes(:,1)<=t2) | (windowTimes(:,2)>=t1 & windowTimes(:,2)<=t2) | (windowTimes(:,1)<t1 & windowTimes(:,2)>t2));
        windowTimes = windowTimes(~rmvIdx,:);
        S = S(:,~rmvIdx);
    end
    if size(S,2)<3
        error('Not enough windows to get baseline power')
    end
    meanPow = mean(mean(S,2)); % mean power of whole freqRange
    stdPow = sqrt(mean(std(S).^2));  % propogate standard deviation
