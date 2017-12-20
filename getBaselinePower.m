function [meanPow,stdPow] = getBaselinePower(eegDat,eegTime,fs,freqRange,winSize,excludeTimes)
    % [meanPow,stdPow] = getBaselinePower(eegDat,eegTime,freqRange,winSize,excludeTimes)

    nFreqs = 128;
    freqs = linspace(freqRange(1),freqRange(2),nFreqs);
    params = struct('fpass',freqRange,'Fs',fs,'tapers',[3 4],'err',[2 0.05]);
    [S,f,~,~,Serr] = mtspectrumsegc(eegDat,winSize,params,0);
    windowTimes = (0:winSize:(numel(eegDat)/fs)-winSize)'+[0 winSize];
    windowTimes = windowTimes + eegTime(1);
    for k=1:size(excludeTimes,1)
        t1 = excludeTimes(k,1);
        t2 = excludeTimes(k,2);
        rmvIdx = ((windowTimes(:,1)>=t1 & windowTimes(:,1)<=t2) | (windowTimes(:,2)>=t1 & windowTimes(:,2)<=t2));
        windowTimes = windowTimes(~rmvIdx,:);
        S = S(:,~rmvIdx);
        Serr = Serr(:,:,~rmvIdx);
    end
    S1 = interp1(f,S,freqs,'linear','extrap');
    meanPow = mean(mean(S1)); % mean power of whole freqRange
    stdPow = sqrt(mean(std(S1).^2));  % propogate standard deviation
