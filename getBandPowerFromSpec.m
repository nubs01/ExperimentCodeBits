function [powMean,powSTD,baselineMean,baselineSTD] = getBandPowerFromSpec(S,t,f,freqBand,normFlag,baselineTime)
% If normFlag = 1, the returned power trace will be normalized to the mean and std during baselineTime
    if ~exist('normFlag','var')
        normFlag = 0;
    end
    if ~exist('baselineTime','var')
        baselineTime = [t(1) t(end)];
    end
    freqIdx = find(f>=freqBand(1) & f<=freqBand(2));
    baselineIdx = find(t>=baselineTime(1) & t<=baselineTime(2));
    
    % make sure S is orientated so that each row is a time point and each column is a frequency
    if size(S,2) == numel(t)
        S = S';
    end

    if size(S,1) ~= numel(t) || size(S,2) ~= numel(f)
        error('S must have dimensions matching t and f')
    end

    powBand = S(:,freqIdx);
    powMean = mean(powBand,2);
    powSTD = std(powBand,0,2);
    
    baselineMean = mean(powBand(baselineIdx,:),2);
    baselineSTD = std(powBand(baselineIdx,:),0,2);

    % get pooled standard deviation and mean for baseline
    ncol = size(powBand,2);
    nrow = size(powBand,2);
    baselineMean = mean(baselineMean); % since all rows are same length
    baselineSTD = sqrt(mean(baselineSTD.^2));
    
    if normFlag % now z-score powMean and propogate error into powSTD
        powMean = (powMean-baselineMean)/baselineSTD;
        powSTD = sqrt(powSTD.^2+baselineSTD^2)./baselineSTD;
    end


