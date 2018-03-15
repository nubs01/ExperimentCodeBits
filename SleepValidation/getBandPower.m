function [bandPow,bandErr] = getBandPower(eegTraces,fs,freqBand,varargin)
    % [bandPow,bandErr] = getBandPower(eegTraces,fs,freqBand,NAME,VALUE)
    % eegTraces should be a vector or a matrix with an eegTrace as each row.
    % Retruns the spectral power in the band freqBand.  
    % NAME-VALUE pairs:
    %   - powerType    : string designating a method for obtaining spectral power. options are 'hilbert' (default) or 'stft'
    %   - winSize      : band power will be computed in sliding windows of winSize (seconds) and averaged. [default=2]
    %   - winStep      : size of steps to take with the sliding window. [default=winSize/2]
    %   - eegTime      : vector or matrix matching eegTraces, with the time at
    %                    each point (in seconds). only used if time periods should be excluded
    %                    from the power calculation. 
    %   - excludeTimes : a 2-column matrix with the start and end times of periods to exclude when computing band power. 
    %   - excludeFunc  : a function handle that can be passed an eegTrace and
    %                    return a 2-column matrix with start and end indices of periods to
    %                    exclude from power calculation. such as an artifact detection function.
    %                    [default=[]]
    %   - method       : mean or total power in freqBand. default='mean'
    %   - nboot        : number of times to bootstrap for SEM and CI

    powerType = 'hilbert';
    winSize = 2;
    winStep = [];
    eegTime = 0:1/fs:(size(eegTraces,2)-1)/fs;
    excludeTimes = [];
    excludeFunc = [];
    method = 'mean';
    nboot = 10000;

    assignVars(varargin{:})

    if isempty(winStep)
        winStep = winSize/2;
    end

    isContained = @(x,y,w,v) (w<=x & v>=x) | (w<=y & v>=y) | (w>=x & v<=y);

    [b,a] = butter(3,freqBand/(fs/2),'bandpass');
    Ntraces = size(eegTraces,1);
    winSize = fix(winSize*fs); % window size in pts
    winStep = fix(winStep*fs); % window step in pts
    bandPow = zeros(Ntraces,1);
    bandErr = zeros(Ntraces,1);

    for k=1:Ntraces
        tmp1 = filtfilt(b,a,eegTraces(k,:));
        if any(isnan(tmp1))
            error('Filtering error...nans returned')
        end
        
        [tmp2,tmpIdx] = slidingWindow(tmp1,winSize,winStep);
        if ~isempty(excludeFunc)
            eT = excludeFunc(tmp1);
            rmv = isContained(eT(:,1),eT(:,2),tmpIdx(1,:),tmpIdx(2,:));
            rmv = any(rmv);
            tmp2(:,rmv) = [];
            tmpIdx(:,rmv) = [];
        end
        if ~isempty(excludeTimes)
            tmpTimeIdx = eegTime(tmpIdx);
            rmv = isContained(excludeTimes(:,1),excludeTimes(:,2),tmpTimeIdx(1,:),tmpTimeIdx(2,:));
            rmv = any(rmv);
            tmp2(:,rmv) = [];
            tmpIdx(:,rmv) = [];
        end
        switch lower(powerType)
            case 'hilbert'
                [tmpPow,tmpErr] = hilbertPower(tmp2,method);
            case 'stft'
                params = struct('tapers',[3 4],'Fs',fs,'fpass',freqBand);
                S = mtspectrumc(tmp2,params);
                if strcmp(method,'mean')
                    tmpPow = mean(S);
                    tmpErr = std(S);
                elseif strcmp(method,'total')
                    tmpPow = sum(S);
                    tmpErr = [];
                else
                    error('Invalid method');
                end
            otherwise
                error('Invalid choice for power calculation')
        end
        bandPow(k) = mean(tmpPow);
        if numel(tmpPow)>10
            bootstat = bootstrp(nboot,@mean,tmpPow);
            bandErr(k) = std(bootstat);
        elseif numel(tmpPow)>=3
            bandErr(k) = std(tmpPow);
        else
            bandErr(k) = nan;
        end
    end
 


function [P,sigma] = hilbertPower(eegWin,method)
    a = envelope(eegWin);
    if strcmp(method,'mean')
        P = mean(a.^2);
        sigma = std(a.^2);
    elseif strcmp(method,'total')
        P = sum(a.^2);
        sigma = [];
    else
        error('Invalid Method')
    end
