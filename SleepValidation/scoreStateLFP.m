function [ratioVectors,states,epochTimes] = scoreStateLFP(lfpTime,lfpDat,fs,varargin)
    % [states,epochTimes] = scoreStateLFP(lfpTime,lfpDat,NAME,VALUE)
    % break lfp in to 10 second epochs and scores each as
    % 'sleep','wake','artifact' based on the spectral power ratios used in
    % Gervasoni et al 2004. Ratio 1: 0.5-20Hz / 0.5-55Hz, Ratio 2: 0.5-4.5Hz /
    % 0.5-9Hz. 
    % * Could possibly separate Wake, SWS & REM 

    % Setup initial Variables
    freqBands = {[0.5 20],[0.5 55],[0.5 4.5],[0.5 9]};
    compRatios = {[1 2],[3 4]};
    epochSize = 10;
    winSize = 2;
    winStep = 1;
    artifactTimes = detectArtifacts(lfpTime,lfpDat,fs,'stdThresh',5);
    params = struct('tapers',[3 4],'Fs',fs,'trialave',1);

    % Assign Name, Values pairs from input to overwrite any initial vars
    if ~isempty(varargin)
        assignVars(varargin{:});
    end

    % Setup other variables
    winSize = fix(winSize*fs); % convert from sec to points
    winStep = fix(winStep*fs);
    epochSize = fix(epochSize*fs);
    N = numel(lfpTime);
    epochIdxs = (1:epochSize:N-epochSize+1)';
    epochIdxs = [epochIdxs epochIdxs+epochSize-1];
    Nepochs = size(epochIdxs,1);
    minFreq = min(cellfun(@min,freqBands));
    maxFreq = max(cellfun(@max,freqBands));
    lfpDat = interpInf(lfpDat);
    params.fpass = [minFreq maxFreq];

    % Define Outputs
    epochTimes = lfpTime(epochIdxs);
    states = cell(Nepochs,1);

    % Collect band power for each epoch
    bandPower = zeros(Nepochs,numel(freqBands));
    badEpochs = [];
    for l=1:Nepochs
        i1 = epochIdxs(l,1);
        i2 = epochIdxs(l,2);
        tmpDat = lfpDat(i1:i2);
        tmpTime = lfpTime(i1:i2);
        winIdx = 1:winStep:numel(tmpDat)-winSize+1;
        winIdx = arrayfun(@(x,y) (x:y)',winIdx,winIdx+winSize-1,'UniformOutput',false);
        winIdx = cell2mat(winIdx);
        winDat = tmpDat(winIdx);
        winTime = tmpTime(winIdx);

        % remove any windows containing artifacts
        rmv = arrayfun(@(x,y) any(artifactTimes>=x & artifactTimes<=y),winTime(1,:),winTime(end,:));
        winDat(:,rmv) = [];

        if size(winDat,2)<3
            bandPower(l,:) = nan(1,numel(freqBands));
            badEpochs = [badEpochs l];
            continue;
        end

        [S,f] = mtspectrumc(winDat,params);

        for k=1:numel(freqBands)
            fidx = (f>=freqBands{k}(1) & f<=freqBands{k}(2));
            % bandPower(l,k) = mean(S(fidx));
            bandPower(l,k) = sum(S(fidx));
        end
        clear S f tmpDat tmpTime winDat winTime winIdx i1 i2 rmv fidx k
    end

    ratioVectors = cell(1,numel(compRatios));
    for l=1:numel(compRatios)
        a = compRatios{l}(1);
        b = compRatios{l}(2);
        ratioVectors{l} = bandPower(:,a)./bandPower(:,b);
    end





