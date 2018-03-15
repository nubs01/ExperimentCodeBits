function [outDataStruct,dataParams] = collectPrelimDataset(varargin)
% collects data from RW2 RW3, RZ2 and RZ3 
% breaks CA1 riptet LFP and Velocity data in to winSize second windows and for
% each 10 sec window creates a structure entry with fields: animal, genotype,
% day, epoch, tetrode, behavioral_state ('sleeping','resting','moving','transition'),
% band_power (array with col for each band in brandFreqs), band_freqs (cell
% array with 2-element frequency bands [band_min band_max]), band_names,
% normalization_power, peak_frequency, preak_freq_range
% Also returns a structure dataParams with fields: delta_freq, theta_freq,
% normalization_freq, power_calc, power_method, min_of_data, behavioralParams,
% eegFs, velFs, window_size, power_window ([window step] : to average power
% spectra over for each window)


% animal metadata
animIDs = {'RW2','RW3','RZ2','RZ3'};
genotypes = {'WT','WT','Df1','Df1'};
projDir = '~/Projects/rn_Schizophrenia_Project/';
expDirs = strcat(projDir,animIDs,'_',genotypes,'_Experiment',filesep);
dataDirs = strcat(expDirs,animIDs,'_direct',filesep);
epochs = {};


% analysis metadata
bandFreqs = {[2 5],[5 10],[1 12]}; % Hz
bandNames = {'delta','theta','delta-theta'};
power_calc = 'stft';
power_method = 'total';
normBand = [0.5 100]; % Hz
eegFs = 1500; % Hz
velFs = 25; % Hz
winSize = 10; % seconds
powerWin = [2 1]; % seconds
minOfData = 28; % minutes
peakFreqRange = [0.5 12]; % Hz
nboot = 1000; % for SEM calculations
eegRefCutoff = 50; % Hz; any freq band greater than this will use eegRef instead off eeg for power analysis 

% behavioral state parameters
behavioralStates = {'sleeping','resting','moving'};
stateVelRanges = {[0 .3],[0 .3],[.4 80]};
stateDurationRanges = {[40 1800],[10 40],[10 1800]};
allowedInterrupts = {0,0,2};
stateParams = struct('state',behavioralStates,'velRange',stateVelRanges,'durRange',stateDurationRanges,'allowedInterrupt',allowedInterrupts);

% Spectrum parameters
specFreqRange = [0.5 55];
specFreqRes = 0.25;

% override vars with input
assignVars(varargin)

% setup final analysis variables from previous varaibles
Nanim = numel(animIDs);
Nstates = numel(stateParams);
Nepochs = numel(epochs);
Nbands = numel(bandFreqs);

specFreq = specFreqRange(1):specFreqRes:specFreqRange(2);

% create empty output variables
dataStruct = cell(1,Nanim);
dataParams = struct('band_freqs',{bandFreqs},'band_names',{bandNames},...
                    'normalization_freq',normBand,'power_calc',power_calc,...
                    'power_method',power_method,'eegFs',eegFs,...
                    'velFs',velFs,'window_size',winSize,...
                    'power_window',powerWin,'minutes_per_epoch',minOfData,...
                    'behavioral_params',stateParams,'peak_freq_range',peakFreqRange,...
                    'specFreqRange',specFreqRange,'specFreqVec',specFreq);

% loop through animals

for k=1:Nanim
    fprintf('Collecting data for %s ...\n',animIDs{k});
    % load task and tetinfo
    tetinfo = load(sprintf('%s%stetinfo.mat',dataDirs{k},animIDs{k}));
    tetinfo = tetinfo.tetinfo;
    rntask = load(sprintf('%s%srntask.mat',dataDirs{k},animIDs{k}));
    rntask = rntask.rntask;

    % get valid indices for data
    tetIdx = getCellStructIdx(tetinfo,'descrip','riptet','area','CA1');
    if Nepochs>0
        taskIdx = [];
        for e=1:Nepochs
            tmp = getCellStructIdx(rntask,'epoch',epochs{e});
            taskIdx = [taskIdx;tmp];
        end
    else
        taskIdx = getCellStructIdx(rntask,'animal',animIDs{k});
    end
    tmp = arrayfun(@(x,y) any(taskIdx(:,1)==x) & any(taskIdx(:,2)==y),tetIdx(:,1),tetIdx(:,2));
    tetIdx = tetIdx(tmp,:);
    clear taskIdx tmp 

    % empty variables to fill
    Ntraces = size(tetIdx,1);
    tmpDataStruct = cell(1,Ntraces);
    fprintf('   Found %g LFP traces to process\n    Collecting Metrics...\n',Ntraces)

    % loop through all [day epoch tet] and compute variables 
    for l=1:Ntraces

        fprintf('       Processing trace %g\n',l)
        % get LFP and Velocity
        eegTrace = getAnimEEG(dataDirs{k},animIDs{k},tetIdx(l,:),'truncateDat',minOfData);
        if iscell(eegTrace)
            eegTrace = eegTrace{1};
        end
        if any(strcmp(bandNames,'gamma')) || any(strcmp(bandNames,'ripple'))
            eegRef = getAnimEEG(dataDirs{k},animIDs{k},tetIdx(l,:),'truncateDat',minOfData,'dataType','eegref');
            if iscell(eegRef)
                eegRef = eegRef{1};
            end
        else
            eegRef = [];
        end
        velTrace = getAnimVel(dataDirs{k},animIDs{k},tetIdx(l,1:2),'velFs',velFs,'truncateDat',minOfData);
        if iscell(velTrace)
            velTrace = velTrace{1};
        end
        pos = load(sprintf('%s%spos%02i.mat',dataDirs{k},animIDs{k},tetIdx(l,1)));
        pos = pos.pos{tetIdx(l,1)}{tetIdx(l,2)};
        xTrace = pos.data(:,2);
        yTrace = pos.data(:,3);
        posTime = pos.data(:,1);
        posTime = posTime-posTime(1);
        clear pos
        velTime = 0:1/velFs:(size(velTrace,2)-1)/velFs;
        eegTime = 0:1/eegFs:(size(eegTrace,2)-1)/eegFs;

        % get behavioral states
        tmpStateMat = getBehavioralEpisodes(velTime,velTrace,stateParams);

        % break into windows
        [winEEG,winIdx] = slidingWindow(eegTrace,winSize*eegFs,winSize*eegFs);
        [winVel,winVIdx] = slidingWindow(velTrace,winSize*velFs,winSize*velFs);
        winEEG = winEEG';
        winVel = winVel';

        Nwin = size(winEEG,1);
        tmpPow = zeros(Nwin,Nbands);
        tmpPeak = zeros(1,Nwin);
        tmpState = zeros(1,Nwin);
        tmpVelSeg = zeros(1,Nwin);
        tmpVelSegErr = zeros(1,Nwin);
        tmpDistSeg = zeros(1,Nwin);
        tmpWinTime = cell(1,Nwin);
        tmpMaxVel = zeros(1,Nwin);
        tmpSpectrum = zeros(Nwin,numel(specFreq));

        % get normalization power
        [tmpNorm,tmpNormErr] = getBandPower(eegTrace,eegFs,normBand,'winSize',powerWin(1),...
                                            'winStep',powerWin(2),'powerType',power_calc,...
                                            'method',power_method);
        % get epoch velocity & distance
        tmpVelEpo = mean(velTrace);
        tmpVelEpoErr = std(bootstrp(nboot,@mean,velTrace));
        tmpDistEpo = getTotalDistance(xTrace,yTrace);

        % loop of winSize segments
        for m=1:Nwin

            %get distance and velocity metrics
            tv = velTime(winVIdx(:,m));
            te = eegTime(winIdx(:,m));
            tmpWinTime{m} = te';
            idx = (posTime>=tv(1) & posTime<=tv(2));
            tmpDistSeg(m) = getTotalDistance(xTrace(idx),yTrace(idx));
            tmpVelSeg(m) = mean(winVel(m,:));
            tmpVelSegErr(m) = std(bootstrp(nboot,@mean,winVel(m,:)));
            tmpMaxVel(m) = max(winVel(m,:));
            % get state of segment
            idx = (tmpStateMat(:,1)<=te(1) & tmpStateMat(:,2)>=te(2));
            if ~any(idx)
                tmpState(m) = -1;
            else
                tmpState(m) = tmpStateMat(idx,3);
            end

            % get power for each freq band
            for n=1:Nbands
                if all(bandFreqs{n}>=eegRefCutoff)
                    tmpBandTrace = eegRef(winIdx(1,m):winIdx(2,m));
                else
                    tmpBandTrace = winEEG(m,:);
                end
                tmpPow(m,n) = getBandPower(tmpBandTrace,eegFs,bandFreqs{n},...
                                           'winSize',powerWin(1),'winStep',powerWin(2),...
                                           'powerType',power_calc,'method',power_method);

                % normalize high-freq band power to eegRef trace
                if all(bandFreqs{n}>eegRefCutoff)
                    np = getBandPower(eegRef,eegFs,normBand,'winSize',powerWin(1),...
                                      'winStep',powerWin(2),'powerType',power_calc,...
                                      'method',power_method);
                    tmpPow(m,n) = tmpPow(m,n)/np;
                end
            end
            clear tmpBandTrace

            % get peak frequency
            [tmpSpectrum(m,:),tmpF] = getPSD(winEEG(m,:),eegFs,'psdType',power_calc,...
                             'winSize',powerWin(1),'noverlap',powerWin(1)-powerWin(2),...
                             'freqRange',specFreqRange,'freqRes',specFreqRes);
            [~,idx] = max(tmpSpectrum(m,:));
            tmpPeak(m) = tmpF(idx); % TODO: Fit curve and find peak

        end

        % create struct array for this [anim day epoch tet]

        tpow = num2cell(tmpPow,2)';
        tpeak = num2cell(tmpPeak);
        tstate = cell(1,Ntraces);
        tstate(tmpState>0) = behavioralStates(tmpState(tmpState>0));
        tstate(tmpState==-1) = {'transition'};
        tvelseg = num2cell([tmpVelSeg' tmpVelSegErr'],2)';
        tdistseg = num2cell(tmpDistSeg);
        tmaxv = num2cell(tmpMaxVel);
        tspec = num2cell(tmpSpectrum,2)';
        tmpStruct = struct('animal',animIDs{k},'genotype',genotypes{k},...
                           'day',tetIdx(l,1),'epoch',rntask{tetIdx(l,1)}{tetIdx(l,2)}.epoch,...
                           'tet',tetIdx(l,3),'epochIdx',tetIdx(l,2),...
                           'band_names',{bandNames},'band_freqs',{bandFreqs},...
                           'band_power',tpow,'normalization_power',[tmpNorm tmpNormErr],...
                           'epoch_mean_vel',[tmpVelEpo tmpVelEpoErr],'epoch_total_distance',tmpDistEpo,...
                           'peak_frequency',tpeak,'behavioral_state',tstate,...
                           'segment_mean_vel',tvelseg,'segment_total_distance',tdistseg,...
                           'peak_frequency_range',peakFreqRange,'segment_time',tmpWinTime,...
                           'segment_max_vel',tmaxv,'spectrum',tspec);
        tmpDataStruct{l} = tmpStruct;
    end
    dataStruct{k} = cell2mat(tmpDataStruct);
end

outDataStruct = cell2mat(dataStruct);

end



