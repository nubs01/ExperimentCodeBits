% metadata
animIDs = {'RW2','RW3','RZ2','RZ3'};
genotypes = {'WT','WT','Df1','Df1'};
projDir = '~/Projects/rn_Schizophrenia_Project/';
expDirs = strcat(animIDs,'_');
expDirs = strcat(expDirs,genotypes);
expDirs = strcat(expDirs,'_Experiment/');
dataDirs = strcat(projDir,strcat(expDirs,strcat(animIDs,'_direct/')));
baseFigDir = [projDir 'RZ-RW_PrelimAnalysisFigs/Delta-Theta_Distributions/'];
mkdir(baseFigDir)
% testEpochs = {'sleep','saline','ketamine'};
testEpoch = 'sleep';
Nanim = numel(animIDs);
normBand = [0.5 100];
power_calc = 'stft';
power_method = 'total';
% I can either take mean_band_power/mean_spectrum_power or total_band_power/total_spectrum_power with hilbert power or stft power
%PCs = {'hilbert','stft'};
%PMs = {'mean','total'};
%NBs = {[0.5 55],[0.5 100]};
%NMVARS = {'hilbert','mean',[0.5 55];'hilbert','total',[0.5 55];'hilbert','mean',[0.5 100];'hilbert','total',[0.5 100];...
%          'stft','mean',[0.5 55];'stft','total',[0.5 55];'stft','mean',[0.5 100];'stft','total',[0.5 100];'stft','total',[100 250]};
%for NT = 9:size(NMVARS,1)
   % normBand = NMVARS{NT,3};
   % power_method = NMVARS{NT,2};
   % power_calc = NMVARS{NT,1};


    figDir = [baseFigDir 'Raw2' filesep];
    mkdir(figDir)

    analysisParams = struct('normalization','raw: 10*log10','power_calc',power_calc,'normBand',[],'power_method',power_method,'freqBand',{[1 5],[6 12],[1 12]},'freqBandName',{'delta','theta','delta&theta'},'epoch',testEpoch);
    behavioralStates = {'sleeping','resting','moving'};
    stateVelRanges = {[0 .3],[0 .3],[.4 80]}; % cm/sec - need to look at states found to check limits
    stateDurationRanges = {[40 1800],[10 40],[10 1800]}; % seconds
    allowedInterrupts = {0,0,2}; % seconds outside of velocity thresh allowed in an  episode (to account for breif stops while moving and lags in pos tracking)
    stateParams = struct('state',behavioralStates,'velRange',stateVelRanges,'durRange',stateDurationRanges,'allowedInterrupt',allowedInterrupts);

    minOfData = 28; % cutoff first and last 2 minutes of data
    velFs = 25;
    eegFs = 1500;
    winSize = 10;
    genoGroups = {'WT','Df1'};
    eegTraces = cell(Nanim,1);
    % eegrefTraces = cell(Nanim,1);
    velTraces = cell(Nanim,1);
    datIdx = cell(Nanim,1);
    tenSecEEG = cell(Nanim,1);
    tenSecIdx = cell(Nanim,1);
    tenSecNorm = cell(Nanim,1);

    for k=1:Nanim
        tetinfo = load(sprintf('%s%stetinfo.mat',dataDirs{k},animIDs{k}));
        rntask  = load(sprintf('%s%srntask.mat',dataDirs{k},animIDs{k}));
        tetinfo = tetinfo.tetinfo;
        rntask  = rntask.rntask;
        tetIdx = getCellStructIdx(tetinfo,'descrip','riptet','area','CA1');
        taskIdx = getCellStructIdx(rntask,'epoch',testEpoch);
        days = taskIdx(:,1);
        epochs = taskIdx(:,2);
        tmp1 = arrayfun(@(x,y) any(days==x) & any(epochs==y),tetIdx(:,1),tetIdx(:,2));
        tetIdx = tetIdx(tmp1,:);
        
        eegTraces{k} = getAnimEEG(dataDirs{k},animIDs{k},tetIdx,'truncateDat',minOfData);
    %    eegrefTraces{k} = getAnimEEG(dataDirs{k},animIDs{k},tetIdx,'truncateDat',minOfData);
        datIdx{k} = [repmat(k,size(tetIdx,1),1) tetIdx];
        velTraces{k} = getAnimVel(dataDirs{k},animIDs{k},tetIdx(:,1:2),'velFs',velFs,'truncateDat',minOfData);
        velTime = 0:1/velFs:(size(velTraces{k},2)-1)/velFs;
        eegTime = 0:1/eegFs:(size(eegTraces{k},2)-1)/eegFs;

        tmptenSec = cell(size(eegTraces{k},1),1);
        tmpStates = tmptenSec;
        tmpNormPow = tmptenSec;

        for l = 1:size(eegTraces{k},1)
            tmpStateMat = getBehavioralEpisodes(velTime,velTraces{k}(l,:),stateParams);
            [tmptenSec{l},tmpIdx] = slidingWindow(eegTraces{k}(l,:),winSize*eegFs,winSize*eegFs);
            tmptenSec{l} = tmptenSec{l}';
            tmpTime = eegTime(tmpIdx);
            tmpStates{l} = zeros(size(tmpIdx,2),2);
            tmpNormPow{l} = zeros(size(tmpIdx,2),2);
            tmpStates{l}(:,1) = find(strcmp(genoGroups,genotypes{k}));
            [tmpNormPow{l}(:,1),tmpNormPow{l}(:,2)] = getBandPower(eegTraces{k}(l,:),eegFs,normBand,'powerType',power_calc,'method',power_method);
            for m=1:size(tmpIdx,2)
                a = find(tmpStateMat(:,1)<=tmpTime(1,m) & tmpStateMat(:,2)>=tmpTime(2,m));
                if isempty(a)
                    tmpStates{l}(m,2) = -1;
                else
                    tmpStates{l}(m,2) = tmpStateMat(a,3);
                end
            end
        end
        tenSecEEG{k} = cell2mat(tmptenSec);
        tenSecNorm{k} = cell2mat(tmpNormPow);
        tenSecIdx{k} = cell2mat(tmpStates);
        clear tmp*
    end


    eegDat = cell2mat(eegTraces);
    velDat = cell2mat(velTraces);
    % eegrefDat = cell2mat(eegrefTraces);
    eegFs = 1500; % TODO: fix hardcoding
    velFs = 25; % TODO: fix hardcoding
    dataIdx = cell2mat(datIdx);
    velTime = 0:1/velFs:(size(velDat,2)-1)/velFs;
    eegTime = 0:1/eegFs:(size(eegDat,2)-1)/eegFs;
    tenSecDat = cell2mat(tenSecEEG);
    tenSecIdx = cell2mat(tenSecIdx);
    tenSecNormPow = cell2mat(tenSecNorm);
    winSize = winSize*eegFs;

    for O = 1:numel(analysisParams)
        fBand = analysisParams(O).freqBand;
        fBName = analysisParams(O).freqBandName;
        normBand = analysisParams(O).normBand;
        power_calc = analysisParams(O).power_calc;
        epoch = analysisParams(O).epoch;
        datasetName = [fBName '_' epoch '_dataset'];
        dataParams = analysisParams(O);
        power_method = analysisParams(O).power_method;

        powDat = cell(1,numel(genoGroups));
        powErr = powDat;
        fprintf('Getting %s %s power using %s...\n',power_method,fBName,power_calc);
        for gg = 1:numel(genoGroups)
            %animsInGeno = find(strcmp(genotypes,genoGroups{gg}));
            %traceIdx = find(arrayfun(@(x) any(animsInGeno==x),dataIdx(:,1)));
            powDat{gg} = cell(1,numel(stateParams));
            powErr{gg} = powDat{gg};
            for ss = 1:numel(stateParams)
                stateIdx = (tenSecIdx(:,1)==gg & tenSecIdx(:,2)==ss);
                tmpTraces = tenSecDat(stateIdx,:);
                [tmpPow,tmpErr] = getBandPower(tmpTraces,eegFs,fBand,'powerType',power_calc,'method',power_method);
                powDat{gg}{ss} = tmpPow;%./tenSecNormPow(stateIdx,1);
                powErr{gg}{ss} = tmpErr;%abs(powDat{gg}{ss}).*sqrt((tmpErr./tmpPow).^2 + (tenSecNormPow(stateIdx,2)./tenSecNormPow(stateIdx,1)).^2);
            end
        end
        save([figDir datasetName '.mat'],'powDat','powErr','dataParams');

        % make plots
        distribFig = figure();
        histColors = 'bmcrgky';
        setFigureProperties(distribFig,'Position',[1 1 1440 600])
        for ss=1:numel(stateParams)
            subplot(1,numel(stateParams),ss)
            hold on
            binWidth = [];
            for gg=1:numel(genoGroups)
                if isempty(binWidth)
                    h = histogram(10*log10(powDat{gg}{ss}),'FaceColor',histColors(gg));
                    binWidth = h.BinWidth;
                else
                    h = histogram(10*log10(powDat{gg}{ss}),'FaceColor',histColors(gg),'BinWidth',binWidth);
                end
            end
            if ss==1
                ylabel('Counts')
            end
            if ss==ceil(numel(stateParams)/2)
                xlabel('Raw Power')
                legend(genoGroups)
            end
            title(stateParams(ss).state)
        end
        suptitle({sprintf('Raw %s band power: %g - %g Hz',fBName,fBand),sprintf('%s %s power: log scale',capFirst(power_method),power_calc)})
        pause(2)
        saveas(distribFig,[figDir fBName '_Distributions'],'svg')
    end    
    close all
%end
