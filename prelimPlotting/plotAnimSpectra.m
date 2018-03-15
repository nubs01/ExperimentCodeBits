% metadata
animIDs = {'RW2','RW3','RZ2','RZ3'};
genotypes = {'WT','WT','Df1','Df1'};
projDir = '~/Projects/rn_Schizophrenia_Project/';
expDirs = strcat(animIDs,'_');
expDirs = strcat(expDirs,genotypes);
expDirs = strcat(expDirs,'_Experiment/');
dataDirs = strcat(projDir,strcat(expDirs,strcat(animIDs,'_direct/')));
baseFigDir = [projDir 'RZ-RW_PrelimAnalysisFigs/NormSpectra/'];
mkdir(baseFigDir)
testEpochs = {'sleep','saline','ketamine'};
%testEpoch = 'sleep';
Nanim = numel(animIDs);
freqRange = [0.5 55];
freqRes = 0.25;
normBand = [0.5 100];
power_calc = 'stft';
power_method = 'total';
behavioralStates = {'sleeping','resting','moving'};
stateVelRanges = {[0 .3],[0 .3],[.4 80]}; % cm/sec - need to look at states found to check limits
stateDurationRanges = {[40 1800],[10 40],[10 1800]}; % seconds
allowedInterrupts = {0,0,2}; % seconds outside of velocity thresh allowed in an  episode (to account for breif stops while moving and lags in pos tracking)
stateParams = struct('state',behavioralStates,'velRange',stateVelRanges,'durRange',stateDurationRanges,'allowedInterrupt',allowedInterrupts);

for M=1:numel(testEpochs)
    testEpoch = testEpochs{M};

    analysisParams = struct('normalization','mean stft over 0.5-100Hz','power_calc',power_calc,'normBand',normBand,'power_method',power_method,'freqRange',freqRange,'freqRes',freqRes,'epoch',testEpoch);

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

        % TODO: FIX: Hard-coded to ignore day 1 epoch 2 for RW3 since it is shorter than the rest for some reason
        if strcmp(animIDs{k},'RW3') && strcmp(testEpoch,'saline')
            tetIdx = tetIdx(tetIdx(:,1)~=1,:);
        end
        
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
            tmpVel = velTraces{k}(l,:);
            if tmpVel(end)==-1
                tmpVel = tmpVel(tmpVel~=-1);
            end
            tmpStateMat = getBehavioralEpisodes(velTime,tmpVel,stateParams);
            [tmptenSec{l},tmpIdx] = slidingWindow(eegTraces{k}(l,:),winSize*eegFs,winSize*eegFs);
            tmptenSec{l} = tmptenSec{l}';
            tmpTime = eegTime(tmpIdx);
            tmpStates{l} = zeros(size(tmpIdx,2),3);
            tmpNormPow{l} = zeros(size(tmpIdx,2),2);
            tmpStates{l}(:,1) = k;
            tmpStates{l}(:,2) = find(strcmp(genoGroups,genotypes{k}));
            [tmpNormPow{l}(:,1),tmpNormPow{l}(:,2)] = getBandPower(eegTraces{k}(l,:),eegFs,normBand,'powerType',power_calc,'method',power_method);
            for m=1:size(tmpIdx,2)
                a = find(tmpStateMat(:,1)<=tmpTime(1,m) & tmpStateMat(:,2)>=tmpTime(2,m));
                if isempty(a)
                    tmpStates{l}(m,3) = -1;
                else
                    tmpStates{l}(m,3) = tmpStateMat(a,3);
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
    tenSecIdx = cell2mat(tenSecIdx); % [anim,geno,state]
    tenSecNormPow = cell2mat(tenSecNorm); % [power,error]

    normBand = analysisParams.normBand;
    power_calc = analysisParams.power_calc;
    epoch = analysisParams.epoch;
    datasetName = ['NormalizedSpectra_' epoch '_dataset'];
    power_method = analysisParams.power_method;
    freqRange = analysisParams.freqRange;
    freqRes = analysisParams.freqRes;

    % get normalized spectrum for each 10sec window
    freqVec = freqRange(1):freqRes:freqRange(2);
    tenSecSpec = zeros(size(tenSecDat,1),numel(freqVec));
    for k=1:size(tenSecDat,1)
        tenSecSpec(k,:) = getPSD(tenSecDat(k,:),eegFs,'psdType',power_calc,'winSize',2,'noverlap',1,'freqRange',freqRange,'freqRes',freqRes);
    end
    tenSecNormSpec = tenSecSpec./tenSecNormPow(:,1);
    [tenSecPeak,tmpIdx] = max(tenSecNormSpec,[],2);
    tenSecPeakFreq = freqVec(tmpIdx);
    if isrow(tenSecPeakFreq)
        tenSecPeakFreq = tenSecPeakFreq';
    end
    
    save([baseFigDir datasetName '.mat'],'tenSecNormSpec','tenSecSpec','tenSecIdx','tenSecDat','analysisParams','tenSecPeak','tenSecPeakFreq')

    winSize = winSize*eegFs;
    stateGroups = behavioralStates;
    animGroups = animIDs;
    gSpecDir = [baseFigDir 'GenoSpectra' filesep];
    mkdir(gSpecDir)
    aSpecDir = [baseFigDir 'AnimSpectra' filesep];
    mkdir(aSpecDir)
    pSpecDir = [baseFigDir 'PeakDistributions' filesep];
    mkdir(pSpecDir)
    Nstates = numel(stateGroups);

    % Make Plot 1: spectra for each geno in each state
    genoSpecFig = figure();
    plotColors = 'bmcgky';
    setFigureProperties(genoSpecFig,'Position',[1 1 1440 600])
    for k=1:Nstates
        subplot(2,Nstates,k)
        hold on
        genoH = gobjects(2,1);
        binWidth = [];
        for l=1:numel(genoGroups)
            subplot(2,Nstates,k)
            hold on
            idx = (tenSecIdx(:,2)==l & tenSecIdx(:,3)==k);
            [~,genoH(l)] = shadedErrorPlot(freqVec,mean(tenSecNormSpec(idx,:)),std(tenSecNormSpec(idx,:)),[plotColors(l) '-']);
            yl = get(gca,'ylim');
            yl(1) = 0;
            set(gca,'ylim',yl)
            set(gca,'xlim',freqRange)

            subplot(2,Nstates,k+Nstates)
            hold on
            if isempty(binWidth)
                h = histogram(tenSecPeakFreq(idx),'FaceColor',plotColors(l));
                binWidth = h.BinWidth;
            else
                h = histogram(tenSecPeakFreq(idx),'FaceColor',plotColors(l),'BinWidth',binWidth);
            end
        end
        if k==ceil(Nstates/2)
            subplot(2,Nstates,k)
            xlabel('Frequency (Hz)')
            legend(genoH,genoGroups)
            subplot(2,Nstates,Nstates+k)
            xlabel('Peak Frequency (Hz)')
        end
        if k==1
            subplot(2,Nstates,k)
            ylabel('Normalized Power: log scaled')
            subplot(2,Nstates,Nstates+k)
            ylabel('Counts')
        end
        subplot(2,Nstates,k)
        title(stateGroups{k})
    end
    suptitle({['Normalized Power Spectra: ' epoch ' epoch'],sprintf('Norm to %s %s power: %g-%g Hz',power_method,power_calc,normBand)})
    pause(2)
    saveas(genoSpecFig,[gSpecDir epoch '_spectra'],'svg')
    pause(2)
    close(genoSpecFig)

    animSpecFig = figure();
    setFigureProperties(animSpecFig,'Position',[1 1 2100 1200])
    for k=1:Nstates
        for l=1:numel(animGroups)
            subplot(Nstates,numel(animGroups),(k-1)*numel(animGroups)+l)
            idx = (tenSecIdx(:,1)==l & tenSecIdx(:,3)==k);
            plot(freqVec,10*log10(tenSecNormSpec(idx,:)),'Color',[.7 .7 .7 .75],'LineWidth',2)
            hold on
            plot(freqVec,10*log10(mean(tenSecNormSpec(idx,:))),plotColors(k),'LineWidth',2)
            if k==1
                title(animGroups{l})
            end
            if l==1
                ylabel({'',stateGroups{k}})
            end
            if l==1 && k==ceil(Nstates/2)
                ylabel({'Normalized Power: log scaled',stateGroups{k}})
            end
            if k==Nstates && l==ceil(numel(animGroups)/2)
                xlabel('Frequency (Hz)')
            end
            xlim(freqRange)
        end
    end
    suptitle({['Normalized Power Spectra: ' epoch ' epoch'],sprintf('Norm to %s %s power: %g-%g Hz',power_method,power_calc,normBand)})
    pause(2)
    saveas(animSpecFig,[aSpecDir epoch '_spectra'],'svg')
    pause(2)
    close(animSpecFig)
end
