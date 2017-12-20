% metadata
posFs = 30;
 animIDs = {'RZ3','RW3','RZ2','RW2'};
 genotypes = {'Df1','WT','Df1','WT'};
winSize = 10; % sec for analyzing spectra
epochTypes = {'sleep','clear','saline','ketamine'};

% Define behavioral states
behavioralStates = {'sleeping','resting','moving'};
stateVelRanges = {[0 .3],[0 .3],[.3 30]}; % cm/sec - need to look at states found to check limits
% try to get dataset from paper and test sleep scoring method
stateDurationRanges = {[40 1800],[10 40],[10 1800]}; % seconds
allowedInterrupts = {0,0,2}; % seconds outside of velocity thresh allowed in an  episode (to account for breif stops while moving and lags in pos tracking)
stateParams = struct('state',behavioralStates,'velRange',stateVelRanges,'durRange',stateDurationRanges,'allowedInterrupt',allowedInterrupts);
[vfB,vfA] = butter(5,5/(posFs/2),'low');

% Experiment Directories
projDir = '/home/roshan/Projects/rn_Schizophrenia_Project/';
expDirs = strcat(projDir,animIDs);
expDirs = strcat(expDirs,strcat('_',genotypes));
expDirs = strcat(expDirs,['_Experiment' filesep]);
dataDirs = strcat(expDirs,strcat(animIDs,['_direct' filesep]));
analysisDirs = strcat(expDirs,strcat(animIDs,['_analysis' filesep]));
figSubDirs = {'Velocity_Histograms','Duration_Histograms','State_Spectra','State_Heatmaps_Hi','State_Heatmaps_Lo'};
prelimDir = [projDir 'RZ-RW_PrelimAnalysisFigs' filesep];
if ~exist(prelimDir,'dir')
    mkdir(prelimDir);
end
prelimSubDirs = {'Foldchange_Plots','NormalizedBaseline_Plots'};
% setup variables to keep
stateMats = cell(1,numel(animIDs));
lowSpectra = cell(1,numel(animIDs));
highSpectra = cell(1,numel(animIDs));
highFreqBaseline = cell(1,numel(animIDs));
hiFreqBand = [40 150];
loFreqBand = [0 15];
baselineBand = [100 150];
dataIdxs = [];
dataEpochs = {};
prelimDataFile = [prelimDir 'prelimDataset-171213.mat'];
skipPlot = 0;

if ~exist(prelimDataFile,'file')
    % Begin looping over animals
    for k=1:numel(animIDs)
        fprintf('Processing animal %s\n',animIDs{k});
        
        % load task and tetinfo structures to get days and tets and epochs
        dDir = dataDirs{k};
        aDir = analysisDirs{k};
        anim = animIDs{k};
        geno = genotypes{k};
        task = load([dDir anim 'rntask.mat']);
        task = task.rntask;
        tetinfo = load([dDir anim 'tetinfo.mat']);
        tetinfo = tetinfo.tetinfo;
        
        % Make figure sub-directories
        for l=1:numel(figSubDirs)
            if ~exist([aDir figSubDirs{l}],'dir')
                mkdir([aDir figSubDirs{l}])
            end
        end
        epochDat = cellfetch(task,'epoch');
        days = unique(epochDat.index(:,1))';
        tetDes = cellfetch(tetinfo,'descrip');
        tetArea = cellfetch(tetinfo,'area');
        stateMats{k} = cell(1,max(days));
        lowSpectra{k} = cell(1,max(days));
        highSpectra{k} = cell(1,max(days));
        highFreqBaseline{k} = cell(1,max(days));
        
        t1 = strcmp(tetDes.values,'riptet');
        t2 = strcmp(tetArea.values,'CA1');
        t3 = intersect(tetDes.index(t1,:),tetArea.index(t2,:),'rows');
        t4 = arrayfun(@(x) any(days==x),t3(:,1));
        tetIdx = t3(t4,:);
        clear t1 t2 t3 t4
        
        for dd=days
            fprintf('   - Day %02i\n',dd);
            tets = unique(tetIdx(tetIdx(:,1)==dd,3))';
            
            % load position data & get behavioral states
            pos = load(sprintf('%s%spos%02i.mat',dDir,anim,dd));
            pos = pos.pos;
            ripples = load(sprintf('%s%sripples%02i.mat',dDir,anim,dd));
            ripples = ripples.ripples{dd};
            epochs = unique(tetIdx(tetIdx(:,1)==dd,2))';
            stateMats{k}{dd} = cell(1,numel(epochTypes));
            lowSpectra{k}{dd} = cell(1,numel(epochTypes));
            highSpectra{k}{dd} = cell(1,numel(epochTypes));
            highFreqBaseline{k}{dd} = cell(1,numel(epochTypes));
            
            % get behavioral states and plot velocity and duration histograms
            for ee=epochs
                epochName = task{dd}{ee}.epoch;
                fprintf('       - Epoch %02i: %s\n',ee,epochName);
                ei = find(strcmp(epochTypes,epochName));
                posDat = pos{dd}{ee};
                fVel = filtfilt(vfB,vfA,interpInf(posDat.data(:,5)));
                stateMat = getBehavioralEpisodes(posDat,stateParams,fVel);
                
                
                if ~skipPlot
                    % plot velocity and duration histograms
                    durHist = plotDurationHistogram(stateMat(:,1:2),behavioralStates(stateMat(:,3)));
                    title({'Duration of behavioral episodes',sprintf('%s Day %02i - %s epoch',anim,dd,epochName)})
                    saveas(durHist,sprintf([aDir 'Duration_Histograms' filesep anim 'epDurHist%02i-%02i'],dd,ee),'svg')
                    close(durHist)
                    
                    posTime = posDat.data(:,1);
                    avgVel = arrayfun(@(x,y) mean(fVel(posTime>=x & posTime<=y)),stateMat(:,1),stateMat(:,2));
                    velHist = plotVelocityHistogram(avgVel,behavioralStates(stateMat(:,3)));
                    title({'Average velocity of behavioral episodes',sprintf('%s Day %02i - %s epoch',anim,dd,epochName)})
                    saveas(velHist,sprintf([aDir 'Velocity_Histograms' filesep anim 'epVelHist%02i-%02i'],dd,ee),'svg')
                    close(velHist)
                end
                lowSpectra{k}{dd}{ei} = cell(1,max(tets));
                highSpectra{k}{dd}{ei} = cell(1,max(tets));
                highFreqBaseline{k}{dd}{ei} = cell(1,max(tets));
                
                stateMats{k}{dd}{ei}.stateTimes = stateMat(:,1:2);
                stateMats{k}{dd}{ei}.stateIdx = stateMat(:,3);
                stateMats{k}{dd}{ei}.states = behavioralStates(stateMat(:,3));
                clear fVel avgVel posTime velHist durHist posDat
                
                % Loop through tetrodes and make spectrum plots
                stateSpectraPlot = figure();
                setFigureProperties(stateSpectraPlot,'Position',[1 1 1440 805]);
                plotRow = 1;
                
                for tt=tets
                    fprintf('           - Tet %02i\n',tt)
                    eeg = load(sprintf('%s%seeg%02i-%02i-%02i.mat',[dDir 'EEG' filesep],anim,dd,ee,tt));
                    eeg = eeg.eeg;
                    eeg = eeg{dd}{ee}{tt};
                    fs = eeg.samprate;
                    eegTime = eeg.starttime:1/fs:eeg.endtime;
                    ripDat = ripples{ee}{tt};
                    
                    [loSpec,loFreq] = getStateSpectra(eeg.data,eegTime,fs,stateMat,loFreqBand,winSize);
                    [hiSpec,hiFreq] = getStateSpectra(eeg.data,eegTime,fs,stateMat,hiFreqBand,winSize);
                    [baselineMean,baselineSTD] = getBaselinePower(eeg.data,eegTime,fs,baselineBand,winSize,[ripDat.starttime ripDat.endtime]);
                    if ~skipPlot
                        % Plot Heatmap
                        figName = @(x,y) sprintf('%s%s%s%s%s%02i-%02i-%02i',aDir,x,filesep,anim,y,dd,ee,tt);
                        loSpecHeatmap = plotAllStateSpectra(loSpec,loFreq,stateMats{k}{dd}{ei}.states,behavioralStates);
                        title({'Power Spectra for All Behavioral Episodes',sprintf('%s Day %02i Tet %02i - %s: Low Freq',anim,dd,tt,epochName)})
                        saveas(loSpecHeatmap,figName('State_Heatmaps_Lo','loStateHeatmap'),'svg')
                        close(loSpecHeatmap)
                        
                        hiSpecHeatmap = plotAllStateSpectra(hiSpec,hiFreq,stateMats{k}{dd}{ei}.states,behavioralStates);
                        title({'Power Spectra for All Behavioral Episodes',sprintf('%s Day %02i Tet %02i - %s: High Freq',anim,dd,tt,epochName)})
                        saveas(hiSpecHeatmap,figName('State_Heatmaps_Hi','hiStateHeatmap'),'svg')
                        close(hiSpecHeatmap)
                    end
                    lowSpectra{k}{dd}{ei}{tt}.spectra = loSpec;
                    lowSpectra{k}{dd}{ei}{tt}.freqs = loFreq;
                    
                    highSpectra{k}{dd}{ei}{tt}.spectra = hiSpec;
                    highSpectra{k}{dd}{ei}{tt}.freqs = hiFreq;
                    
                    highFreqBaseline{k}{dd}{ei}{tt} = struct('mean',baselineMean,'std',baselineSTD);
                    
                    dataIdxs = [dataIdxs;k dd ei tt];
                    dataEpochs{end+1} = epochName;
                    
                    % Plot state spectra for tet
                    if ~skipPlot
                        figure(stateSpectraPlot)
                        plotStateSpectra(loSpec,loFreq,stateMats{k}{dd}{ei}.states,stateSpectraPlot,plotRow,numel(tets),tt,behavioralStates);
                        plotRow = plotRow+1;
                    end
                    
                    clear eeg fs eegTime ripDat loSpec loFreq hiSpec hiFreq baselineMean baselineSTD figName loSpecHeatmap hiSpecHeatmap
                end
                if ~skipPlot
                    figure(stateSpectraPlot)
                    suptitle({'Average Spectra for Each State',sprintf('%s Day %02i - %s',anim,dd,epochName)})
                    saveas(stateSpectraPlot,sprintf('%sState_Spectra%s%sstateSpectra%02i-%02i',aDir,filesep,anim,dd,ee),'svg')
                    close(stateSpectraPlot)
                end
                clear plotRow stateSpectraPlot epochName ei
            end
            
            
            clear dd tets pos
        end
    end
    
    % Save consolidated prelim data
    save([prelimDir 'prelimDataset-' datestr(datenum(date),'yymmdd') '.mat'],'lowSpectra','highSpectra','stateMats','dataIdxs','dataEpochs','highFreqBaseline','animIDs','genotypes','winSize','epochTypes','behavioralStates','stateParams','projDir','expDirs','dataDirs','analysisDirs','figSubDirs','prelimDir','prelimSubDirs','hiFreqBand','loFreqBand','baselineBand')
else
    load(prelimDataFile)
    fprintf('Loaded Prelim Data\n')
end

fprintf('Making Consolidated Figures\n')
%% Make delta and theta box plots
comparisonGenotytpes = {'WT','Df1'};
% Make figure sub-directories
for l=1:numel(prelimSubDirs)
    if ~exist([prelimDir prelimSubDirs{l}],'dir')
        mkdir([prelimDir prelimSubDirs{l}])
    end
end

% Plot and compare change in power between behavioral states, comparing between genotype within epoch type
epochGroups = {'sleep','saline','ketamine','clear'};
stateGroups = {'sleeping','resting','moving'};
diffStates = {{'sleeping','resting'},{'moving','resting'}};
bandNames = {'delta','theta'};
genoGroups = {'WT','Df1'};
freqBands = {[2 5],[5 10]};

getMatchFunc = @(x,y) arrayfun(@(z) any(y==z),x);
XcontainsY = @(x,y) (all(x(1)<=y) & all(x(2)>=y));
capFirst = @(x) [upper(x(1)) lower(x(2:end))];

for k=1:numel(diffStates)
    s1 = diffStates{k}{1};
    % si1 = find(strcmp(behavioralStates,s1));
    s2 = diffStates{k}{2};
    % si2 = find(strcmp(behavioralStates,s2));
    
    for bb = 1:numel(bandNames)
        % Collect data to plot and compare
        bN = bandNames{bb};
        freqRange = freqBands{bb};
        nGroups = numel(genoGroups)*numel(epochGroups);
        groupDat = cell(numel(epochGroups),numel(genoGroups)); % significance test between genotypes for each epochType
        groupStats = cell(numel(epochGroups),1); % to store p-values and such
        flatDat = [];
        flatDatGroups = {};
        fName = [prelimDir 'Foldchange_Plots' filesep s1 'V' s2 '_' bN];
        
        for ee=1:numel(epochGroups)
            epochName = epochGroups{ee};
            ei = find(strcmp(epochTypes,epochName));
            % eIdx = dataIdxs(:,3)==ei;
            
            for gg=1:numel(genoGroups)
                geno = genoGroups{gg};
                anims = find(strcmp(genotypes,geno));
                % gIdx = getMatchFunc(dataIdxs(:,1),anims);
                % dIdx = dataIdxs((gIdx & eIdx),:);
                
                % Get spectra for each state to be differenced
                if XcontainsY(loFreqBand,freqRange)
                    S = cellfetch(lowSpectra,'spectra');
                    F = cellfetch(lowSpectra,'freqs');
                elseif XcontainsY(hiFreqBand,freqRange)
                    S = cellfetch(highSpectra,'spectra');
                    F = cellfetch(highSpectra,'freqs');
                end
                if isempty(S.values) || all(cellfun(@isempty,S.values))
                    error('No appropriate spectra found...Ahhhh!')
                end
                freqs = F.values{1};
                if isempty(freqs)
                    error('Well thats fucked up....no freqs')
                end
                sIdx = S.index;
                spectra = S.values;
                ri = (getMatchFunc(sIdx(:,1),anims) & sIdx(:,3)==ei);
                sIdx = sIdx(ri,:);
                spectra = spectra(ri);
                
                clear S F
                
                spec1 = [];
                spec2 = [];
                freqIdx = find(freqs>=freqRange(1) & freqs<=freqRange(2));
                
                for l=1:numel(spectra)
                    states = stateMats{sIdx(l,1)}{sIdx(l,2)}{sIdx(l,3)}.states;
                    s1Idx = find(strcmp(states,s1));
                    s2Idx = find(strcmp(states,s2));
                    spec1 = [spec1;mean(spectra{l}(s1Idx,freqIdx),1)];
                    spec2 = [spec2;mean(spectra{l}(s2Idx,freqIdx),1)];
                end
                groupDat{ee}{gg} = (mean(spec1,2)-mean(spec2,2))./mean(spec2,2);
                flatDat = [flatDat;groupDat{ee}{gg}];
                flatDatGroups = [flatDatGroups;repmat({[geno '_' epochName]},numel(groupDat{ee}{gg}),1)];
            end
            [~,p,ci,stats] = ttest2(groupDat{ee}{1},groupDat{ee}{2}); % HARDCODED: for 2 genotypes to compare
            groupStats{ee}.P = p;
            groupStats{ee}.CI = ci;
            groupStats{ee}.Stats = stats;
        end
        
        % Make Figure
        sigPairs = arrayfun(@(x) [2*x-1 2*x],1:numel(groupStats),'UniformOutput',false);
        P = cellfetch(groupStats,'P');
        P = [P.values{sort(P.index)}];
        P(P>0.05) = nan;
        foldChangePlot = figure();
        setFigureProperties(foldChangePlot,'Position',[1 1 1440 400])
        boxplot(flatDat,flatDatGroups)
        set(gca,'TickLabelInterpreter','tex','XTickLabel',strcat(strrep(get(gca,'XTickLabel'),'_','_{'),'}'));
        xlabel('Genotype x Epoch Type')
        ylabel(['Fold change from ' s2])
        title({sprintf('%s vs %s',capFirst(s1),capFirst(s2)),sprintf('Fold Change in %s (%i-%i Hz)',capFirst(bN),freqRange(1),freqRange(2)),strjoin(animIDs,', ')})
        sigstar(sigPairs,P)
        pause(2)
        saveas(foldChangePlot,fName,'svg')
        save([fName '.mat'],'freqRange','groupDat','groupStats','bN','s1','s2','epochGroups','genoGroups','behavioralStates')
        close(foldChangePlot)
        clear sigPair foldChangePlot bN flatDat flatDatGroups groupDat groupStats
    end
    clear s1 s2
end


% Now make boxplot comparing delta & theta normalized to high freq baseline.
% Only use sleep epoch data. Comparison between genotypes within each
% behavioral state
epochName = 'sleep';
ei = 1;
for ei=1:numel(epochGroups)
    epochName = epochTypes{ei};
    for bb=1:numel(bandNames)
        bN = bandNames{bb};
        freqRange = freqBands{bb};
        groupDat = cell(numel(stateGroups),numel(genoGroups));
        groupStats = cell(numel(stateGroups),1);
        flatDat = [];
        flatDatGroups = {};
        fName = [prelimDir 'NormalizedBaseline_Plots' filesep epochName '_norm' capFirst(bN)];
        for ss = 1:numel(stateGroups)
            sg = stateGroups{ss};
            for gg=1:numel(genoGroups)
                geno = genoGroups{gg};
                anims = find(strcmp(genotypes,geno));
                if XcontainsY(loFreqBand,freqRange)
                    S = cellfetch(lowSpectra,'spectra');
                    F = cellfetch(lowSpectra,'freqs');
                elseif XcontainsY(hiFreqBand,freqRange)
                    S = cellfetch(highSpectra,'spectra');
                    F = cellfetch(highSpectra,'freqs');
                end
                if isempty(S.values) || all(cellfun(@isempty,S.values))
                    error('No appropriate spectra found...Ahhhh!')
                end
                freqs = F.values{1};
                if isempty(freqs)
                    error('Well thats fucked up....no freqs')
                end
                sIdx = S.index;
                spectra = S.values;
                baseline = cellfetch(highFreqBaseline,'mean');
                baselineSTD = cellfetch(highFreqBaseline,'std');
                [bIdx,iB1,iB2] = intersect(baseline.index,baselineSTD.index,'rows');
                baseline = baseline.values(iB1);
                baselineSTD = baselineSTD.values(iB2);
                [bsIdx,iS,iB] = intersect(sIdx,bIdx,'rows');
                spectra = spectra(iS);
                baseline = baseline(iB);
                baselineSTD = baselineSTD(iB);
                
                ri = (getMatchFunc(bsIdx(:,1),anims) & bsIdx(:,3)==ei);
                spectra = spectra(ri);
                baseline = baseline(ri);
                baselineSTD = baselineSTD(ri);
                bsIdx = bsIdx(ri,:);
                
                spec = [];
                base = [];
                baseSTD = [];
                freqIdx = find(freqs>=freqRange(1) & freqs<=freqRange(2));
                
                for l=1:numel(spectra)
                    states = stateMats{bsIdx(l,1)}{bsIdx(l,2)}{bsIdx(l,3)}.states;
                    s1Idx = find(strcmp(states,sg));
                    spec = [spec;spectra{l}(s1Idx,freqIdx)];
                    base = [base;repmat(baseline{l},numel(s1Idx),1)];
                    baseSTD = [baseSTD;repmat(baselineSTD{l},numel(s1Idx),1)];
                end
                groupDat{ss}{gg}.values = mean(spec,2)./base;
                groupDat{ss}{gg}.errors = abs(groupDat{ss}{gg}.values).*sqrt((std(spec,0,2)./mean(spec,2)).^2+(baseSTD./base).^2);
                flatDat = [flatDat;groupDat{ss}{gg}.values];
                flatDatGroups = [flatDatGroups;repmat({[geno '_' sg]},numel(groupDat{ss}{gg}.values),1)];
            end
            [~,p,ci,stats] = ttest2(groupDat{ss}{1}.values,groupDat{ss}{2}.values); % HARDCODED
            groupStats{ss}.P = p;
            groupStats{ss}.CI = ci;
            groupStats{ss}.Stats = stats;
        end
        
        % Make Figure
        sigPairs = arrayfun(@(x) [2*x-1 2*x],1:numel(groupStats),'UniformOutput',false);
        P = cellfetch(groupStats,'P');
        P = [P.values{sort(P.index)}];
        P(P>0.05) = nan;
        normPlot = figure();
        setFigureProperties(normPlot,'Position',[1 1 1440 800])
        ax1 = axes('Position',[0 0 1 1],'Visible','off');
        ax2 = axes('Position',[.1 .5 .8 .3]);
        boxplot(flatDat,flatDatGroups)
        set(gca,'TickLabelInterpreter','tex','XTickLabel',strcat(strrep(get(gca,'XTickLabel'),'_','_{'),'}'));
        xlabel('Genotype x Behavioral State Type')
        ylabel('Normalized Band Power')
        title({sprintf('%s (%i-%i Hz) Power normalized to Hi-Frequency (%i-%i Hz) Power',capFirst(bN),freqRange(1),freqRange(2),baselineBand(1),baselineBand(2)),strjoin(animIDs,', '),sprintf('^{%s Epoch}',capFirst(epochName))})
        sigstar(sigPairs,P)
        % TODO: Insert T-Test stats string & test
        pause(2)
        saveas(normPlot,fName,'svg')
        save([fName '.mat'],'groupDat','groupStats','bN','epochName','sg','stateGroups','genoGroups','freqRange','baselineBand')
        close(normPlot)
    end
end
