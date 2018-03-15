nboot = 1000;
animIDs = {'RW2','RW3','RZ2','RZ3'};
genotype = {'WT','WT','Df1','Df1'};
epoch = 'sleep';
states = {'sleeping','resting','moving'};
groupNames = animIDs;
alpha = 0.05;

projDir = '~/Projects/rn_Schizophrenia_Project/';
baseFigDir = [projDir 'RZ-RW_PrelimAnalysisFigs/'];
baseFigDir = [baseFigDir '2s_stats_ShiftedFreq' filesep];
datasetFN = [baseFigDir 'allPrelimMetrics-2s.mat'];
if ~exist(datasetFN,'file')
    [dataStruct,dataParams] = collectPrelimDataset();
    save(datasetFN,'dataStruct','dataParams')
elseif ~exist('dataStruct','var')
    load(datasetFN)
end
baseFigDir = [baseFigDir 'Kruskal-Wallis-Test' filesep 'NoSigBars' filesep epoch filesep];
mkdir(baseFigDir)

bandFreqs = dataParams.band_freqs;
bandNames = dataParams.band_names;
Nbands = numel(bandNames);

for ss=1:numel(states)
    state = states{ss};
    groupDefs = struct('genotype',genotype,'animal',animIDs,'behavioral_state',state,'epoch',epoch);
    Ngroups = numel(groupDefs);
    grpIdx = cell(Ngroups,1);
    grpPow = cell(Ngroups,1);
    grpBoot = cell(Ngroups,1);
    grpLabels = cell(Ngroups,1);
    for k=1:Ngroups
        [grpIdx{k},grpPow{k},grpBoot{k}] = getGroupMetrics(groupDefs(k),dataStruct);
        grpLabels{k} = repmat(groupNames(k),size(grpPow{k},1),1);
    end

    figDir = [baseFigDir state filesep];
    mkdir(figDir)

    % flatten group stats
    testPow = vertcat(grpPow{:});
    testLabels = vertcat(grpLabels{:});
    P = zeros(1,Nbands);
    ANOV = cell(1,Nbands);
    KWstat = cell(1,Nbands);
    MCstat = cell(1,Nbands);
    MCmeans = cell(1,Nbands);
    MCgnames = cell(1,Nbands);
    for m=1:Nbands
        [P(m),ANOV{m},KWstat{m}] = kruskalwallis(testPow(:,m),testLabels,'off');
        [MCstat{m},MCmeans{m},~,MCgnames{m}] = multcompare(KWstat{m},'Alpha',alpha,'CType','bonferroni');
        F1 = figure();
        F2 = figure();
        setFigureProperties(F1,'defaultAxesFontSize',24)
        setFigureProperties(F2,'defaultAxesFontSize',24)
        plotColors = 'rmbc';%[1 0 0],[.7 .2 0],[0 0 1],[0 .3 .7];
        binWidth = [];
        sigmean = zeros(1,Ngroups);
        for k=1:Ngroups
            figure(F1)
            subplot(2,2,k) % TODO: Fix Hardcoding
            idx = strcmp(testLabels,groupNames{k});
            histogram(testPow(idx,m),'FaceColor',plotColors(k))
            if k==1 || k==3 % TODO: fix hardcoding
                ylabel('Counts')
            end
            if k>2 % TODO: fix hardcoding
                xlabel('Normalized Power')
            end
            title(groupNames{k})
            set(gca,'XLim',[0 1.5])
            figure(F2)
            hold on
            tmpBoot = grpBoot{k}(m);
            if isempty(binWidth)
                h = histogram(tmpBoot.bootstats,'FaceColor',plotColors(k));
                binWidth = h.BinWidth;
            else
                h = histogram(tmpBoot.bootstats,'BinWidth',binWidth,'FaceColor',plotColors(k));
            end
            txtY = max(h.Values)+5;
            plot([1 1]*tmpBoot.mean,[0 txtY],'k','LineWidth',2)
            text(tmpBoot.mean,txtY,sprintf('%.3g \\pm %.3g',tmpBoot.mean,tmpBoot.SEM),...
                'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',24)
            sigmean(k) = tmpBoot.mean;
        end
        figure(F2)
        sigArr = sigmean(MCstat{m}(:,1:2));
        sigP = MCstat{m}(:,end);
        sigArr(sigP>alpha,:) = [];
        sigP(sigP>alpha) = [];
        sigCell = num2cell(sigArr,2);
        %sigstar(sigCell,sigP)
        pause(1)
        figure(F1)
        bn = bandNames{m};
        bn = strrep(bn,'-',' & ');
        suptitle({sprintf('Normalized %s Power',capFirst(bn)),[capFirst(state) ' Episodes']}) %TODO put epoch in title
        figure(F2)
        %title({sprintf('Bootstrapped Mean %s power',capFirst(bn)),[capFirst(state) ' Episodes']})
        xlim([0.1 .25])
        xlabel('Mean Normalized Power')
        ylabel('Counts')

        saveas(F1,sprintf('%sPowHistograms_%s-%s-%s',figDir,epoch,state,bn),'fig')
        saveas(F2,sprintf('%sMeanPow_%s-%s-%s',figDir,epoch,state,bn),'fig')
        close(F1)
        close(F2)
    end
    save([figDir 'KruskalWallis-Stats.mat'],'P','ANOV','KWstat','MCstat','MCmeans','MCgnames','testPow','testLabels')
end   
