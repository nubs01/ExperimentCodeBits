nboot = 1000;
genotype = 'WT';
epoch = 'sleep';
groupDefs = struct('genotype',genotype,'behavioral_state',{'sleeping','resting','moving'},'epoch',epoch);
groupNames = {'Sleeping','Resting','Moving'};

projDir = '~/Projects/rn_Schizophrenia_Project/';
baseFigDir = [projDir 'RZ-RW_PrelimAnalysisFigs/'];
baseFigDir = [baseFigDir '10s_stats' filesep];
datasetFN = [baseFigDir 'allPrelimMetrics-10s.mat'];
if ~exist(datasetFN,'file')
    [dataStruct,dataParams] = collectPrelimDataset();
    save(datasetFN,'dataStruct','dataParams')
elseif ~exist('dataStruct','var')
    load(datasetFN)
end
baseFigDir = [baseFigDir 'WT_stateComparison' filesep 'Raw' filesep];
mkdir(baseFigDir)

Ngroups = numel(groupDefs);
grpIdx = cell(1,Ngroups);
grpPow = cell(1,Ngroups);
grpBoot = cell(1,Ngroups);
for k=1:Ngroups
    [grpIdx{k},grpPow{k},grpBoot{k}] = getGroupMetrics(groupDefs(k),dataStruct);
end

grpPairs = nchoosek(1:Ngroups,2);
Npairs = size(grpPairs,1);
bandNames = dataParams.band_names;
bandFreqs = dataParams.band_freqs;
Nbands = numel(dataParams.band_freqs);
plotColors = 'bmkcrgy';
pairStats = cell(1,Nbands);
for l=1:Nbands
    pairStats{l} = nan(Ngroups,Ngroups,2);
end

for l=1:Npairs
    g1 = grpPairs(l,1);
    g2 = grpPairs(l,2);
    fprintf('   Permutation Testing %s vs %s to see if mean band power is significantly different\n',groupNames{g1},groupNames{g2})
    G1_boot  = grpBoot{g1};
    G1_pow   = grpPow{g1};
    G2_boot  = grpBoot{g2};
    G2_pow   = grpPow{g2};
    gn1 = groupNames{g1};
    gn2 = groupNames{g2};

    statFig = figure();
    setFigureProperties(statFig,'Position',[1 1 1440 1200])

    for m = 1:Nbands
        permStat = permutationTest(G1_pow(:,m),G2_pow(:,m),nboot);
        pairStats{m}(g1,g2,1) = G1_boot(m).mean-G2_boot(m).mean;
        pairStats{m}(g1,g2,2) = permStat.P_two_tail;
        plotRow = (m-1)*3;
        figure(statFig);
        subplot(Nbands,3,plotRow+1)
        h = plotPowHist(G1_pow(:,m),G2_pow(:,m),plotColors(g1),plotColors(g2));
        if m==1
            legend(h,groupNames([g1 g2]))
            title('Total Band Power')
        end
        bn = bandNames{m};
        bn = strrep(bn,'-',' & ');
        ylabel({bn,sprintf('%g - %g Hz',bandFreqs{m})})
        if m==Nbands
            xlabel('Total Power')
        end

        subplot(Nbands,3,plotRow+2)
        h = plotPowHist(G1_boot(m).bootstats,G2_boot(m).bootstats,plotColors(g1),plotColors(g2),G1_boot(m),G2_boot(m),10);
        if m==1
            title(sprintf('Boostrap Means, N=%i',nboot))
        end
        if m==Nbands
            xlabel('Mean Power')
        end

        subplot(Nbands,3,plotRow+3)
        plotPermTest(permStat.permutation_stats,G1_boot(m).mean-G2_boot(m).mean,permStat.P_two_tail);
        permStat.mean_diff = G1_boot(m).mean - G2_boot(m).mean;
        if m==1
            title({sprintf('Permutation Test: %s - %s',groupNames{g1},groupNames{g2}),'two-tailed'})
        end
        if m==Nbands
            xlabel('Difference in Means')
        end
    end
    suptitle({sprintf('%s vs %s',groupNames{g1}, groupNames{g2}),[genotype ': ' epoch ' epoch']})
    figName = [baseFigDir genotype '_SpectralComparison-' groupNames{g1} '_vs_' groupNames{g2} '_' epoch];
    saveas(statFig,figName,'svg')
    close(statFig)
    save([figName '.mat'],'permStat','G1_boot','G1_pow','G2_boot','G2_pow','groupNames','bandNames','bandFreqs','g1','g2')
end

relationMat = figure();
setFigureProperties(relationMat,'Position',[1 1 1800 750])
for l=1:Nbands
    subplot(1,Nbands,l)
    cFlag = 0;
    ylFlag = 1;
    if l==Nbands
        cFlag = 1;
    end
    if l>1
        ylFlag = 0;
    end
    h = plotRelationMatrix(pairStats{l},groupNames,ylFlag,cFlag);
    title({[bandNames{l} ' power'],sprintf('%g - %g Hz',bandFreqs{l})})
end
spt = suptitle(sprintf('Diff Mean Power: %s Epoch: %s',capFirst(epoch),genotype));
set(spt,'Position',[.5 -0.018 0])
figName = [baseFigDir 'RelationMatrix_' genotype '_' epoch];
saveas(relationMat,figName,'svg')
pause(2)
close(relationMat)

function h = plotPowHist(X,Y,colorX,colorY,statX,statY,nbins)
    rX = range(X);
    rY = range(Y);
    minR = min([rX,rY]);
    if ~exist('nbins','var')
        nbins = 30;
    end
    binWidth = minR/nbins; % use 40 bins to cover range
    h1 = histogram(X,'BinWidth',binWidth,'FaceColor',colorX);
    hold on
    h2 = histogram(Y,'BinWidth',binWidth,'FaceColor',colorY);
    set(gca,'XLim',get(gca,'XLim').*[.9 1.1])
    if nargin>4
        yMax = max([max(h1.Values),max(h2.Values)]);
        ylim = [0 max([10,ceil(1.2*yMax)])];
        txtY = ceil(max([3,1.05*yMax]));
        set(gca,'YLim',ylim);
        mX = statX.mean;
        seX = statX.SEM;
        mY = statY.mean;
        seY = statY.SEM;
        txtStr = {sprintf('%0.2g \\pm %2.2g',mX,seX),sprintf('%0.2g \\pm %2.2g',mY,seY)};
        plot([mX mX],[0 txtY],colorX,'LineWidth',2);
        plot([mY mY],[0 txtY],colorY,'LineWidth',2);
        if abs(mX-mY)>0.35*range(get(gca,'XLim'))
            text(mX,txtY,txtStr{1},'HorizontalAlignment','center','VerticalAlignment','bottom')
            text(mY,txtY,txtStr{2},'HorizontalAlignment','center','VerticalAlignment','bottom')
        elseif mX>mY
            text(mX,txtY,txtStr{1},'HorizontalAlignment','left','VerticalAlignment','bottom')
            text(mY,txtY,txtStr{2},'HorizontalAlignment','right','VerticalAlignment','bottom')
        else
            text(mX,txtY,txtStr{1},'HorizontalAlignment','right','VerticalAlignment','bottom')
            text(mY,txtY,txtStr{2},'HorizontalAlignment','left','VerticalAlignment','bottom')
        end
    end
    h = [h1 h2];
end

function h = plotPermTest(X,mX,pVal)
    h = histogram(X);
    yMax = max(h.Values);
    txtY = ceil(max([3,1.05*yMax]));
    yMax = ceil(max([8 1.2*yMax]));
    set(gca,'YLim',[0 yMax]);
    hold on
    plot([mX mX],[0 txtY],'k','LineWidth',2)
    if pVal==0 && mX>max(X)
        hAlign = 'right';
    elseif pVal==0 && mX<max(X)
        hAlign = 'left';
    else
        hAlign = 'center';
    end
    text(mX,txtY,sprintf('%0.3g, p=%g',mX,pVal),'HorizontalAlignment',hAlign,'VerticalAlignment','bottom')
end

function h = plotRelationMatrix(pStats,gNames,ylFlag,cFlag)
    colors = {[1 1 1],[0 0 0],[1 .2 0],[0 .5 .9]}; % untested, not significant, significant lower, significant higher
    a = pStats(:,:,1);
    if any(a==0)
    error('Fix this')
    end
    a(isnan(a)) = 0;
    a = a-a';
    a(a==0) = nan;
    b = pStats(:,:,2);
    b(isnan(b)) = 0;
    b = b+b';
    pMat = ones(numel(gNames),numel(gNames));
    pMat(a>0 & b<0.05) = 4;
    pMat(a<0 & b<0.05) = 3;
    pMat(isnan(a)) = 2;
    cMat = zeros(numel(gNames),numel(gNames),3);
    txtX = zeros(1,numel(a));
    txtY = zeros(1,numel(a));
    txtStr = cell(1,numel(a));
    n=1;
    for k=1:numel(gNames)
    for l=1:numel(gNames)
    cMat(k,l,:) = colors{pMat(k,l)};
    txtX(n) = l;
    txtY(n) = k;
    if isnan(a(k,l))
        txtStr{n} = {'',''};
    else
        txtStr{n} = {sprintf('%2.2g',a(k,l)),sprintf('p=%2.2g',b(k,l))};
    end
    n=n+1;
    end
    end
    h = image(cMat);
    h.AlphaData = 0.5;
    centers = 1:numel(gNames);
    edges = centers(1:end-1)+0.5;
    xl = get(gca,'XLim');
    yl = get(gca,'YLim');
    hold on
    plot(xl,[edges;edges],'k')
    plot([edges;edges],yl,'k')
    text(txtX,txtY,txtStr,'HorizontalAlignment','center','VerticalAlignment','middle')
    set(gca,'XTick',centers,'YTick',centers,'XTickLabels',gNames)
    if ylFlag
    set(gca,'YTickLabels',gNames)
    else
    set(gca,'YTickLabels',cell(1,numel(gNames)))
    end

    % create colorbar
    if cFlag
    mycolors = cell2mat(colors([3 2 4])');
    colormap(mycolors)
    axPos = get(gca,'Position');
    cb = colorbar('Position',[axPos(1)+axPos(3)+0.02 axPos(2) 0.03 axPos(4)]);
    caxis([0 1])
    cb.Ticks = [.1 .5 .9];
    cb.TickLabels = {'Y<X','N.S.','Y>X'};
    end
end
