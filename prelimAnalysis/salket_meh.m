binTimes = [(0:10:14*60+50)' (10:10:15*60)'];
binTimes = num2cell(binTimes,2);
allEpo = {'sleep','saline','ketamine'};
meanDat = cell(3,numel(allEpo));
semDat = cell(3,numel(allEpo));
for k=1:numel(allEpo)
    epoch = allEpo{k};
    groupDefs = struct('genotype','Df1','epoch',epoch,'segment_time',binTimes);
    grpBoot = cell(1,numel(groupDefs));
    for l=1:numel(groupDefs)
        [~,~,grpBoot{l}] = getGroupMetrics(groupDefs(l),dataStruct);
    end
    trace = cell2mat(grpBoot');
    for m=1:3
        meanDat{m,k} = [trace(:,m).mean];
        semDat{m,k} = [trace(:,m).SEM];
    end
end
a = cell2mat(binTimes);
X1 = mean(a,2)/60;
X2 = X1+15;
tracePlot = figure();
setFigureProperties(tracePlot,'Position',[1 1 1800 1200])
plotColors = 'mgk';
for m=1:3
    for k=1:2
        subplot(2,3,m+3*(k-1))
        hold on
        [~,h] = shadedErrorPlot(X1,meanDat{m,1},semDat{m,1},'b');
        [~,h1] = shadedErrorPlot(X2,meanDat{m,k+1},semDat{m,k+1},plotColors(k));
        if m==2
            legend([h,h1],'Sleep',capFirst(allEpo{k+1}))
        end
    end
end
subplot(2,3,1)
yl = get(gca,'YLim');
