groupDefs = [];
grpPow = cell(4,2);
grpBoot = cell(4,2);
foldChange = cell(4,2);
meanFC = cell(4,2);
semFC = cell(4,2);
grpLabels = cell(4,1);
for k=1:4
gD = struct('animal',animIDs{k},'behavioral_state',{'resting','sleeping'},'epoch','sleep');
groupDefs = [groupDefs;gD];
for l=1:2
    [~,grpPow{k,l},grpBoot{k,l}] = getGroupMetrics(gD(l),dataStruct);
end
for m=1:2
    fC = (grpBoot{k,2}(m).bootstats-grpBoot{k,1}(m).bootstats)./grpBoot{k,1}(m).bootstats;
    foldChange{k,m} = fC;
    meanFC{k,m} = (mean(grpPow{k,2}(:,m))-mean(grpPow{k,1}(:,m)))/mean(grpPow{k,1}(:,m));
    semFC{k,m} = std(fC);
    grpLabels{k} = repmat(animIDs(k),size(fC,1),1);
end

end

testLabels = vertcat(grpLabels{:});
testFC = cell2mat(foldChange);
testMeans = cell2mat(meanFC);
testSEM = cell2mat(semFC);

for m=1:2
    F1 = figure();
    setFigureProperties(F1,'defaultAxesFontSize',24)
    bh = errorBarPlot(1:4,testMeans(:,m),testSEM(:,m));
    %boxplot(testFC(:,m),testLabels)
    title(bandNames{m})
end
