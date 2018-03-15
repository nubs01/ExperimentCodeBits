nboot = 1000;
animIDs = {'RW2','RW3','RZ2','RZ3'};
genotype = {'WT','WT','Df1','Df1'};
epoch = 'sleep';

projDir = '~/Projects/rn_Schizophrenia_Project/';
dataDirs = strcat(projDir,animIDs,'_',genotype,'_Experiment',filesep,animIDs,'_direct',filesep);

tetReq = {'descrip','riptet','area','CA1'};
taskReq = {'epoch',epoch};
immThresh = .3; %cm/s

Nanim = numel(animIDs);
ripCounts = cell(Nanim,1);
ripTimes = cell(Nanim,1);
immRipRate = cell(Nanim,1);
grpLabels = cell(Nanim,1);

for k=1:Nanim
    tetIdx = getTetIdx(dataDirs{k},animIDs{k},tetReq,taskReq);
    % only use one tet for each day/epoch
    [~,b] = unique(tetIdx(:,1:2),'rows');
    tetIdx = tetIdx(b,:);
    velTraces = getAnimVel(dataDirs{k},animIDs{k},tetIdx(:,1:2),'velFs',25,'truncateDat',29);
    ripDat = getAnimRipDat(dataDirs{k},animIDs{k},tetIdx);

    Ntraces = numel(ripDat);
    ripCounts{k} = zeros(Ntraces,1);
    ripTimes{k} = zeros(Ntraces,1);
    immRipRate{k} = zeros(Ntraces,1);

    for l=1:Ntraces
        % get counts
        ripCounts{k}(l) = numel(ripDat(l).starttime);
        ripTimes{k}(l) = diff(ripDat(l).timerange);

        % Find number of ripples that occur during immobile times and total time immobile
        vT = velTraces(l,:);
        vTime = 0:1/25:(numel(vT)-1)/25;
        rST = ripDat(l).starttime - ripDat(l).timerange(1);
        rET = ripDat(l).endtime - ripDat(l).timerange(1);

        imm = vT<=immThresh;
        c = contiguous(imm,1);
        c = c{2};
        c(diff(c,[],2)==0,:)=[];
        immTime = vTime(c);
        tmpCount = 0;
        tmpTime = 0;
        for m=1:size(immTime,2)
            tmpTime = diff(immTime(m,:))+ tmpTime;
            tmpCount = sum(rST>=immTime(m,1) & rET<=immTime(m,2))+ tmpCount;
        end
        immRipRate{k}(l) = tmpCount/tmpTime;
    end
    grpLabels{k} = repmat(animIDs(k),Ntraces,1);
end

labels = vertcat(grpLabels{:});
totRipRate = cell2mat(ripCounts)./cell2mat(ripTimes);
RipRateImm = cell2mat(immRipRate);

% make plots
F1 = figure();
setFigureProperties(F1,'Position',[1 1 1800 1200]);
subplot(1,2,1)
boxplot(totRipRate,labels)
title('Total Ripples per sec over whole epoch')
subplot(1,2,2)
boxplot(RipRateImm,labels)
title('Ripples per sec during immobility')

