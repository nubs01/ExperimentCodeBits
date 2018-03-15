% prelim analysis (RZ2 & RW2)

% So I'm trying to see if
% (delta_ket-delta_sleep) > (delta_sal-delta_sleep) -- for both df1 & wt,
% (delta_ket-delta_sal)_df1 > (delta_ket-delta_sal)_wt,
% (theta_sleep)_df1 == (theta_sleep)_wt and
% (delta_sleep/thea_sleep)_df1 > (delta_sleep/theta_sleep)_wt



%% Experiment metadata
animals = {'RZ2','RW2'};
genotypes = {'Df1','WT'};
analysisDays = {[7 8 9 10 11 12],[3 4 5 6 7]};
timeLim = 29; % minutes of each eopch to consider when getting total distance
nEpochs = 3;
dayTypes = {{'sleep','clear','saline'},{'sleep','clear','ketamine'}};
if ~exist('projectDir','var')
    projectDir = uigetdir('Project Directory');
end
currDir = pwd;
cd(projectDir)
expDirs = cell(1,numel(animals));
 for i=1:numel(animals),
     saveFile = sprintf('%s_%s_Experiment/%s_direct/%srntask.mat',animals{i},genotypes{i},animals{i},animals{i}); 
     if exist(saveFile,'file')
         continue;
     end
     days = analysisDays{i};
     rntask = cell(1,max(days));
     for j=1:numel(days),
         rntask{days(j)} = cell(1,nEpochs);
         epochs = dayTypes{mod(days(j),2)+1};
         d = days(j);
         for k=1:nEpochs,
             rntask{d}{k} = struct('animal',animals{i},'genotype',genotypes{i},'epoch',epochs{k},'expDir',[animals{i} '_' genotypes{i} '_Experiment' filesep]);
         end
     end
     save(saveFile,'rntask');
 end

% data groups: WT-Sleep, WT-Clear, WT-Saline, WT-Ketamine, Df1-Sleep, Df1-Clear, Df1-Saline, Df1-Ketamine
% measures: delta power during immobile periods, total distance moved during epoch, average velocity during epoch

% one cell array for all data. Columns: 1-Animal, 2-genotype, 3-day, 4-epoch_type,
% 5-awake_immobile_times (array,use 5Hz low-pass on velocity), 6-immobile_delta
% (array: bin_center,delta_power), 7-sleep_times, 8-sleep_delta, 9-total_dist over
% first 29 minutes (#), 10-[avg_vel std_vel] in first 29 min (2 element vector,
% use low-passed velocity) 
allData = cell(nEpochs*sum(cellfun(@numel,analysisDays)),10); % requires same number of epochs on all days being analyzed
idx = 1;

for a = 1:numel(animals)
    fprintf('Processing animal %s...\n',animals{a})
    anim = animals{a};
    genotype = genotypes{a};
    expDir = [projectDir filesep anim '_' genotype '_Experiment' filesep];
    expDirs{a} = expDir;
    dataDir = [expDir anim '_direct' filesep];
    saveFile = [dataDir anim 'analysisDataframe.mat'];
    if exist(saveFile,'file')
       analysisDataframe = load(saveFile);
       analysisDataframe = analysisDataframe.analysisDataframe;
       allData(idx:idx+size(analysisDataframe,1)-1,:) = analysisDataframe;
       idx = idx+size(analysisDataframe,1);
       clear analysisDataframe dataDir expDir genotype anim
       continue;
   end
    task = load([dataDir anim 'rntask.mat']);
    task = task.rntask;
    tetinfo = load([dataDir anim 'tetinfo.mat']);
    tetinfo = tetinfo.tetinfo;
    riptets = cellfetch(tetinfo,'descrip');
    ripIdx = find(strcmp(riptets.values,'riptet'));
    riptets = riptets.index(ripIdx,:);

    % loop through days and analyze each
    for d = analysisDays{a}
        fprintf('Processing day %02i...\n',d)
        nE = numel(task{d});
        pos = load(sprintf('%s%spos%02i.mat',dataDir,anim,d));
        pos = pos.pos;
        for e=1:nE,
            fprintf('Processing Epoch %02i...\n',e)
            allData(idx,:) = {anim,genotype,d,task{d}{e}.epoch,[],[],[],[],nan,nan};
            posDat = pos{d}{e};

            % Low-pass filter velocity
            [b,a] = butter(4,5/15,'low');
            velVec = posDat.data(:,5);
            if ~all(isfinite(velVec))
                infIdx = ~isfinite(velVec);
                fprintf('Day %02i, epoch %02i: Interpolating to fill %i non-finite velocity points...\n',d,e,sum(infIdx));
                T = 1:numel(velVec);
                velVec(infIdx) = interp1(T(~infIdx),velVec(~infIdx),T(infIdx));
            end
            smoothVel = filtfilt(b,a,velVec);

            % get immobile times, sleep times and awake-immobile times (based on Fisher et al 2012, sleep = immobile for >= 40 sec)
            immTimes = getImmobileTimes(posDat,'velDat',smoothVel,'minDur',5,'thresh',3); % gives 2 col array with start and end times of all immobile (v<=thresh) periods
            diffTimes = diff(immTimes,1,2);
            sleepIdx = diffTimes>=40;
            sleepTimes = immTimes(sleepIdx,:);
            wakeImmTimes = immTimes(~sleepIdx,:);
            allData{idx,5} = wakeImmTimes;
            allData{idx,7} = sleepTimes;

            % get total distance moved and average velocity
            timeIdx =(posDat.data(:,1)<=(posDat.data(1,1)+60*timeLim)); 
            allData{idx,9} = getTotalDistance(posDat.data(timeIdx,2),posDat.data(timeIdx,3));
            allData{idx,10} = [mean(smoothVel(timeIdx)) std(smoothVel(timeIdx))];

            % for each riptet get mean and std of delta power during each
            % immobile time, use delta (col 1&2) and deltaref (col 3&4)
            tets = riptets((riptets(:,1)==d & riptets(:,2)==e),:);
            nTets = size(tets,1);
            for t = tets(:,3)'
                fprintf('Processing tetrode %02i...\n',t)
                delta = load(sprintf('%sEEG%sdelta%02i-%02i-%02i.mat',dataDir,[filesep anim],d,e,t));
                delta = delta.delta{d}{e}{t};
                deltaTime = delta.starttime:1/delta.samprate:delta.endtime;
                deltaPow = double(delta.data(:,3)).^2;
                deltaMat = zeros(size(immTimes,1),2);
                for it = 1:size(immTimes,1)
                    deltaIdx = find(deltaTime>=immTimes(it,1) & deltaTime<=immTimes(it,2));
                    deltaMat(it,:) = [mean(deltaPow(deltaIdx)) std(deltaPow(deltaIdx))];
                end
                wakeImmDat.delta = deltaMat(~sleepIdx,:);
                sleepDat.delta = deltaMat(sleepIdx,:);

                deltaref = load(sprintf('%sEEG%sdeltaref%02i-%02i-%02i.mat',dataDir,[filesep anim],d,e,t));
                deltaref = deltaref.deltaref{d}{e}{t};
                deltarefTime = deltaref.starttime:1/deltaref.samprate:deltaref.endtime;
                deltarefPow = double(deltaref.data(:,3)).^2;
                deltarefMat = zeros(size(immTimes,1),2);
                for it = 1:size(immTimes,1)
                    deltarefIdx = find(deltarefTime>=immTimes(it,1) & deltarefTime<=immTimes(it,2));
                    deltarefMat(it,:) = [mean(deltarefPow(deltarefIdx)) std(deltarefPow(deltarefIdx))];
                end
                wakeImmDat.deltaref = deltarefMat(~sleepIdx,:);
                sleepDat.deltaref = deltarefMat(sleepIdx,:);

                allData{idx,6}{t} = wakeImmDat;
                allData{idx,8}{t} = sleepDat;
                clear sleepDat wakeImmDat delta deltaref deltaMat deltarefMat deltaTime deltarefTime deltaPow deltarefPow deltaIdx deltarefIdx
            end

            % reset loop
            idx = idx+1;
            clear timeIdx sleepIdx immTimes diffTimes sleepTimes wakeImmTimes tets nTets smoothVel posDat
        end
        clear pos nE
    end
    datIdx = find(strcmp(allData(:,1),anim));
    analysisDataframe = allData(datIdx,:);
    fprintf('Saving animal data!\n')
    save(saveFile,'analysisDataframe')
    clear anim riptets ripIdx task tetinfo genotype expDir dataDir analysisDataframe datIdx
end

%% Make plots
set(0,'defaultaxesfontsize',14)
set(0,'defaultaxesfontname','Arial')

% save all plots in RZ2_analysis for now
saveDir = sprintf('%s%sRZ2_Df1_Experiment%sRZ2_analysis%s',projectDir,filesep,filesep,filesep);
genotypes = {'Df1','WT'};
epoch_types = {'sleep','clear','saline','ketamine'};

% comparison groups for delta power
comparisonGroups = {[1 1],[2 1];...
                    [1 2],[2 2];...
                    [1 3],[1 1];...
                    [1 3],[1 4];...
                    [2 3],[2 1];...
                    [2 3],[2 4]};

% for distance pairwise compare all epochs within each genotype; also look at session by session change in move
% compare with t-tests

% Make bar plot of total distance moved during each type of epoch
% make 2 matrices for mean distance and std distance where rows are genotype and columns are epoch type
gt = allData(:,2);
et = allData(:,4);
distMat = zeros(numel(genotypes),numel(epoch_types));
stdDistMat = distMat;
anovaDatMat = [];
anovaGenMat = {};
anovaEpoMat = {};
for i=1:numel(genotypes)
    for j=1:numel(epoch_types)
        idx = find(strcmp(gt,genotypes{i}) & strcmp(et,epoch_types{j}));
        tmpdat = [allData{idx,9}];
        if isrow(tmpdat)
                tmpdat=tmpdat';
        end
        distMat(i,j) = mean(tmpdat);
        stdDistMat(i,j) = std(tmpdat);
        if numel(tmpdat)<3
            stdDistMat(i,j)=0;
        end
        anovaDatMat = [anovaDatMat;tmpdat];
        anovaGenMat = [anovaGenMat;repmat(genotypes(i),numel(tmpdat),1)];
        anovaEpoMat = [anovaEpoMat;repmat(epoch_types(j),numel(tmpdat),1)];
    end
end

% t-tests
pairs = nchoosek(1:numel(epoch_types),2);
distTMat = cell(1,numel(genotypes));
for i=1:numel(genotypes)
    distTMat{i} = ones(size(pairs,1),1);
    for l=1:size(pairs,1)
        j = pairs(l,1);
        k = pairs(l,2);
        xIdx = find(strcmp(gt,genotypes{i}) & strcmp(et,epoch_types{j}));
        xx = [allData{xIdx,9}];
        yIdx = find(strcmp(gt,genotypes{i}) & strcmp(et,epoch_types{k}));
        yy = [allData{yIdx,9}];
        [h,distTMat{i}(l)] = ttest2(xx,yy,'Vartype','unequal');
    end
end

[hb,eb] = errorBarPlot(1:numel(genotypes),distMat,stdDistMat);
a = unique([hb.XData]);
b = [hb.XOffset];
c = [];
for d=1:numel(a)
    c = [c a(d)+b];
end
xCoord = c;
xPairs = cell(1,numel(genotypes)*size(pairs,1));
xPvals = ones(1,numel(xPairs));
k=1;
for i=1:numel(genotypes)
    for j=1:size(pairs,1)
        i1 = pairs(j,1);
        i2 = pairs(j,2);
        xPairs{k} = [xCoord((i-1)*numel(epoch_types)+i1) xCoord((i-1)*numel(epoch_types)+i2)];
        xPvals(k) = distTMat{i}(j);
        k=k+1;
    end
end
validIdx = xPvals<0.05;
sigstar(xPairs(validIdx),xPvals(validIdx))
ax = gca;
f1 = gcf;
ax.XTick = [0 1 1.5 2 2.5];
ax.XTickLabel = {'','Df1_R_Z_2','','WT_R_W_2',''};
ylabel('Total Distance Moved (cm)')
title({'Total Activity in 29 min period','average over days'})
legend(epoch_types)
saveas(f1,[saveDir 'RZ2-RW2_TotalDistance'],'svg')
close all


% Collect awake-immobile delta data
deltaType = 'delta';
groups1 = []; % tags for animal-epoch-TetXX
groups2 = []; % tags for genotype-epoch
deltaPow = [];

for i=1:size(allData,1)
    tmpDat = allData{i,6};
    tetDat = cellfetch(tmpDat,deltaType);
    tetVals = tetDat.values;
    tetIdx = tetDat.index;
    for j=1:size(tetIdx,1)
        tag = sprintf('%s-%s-Tet%02i',allData{i,1},allData{i,4},tetIdx(j));
        tag2 = sprintf('%s-%s',allData{i,2},allData{i,4});
        deltaDat = tetVals{j};
        deltaPow = [deltaPow;deltaDat(:,1)];
        groups1 = char(groups1,repmat(tag,size(deltaDat,1),1));
        groups2 = char(groups2,repmat(tag2,size(deltaDat,1),1));
        if isempty(strtrim(groups1(1,:)))
            groups1 = groups1(2:end,:);
        end
        if isempty(strtrim(groups2(1,:)))
            groups2 = groups2(2:end,:);
        end
    end
end
f2 = figure();
ph = boxplot(deltaPow,groups1,'labelorientation','inline');
ylabel('Delta Power')
title({'Mean Delta Power','during Awake Immobile Periods'})
f3 = figure();
ph2 = boxplot(deltaPow,groups2);
ylabel('Delta Power')
title({'Mean Delta Power','during Awake Immobile Periods'})

saveas(f2,[saveDir 'RZ2-RW2_Delta-perTet'],'svg')
saveas(f3,[saveDir 'RZ2-RW2_Delta'],'svg')
% close all

% deltaType = 'deltaref';
% groups1 = []; % tags for animal-epoch-TetXX
% groups2 = []; % tags for genotype-epoch
% deltaPow = [];
% 
% for i=1:size(allData,1)
%     tmpDat = allData{i,6};
%     tetDat = cellfetch(tmpDat,deltaType);
%     tetVals = tetDat.values;
%     tetIdx = tetDat.index;
%     for j=1:size(tetIdx,1)
%         tag = sprintf('%s-%s-Tet%02i',allData{i,1},allData{i,4},tetIdx(j));
%         tag2 = sprintf('%s-%s',allData{i,2},allData{i,4});
%         delDat = tetVals{j};
%         deltaPow = [deltaPow;deltaDat(:,1)];
%         groups1 = char(groups1,repmat(tag,size(deltaDat,1),1));
%         groups2 = char(groups2,repmat(tag2,size(deltaDat,1),1));
%         if isempty(strtrim(groups1(1,:)))
%             groups1 = groups1(2:end,:);
%         end
%         if isempty(strtrim(groups2(1,:)))
%             groups2 = groups2(2:end,:);
%         end
%     end
% end
% ph = boxplot(deltaPow,groups1);
% f4 = gcf;
% ylabel('Delta Power')
% title({'Mean Referenced Delta Power','during Awake Immobile Periods'})
% ph2 = boxplot(deltaPow,groups2);
% f5 = gcf;
% ylabel('Delta Power')
% title({'Mean Referenced Delta Power','during Awake Immobile Periods'})
% 
% saveas(f4,[saveDir 'RZ2-RW2_Deltaref-perTet'],'svg')
% saveas(f5,[saveDir 'RZ2-RW2_Deltaref'],'svg')
% close all

% Sleep delta
deltaType = 'delta';
groups1 = []; % tags for animal-epoch-TetXX
groups2 = []; % tags for genotype-epoch
sleepDeltaPow = [];

for i=1:size(allData,1)
    tmpDat = allData{i,8};
    tetDat = cellfetch(tmpDat,deltaType);
    tetVals = tetDat.values;
    tetIdx = tetDat.index;
    for j=1:size(tetIdx,1)
        tag = sprintf('%s-%s-Tet%02i',allData{i,1},allData{i,4},tetIdx(j));
        tag2 = sprintf('%s-%s',allData{i,2},allData{i,4});
        deltaDat = tetVals{j};
        sleepDeltaPow = [sleepDeltaPow;deltaDat(:,1)];
        groups1 = char(groups1,repmat(tag,size(deltaDat,1),1));
        groups2 = char(groups2,repmat(tag2,size(deltaDat,1),1));
        if isempty(strtrim(groups1(1,:)))
            groups1 = groups1(2:end,:);
        end
        if isempty(strtrim(groups2(1,:)))
            groups2 = groups2(2:end,:);
        end
    end
end
f2 = figure();
ph = boxplot(sleepDeltaPow,groups1,'labelorientation','inline');
ylabel('Delta Power')
title({'Mean Delta Power','during Awake Immobile Periods'})
f3 = figure();
ph2 = boxplot(sleepDeltaPow,groups2);
ylabel('Delta Power')
title({'Mean Delta Power','during Awake Immobile Periods'})

