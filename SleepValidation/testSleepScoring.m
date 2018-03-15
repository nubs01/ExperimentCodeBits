animIDs = {'RW2'};
projDir = '~/Projects/rn_Schizophrenia_Project/';
expDirs = strcat(animIDs,'_WT_Experiment/');
dataDirs = strcat(projDir,strcat(expDirs,strcat(animIDs,'_direct/')));
figDir = [projDir 'RZ-RW_PrelimAnalysisFigs/WT-SleepScoring/'];
mkdir(figDir)
testEpochs = {'sleep'};
Nanim = numel(animIDs);

minOfData = 26; % cutoff first and last 2 minutes of data
eegTraces = cell(Nanim,1);
eegrefTraces = cell(Nanim,1);
velTraces = cell(Nanim,1);
datIdx = cell(Nanim,1);

for k=1:Nanim
    tetinfo = load(sprintf('%s%stetinfo.mat',dataDirs{k},animIDs{k}));
    rntask  = load(sprintf('%s%srntask.mat',dataDirs{k},animIDs{k}));
    tetinfo = tetinfo.tetinfo;
    rntask  = rntask.rntask;
    tetIdx = getCellStructIdx(tetinfo,'descrip','riptet','area','CA1');
    taskIdx = getCellStructIdx(rntask,'epoch','sleep');
    days = taskIdx(:,1);
    epochs = taskIdx(:,2);
    tmp1 = arrayfun(@(x,y) any(days==x) & any(epochs==y),tetIdx(:,1),tetIdx(:,2));
    tetIdx = tetIdx(tmp1,:);
    
    eegTraces{k} = getAnimEEG(dataDirs{k},animIDs{k},tetIdx,'truncateDat',28);
    eegrefTraces{k} = getAnimEEG(dataDirs{k},animIDs{k},tetIdx,'truncateDat',28);
    datIdx{k} = [repmat(k,size(tetIdx,1),1) tetIdx];
    velTraces{k} = getAnimVel(dataDirs{k},animIDs{k},tetIdx,'truncateDat',28);
end


eegDat = cell2mat(eegTraces);
velDat = cell2mat(velTraces);
eegrefDat = cell2mat(eegrefTraces);
eegFs = 1500; % TODO: fix hardcoding
velFs = 25; % TODO: fix hardcoding
dataIdx = cell2mat(datIdx);


