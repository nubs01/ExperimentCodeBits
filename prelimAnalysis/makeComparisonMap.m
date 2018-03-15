compNames = {{'genotype','genotype','animal','animal','animal','animal'},{'WT','Df1','RZ3','RZ2','RW3','RW2'}};
epoch = 'sleep';
state = 'resting';
nboot = 
comparisons = struct('genotype',{},'animal',{},'epoch',{},'state',{});

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
baseFigDir = [baseFigDir 'KetamineAnalysis' filesep];
mkdir(baseFigDir)

% Compare all anim/geno in  sleep, per state
% Compare WT sleep vs saline vs ketamine per state : focus on first 10 min of ketamine epoch (1 min bins using segment_time)
% Compare Df1 same
% Repeat for each animal individually 

% Ketamine Analysis
binTimes = [(0:60:14*60)' (60:60:15*60)'];
allStates = {'sleeping','resting','moving'};
allGeno = {'WT','DF1'};
allEpochs = {'sleep','saline','ketamine'};
Ngeno = numel(allGeno);
Nstates = numel(allStates);
relationMats = figure();
setFigureProperties(relationMats,'Position',[1 1 1800 1200]);
tracePlot = figure();
setFigureProperties(tracePlot,'Position',[1 1 1800 700]);

for g=1:numel(allGeno)
    genotype = allGeno{g};


    % plot min by min delta power for each epoch (mean and bootstrap SEM)
    
        




    for s=1:numel(allStates)
        state = allStates{s};

        groupDefs = struct('genotype',genotype,'epoch',{'sleep','saline','ketamine'},'state',state);
        Ngroups = numel(groupDefs);
        grpIdx = cell(1,Ngroups);
        grpPow = cell(1,Ngroups);
        grpBoot = cell(1,Ngroups);

        for l=1:Ngroups
            [grpIdx{l},grpPow{l},grpBoot{l}] = getGroupMetrics(groupDefs{l},dataStruct);
        end
        
        % compare each minute of ketamine data to all sleep epoch and all saline epoch
        % repeat with min-by-min saline vs all sleep
        % compare first 5 min or ketamine to first 5 min of saline
        
        % plot min-by-min sleep mean + SEM & ket & sal, shaded error plot
        
