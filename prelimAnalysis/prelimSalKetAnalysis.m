% ketamine saline analysis
if ~exist('dataStruct','var')
    cd ~/Projects/rn_Schizophrenia_Project/RZ-RW_PrelimAnalysisFigs/10s_stats_ShiftedFreq/
    load allPrelimMetrics-10s.mat
end
genotypes = {'WT','Df1'};
state = 'resting';
epochs = {'saline','ketamine'};
animIDs = {'RW2','RW3','RZ2','RZ3'};

testSegment = [5*60 10*60];
Ngeno = numel(genotypes);
bandFreqs = dataParams.band_freqs;
bandNames = dataParams.band_names;
Nbands = numel(bandFreqs);

allFoldChange = cell(1,numel(epochs));
allDat = cell(1,numel(epochs));
allDatLbls = cell(1,numel(epochs));
for e = 1:numel(epochs)
    epoch = epochs{e};
    foldChange = cell(Ngeno,1);
    fcLabels = cell(Ngeno,1);
    for k=1:Ngeno
        geno = genotypes{k};
        groupDef = struct('genotype',geno,'epoch','sleep','behavioral_state',state);
        [grpIdx_sleep,grpPow_sleep,grpBoot_sleep] = getGroupMetrics(groupDef,dataStruct);
        groupDef.segment_time = testSegment;
        groupDef.epoch = epoch;
        [grpIdx,grpPow,grpBoot] = getGroupMetrics(groupDef,dataStruct);
        
        % align sleep periods with ketamine periods and saline periods
        tDat = dataStruct(grpIdx);
        bDat = dataStruct(grpIdx_sleep);
        
        animIdx_baseline = {bDat.animal}';
        animIdx_test = {tDat.animal}';
        animIdx_b = cellfun(@(x) find(strcmp(animIDs,x)),animIdx_baseline);
        animIdx_t = cellfun(@(x) find(strcmp(animIDs,x)),animIdx_test);
        
        tetIdx_baseline = [animIdx_b [bDat.day]' [bDat.tet]'];
        tetIdx_test = [animIdx_t [tDat.day]' [tDat.tet]'];
        tetIdx = unique(tetIdx_test,'rows');
        
        Ntet = size(tetIdx,1);
        bPow = zeros(Ntet,Nbands);
        tPow = zeros(Ntet,Nbands);
        for l=1:Ntet
            x = tetIdx(l,:);
            idx1 = (tetIdx_baseline(:,1)==x(1) & tetIdx_baseline(:,2)==x(2) & tetIdx_baseline(:,3)==x(3));
            idx2 = (tetIdx_test(:,1)==x(1) & tetIdx_test(:,2)==x(2) & tetIdx_test(:,3)==x(3));
            tmpPow1 = cell2mat({bDat(idx1).band_power}');
            tmpPow2 = cell2mat({tDat(idx2).band_power}');
            tmpNorm1 = cell2mat({bDat(idx1).normalization_power}');
            tmpNorm2 = cell2mat({tDat(idx2).normalization_power}');
            
            bPow(l,:) = mean(tmpPow1./tmpNorm1(:,1));
            tPow(l,:) = mean(tmpPow2./tmpNorm2(:,1));
        end
        foldChange{k} = (tPow-bPow)./bPow;
        fcLabels{k} = repmat({[geno '_' epoch]},Ntet,1);
    end
    allDatLbls{e} = vertcat(fcLabels{:});
    allDat{e} = cell2mat(foldChange);
    allFoldChange{e} = foldChange;
end