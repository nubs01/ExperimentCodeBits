grpIdx = cell(1,4);
grpSpec = cell(1,4);
grpSpecSEM = cell(1,4);
for k=1:4
    groupDef = struct('animal',animIDs{k},'epoch','sleep','behavioral_state','sleeping');
    grpIdx{k} = getGroupMetrics(groupDef,dataStruct);
    grpDat = dataStruct(grpIdx{k});
    grpS = {grpDat.spectrum}';
    grpS = cell2mat(grpS);
    grpNorm = {grpDat.normalization_power}';
    grpNorm = cell2mat(grpNorm);
    grpS = grpS./grpNorm(:,1);
    grpSpecBoot = bootstrapDat(grpS,1000);
    grpSpec{k} = [grpSpecBoot.mean];
    grpSpecSEM{k} = [grpSpecBoot.SEM];
end
