function [outIdx,normPow,bootDat] = getGroupMetrics(groupDef,dataStruct,varargin)
    % [outIdx,normPow,bootDat] = gretgroupMetrics(groupDef,dataStruct,NAME,VALUE)
    % groupDef should be a structure with the fields that you want to search
    % for example if groupDef = struct('genotype','Df1','behavioral_state','sleeping') will return data for all Df1 sleeping episodes
    % the segment_time field will be processed slightly differently, you provide a range and any segment times that fall in that range will returned
    % NAME-VALUE
    %   - nboot : number of times to bootstrap [defualt=1000]
    

    nboot = 1000;
    assignVars(varargin);

    fn = fieldnames(groupDef);
    outIdx = [];

    % get indices based on groupDef
    for k=1:numel(fn)
        if strcmp(fn{k},'epoch') && strcmp(groupDef.(fn{k}),'all')
            continue;
        end
        if strcmp(fn{k},'segment_time')
            segR = groupDef.(fn{k});
            segT = [dataStruct.(fn{k})];
            idx = (segT(1,:)>=segR(1) & segT(2,:)<=segR(2));
        else
            idx = strcmp({dataStruct.(fn{k})},groupDef.(fn{k}));
        end
        if isempty(outIdx)
            outIdx = idx;
        else
            outIdx = (outIdx & idx);
        end
    end

    % get normalized power
    dat = dataStruct(outIdx);
    bandPow = cell2mat({dat.band_power}');
    tmpNorm = cell2mat({dat.normalization_power}');
    normPow = bandPow./tmpNorm(:,1);
    bootDat = bootstrapDat(normPow,nboot);
