function allDat = getAnimEEG(dataDir,animID,tetIdx,varargin)
    % allDat = getAnimEEG(dataDir,animID,epochs,varargin)
    % returns all raw eeg traces for the given anim & tetIdx ([session epoch tet])
    % this is a cell array of eegTraces with each element corrsponding to a row of epochs.
    % If truncateDat is set to nan (default) then you get the full data, if
    % truncateDat is set to a number, then the data is truncated to that many
    % minutes and returned as a matrix with an eegTrace on each row. It is assumed eeg data is in dataDir/EEG/
    % NAME-VALUE pairs:
    %       - truncateDat : minutes to truncate data to, nan (default) does not truncate, and -1 truncates to match 
    %                       smallest data set
    %
    %       - dataType    : string of which data structure to get, 'eeg'
    %                       (defualt) can also be changed to any structure in the EEG folder
    %                       such as 'theta' or 'eegref'
    %
    %       - returnField : field of data structure to return
    %
    %       - returnCol   : column of data to return is the  field contains
    %                       more than one column, set to nan (default) to return  all columns,
    %                       multi-column data will only be retruned in a cell array
    
    truncateDat = nan;
    dataType = 'eeg';
    eegDir = [dataDir filesep 'EEG' filesep];
    returnField = 'data';
    returnCol = nan;

    assignVars(varargin)

    fName = @(x) sprintf('%s%s%s%02i-%02i-%02i.mat',eegDir,animID,dataType,x);

    allDat = cell(size(tetIdx,1),1);
    fsDat = zeros(size(tetIdx,1),1);
    truncateTo = inf;

    for k=1:size(tetIdx,1)
        fn = fName(tetIdx(k,:));
        dat = load(fn);
        dat = dat.(dataType);
        dat = dat{tetIdx(k,1)}{tetIdx(k,2)}{tetIdx(k,3)};
        fsDat(k) = dat.samprate;
        retDat = dat.(returnField);
        if all(size(retDat)>1) && ~isnan(returnCol)
            retDat = retDat(:,returnCol);
            truncateTo = min([numel(retDat),truncateTo]);
        elseif all(size(retDat)>1) && isnan(returnCol)
            retDat = retDat';
            truncateTo = min([size(retDat,2),truncateTo]);
        else
            truncateTo = min([numel(retDat) inf]);
        end
        if any(size(retDat)==1) && ~isrow(retDat)
            retDat = retDat';
        end
        allDat{k} = retDat;
    end

    if isnan(truncateDat)
        return;
    end

    if truncateDat>0
        if numel(unique(fsDat))==1
            truncateTo = min([truncateTo unique(fsDat)*truncateDat*60]);
        else
            truncateTo = fsDat*truncateDat*60;
        end
    end
    if numel(truncateTo)==1
        truncateTo = truncateTo*ones(numel(allDat),1);
    end
    toFlatten = 1;
    for k=1:numel(allDat)
        allDat{k} = allDat{k}(:,1:truncateTo(k));
        if size(allDat{k},1)>1
            toFlatten=0;
        end
    end
    if toFlatten
        allDat = cell2mat(allDat);
    end

