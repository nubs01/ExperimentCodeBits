function allVel = getAnimVel(dataDir,animID,epochs,varargin)
    % allVel = getAnimVel(dataDir,animID,epochs,varargin)
    % returns all velocity traces for the given anim & epochs ([session epoch])
    % this is a cell array of velTraces with each element corrsponding to a row of epochs.
    % If truncateDat is set to nan (default) then you get the full data, if
    % truncateDat is set to a number, then the data is truncated to that many
    % minutes and returned as a matrix with a velTrace on each row.
    % each velocity trace is passed through cleanVelocity to low-pass at 5Hz and resample at velFs Hz
    % The raw trace is assumed to at about 30 Hz.
    % NAME-VALUE pairs:
    %       - truncateDat : minutes to truncate data to, nan (default) does not truncate, and -1 truncates to match 
    %       - velFs       : desired sampling rate of velocity
    %       - forceTrunc  : whether to force truncation length. if 1, epochs shorter than truncateDat will be padded with -1, otherwise all data will be shortened to the minimum length

    truncateDat = nan;
    velFs = 25;
    forceTrunc = 0;

    assignVars(varargin)

    fName = @(x) sprintf('%s%spos%02i.mat',dataDir,animID,x);

    allVel = cell(size(epochs,1),1);
    truncateMin = inf;

    for k=1:size(epochs,1)
        fn = fName(epochs(k,1));
        dat1 = load(fn);
        dat = dat1.pos;
        dat = dat{epochs(k,1)}{epochs(k,2)};
        tmp1 = strsplit(dat.fields);
        r1 = strcmp(tmp1,'time');
        r2 = strcmp(tmp1,'vel');
        retDat = cleanVelocity(dat.data(:,r2),dat.data(:,r1),'newFs',velFs);
        truncateMin = min([numel(retDat) truncateMin]);
        if ~isrow(retDat)
            retDat = retDat';
        end
        allVel{k} = retDat;
    end

    if isnan(truncateDat)
        return;
    end

    if truncateDat>0 && ~forceTrunc
        truncateTo = min([truncateMin velFs*truncateDat*60]);
    elseif truncateDat>0
        truncateTo = velFs*truncateDat*60;
    end
    for k=1:numel(allVel)
        if numel(allVel{k})>=truncateTo
            allVel{k} = allVel{k}(:,1:truncateTo);
        else
            allVel{k} = paddarray(allVel{k},[0 truncateTo-numel(allVel{k})],-1,'post');
        end
    end
    allVel = cell2mat(allVel);

