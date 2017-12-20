function immTimes = getImmobileTimes(pos,varargin)
    % returns an array with start and end times of immobile periods. The default
    % settings for an immobile period are a minimum duration of 1 second and a
    % max speed of 4 px/sec. To change these settings simply pass the arguments
    % 'thresh',# or 'minDur',#. This assumes pos is the unpacked pos data
    % structure with fields: units, cmperpixel, data (5 columns: time x y dir
    % vel). Can also pas filtered velocity data ('velDat')  and time data
    % (timeDat) as long as the new vectors are the same units as the previous
    % vectors.

    thresh = 4; % px/sec
    minDur = 1; % sec
    timeDat = pos.data(:,1);
    velDat = pos.data(:,5);
    if ~isempty(varargin)
        assignVars(varargin);
    end
    immTimes=[];
    if strcmp(pos.units,'cm')
        velDat = velDat/pos.cmperpixel;
    end
    immDat = velDat<=thresh;
    immIdx = contiguous(immDat,1);
    immIdx = immIdx{2};
    immTimes = timeDat(immIdx);
    diffTime = diff(immTimes,1,2);
    keepIdx = find(diffTime>=minDur);
    immTimes = immTimes(keepIdx,:);


