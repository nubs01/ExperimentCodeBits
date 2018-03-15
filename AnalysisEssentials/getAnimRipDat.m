function ripDat = getAnimRipDat(dataDir,animID,tetIdx,varargin)
    ripDat = [];
    for k=1:size(tetIdx,1)
        ripples = load(sprintf('%s%sripples%02i.mat',[dataDir filesep],animID,tetIdx(k,1)));
        rD = ripples.ripples{tetIdx(k,1)}{tetIdx(k,2)}{tetIdx(k,3)};
        if isempty(ripDat)
            ripDat = rD;
        else
            ripDat(k) = rD;
        end
        clear ripples
    end
