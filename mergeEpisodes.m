function out = mergeEpisodes(epTimes,allowedGap)
    
    gaps = [epTimes(2:end,1);inf] - epTimes(:,2);
    B = gaps<=allowedGap;
    if sum(B)==0
        out = epTimes;
        return;
    end
    C = contiguous(B,1);
    C = C{2};
    out = epTimes;
    out(C(:,1),2) = out(C(:,2)+1,2);
    out(C(:,2)+1,:) = [];

