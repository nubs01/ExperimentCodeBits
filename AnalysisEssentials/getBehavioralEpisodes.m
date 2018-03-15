function stateMat = getBehavioralEpisodes(posTime,vel,stateParams)
    % stateMat = getBehavioralEpidsodes(posTime,vel,stateParams) returns a matrix
    % with the start time, end time, and behavioral state index as columns.  
    % stateParams is a structure array where each element is a definition of a
    % behavioral state with fields: state,velRange,and durRange and allowedInterrupt
    
    stateMat = [];
    for k=1:numel(stateParams)
        vLim = stateParams(k).velRange;
        dLim = stateParams(k).durRange;
        sIdx = (vel>=vLim(1) & vel<=vLim(2));
        if ~any(sIdx==1)
            tmpTimes = [];
        else
            sIdx = contiguous(sIdx,1);
            sIdx = sIdx{2};
            tmpTimes = posTime(sIdx);
        end

        % merge episodes if time in-between is less than allowed interrupt
        if ~isempty(tmpTimes) && size(tmpTimes,1)>1
            tmpTimes = mergeEpisodes(tmpTimes,stateParams(k).allowedInterrupt);
        end
        stateDur = diff(tmpTimes,1,2);
        sIdx = (stateDur>=dLim(1) & stateDur<=dLim(2));
        tmpTimes = tmpTimes(sIdx,:);
        stateMat = [stateMat;[tmpTimes repmat(k,size(tmpTimes,1),1)]];
    end
        

