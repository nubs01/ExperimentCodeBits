function fh = plotDurationHistogram(stateTimes,states,binEdges,varargin)
    % fh = plotDuartionHistogram(stateTimes,states,binEdges) - binEdges is optional, default is 0:2:100
    % fh = plotDuartionHistogrtam(fh,stateTimes,states,binEdges)
    % stateTimes is an nx2 array with start and end times of each behavioral episode
    if ishandle(stateTimes)
        fh = stateTimes;
        stateTimes = states;
        states = binEdges;
        if ~isempty(varargin)
            binEdges = varargin{1};
        end
    else
        fh = figure();
    end
    setFigureProperties(fh);
    stateDur = diff(stateTimes,1,2);
    if ~exist('binEdges','var')
        binEdges = 0:2:max([ceil(stateDur);100]);
    end

    allStates = unique(states);
    hold on
    hh = nan(1,numel(allStates));
    for k = 1:numel(allStates)
        hh(k) = histogram(stateDur(strcmp(states,allStates{k})),binEdges);
    end
    legend(hh,allStates)
    xlabel('Episode Duration (sec)')
    ylabel('# of episodes')

