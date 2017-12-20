function fh = plotVelocityHistogram(vel,states,binEdges,varargin)
    % fh = plotVelocityHistogram(vel,states,binEdges) - binEdges is optional, default is 0:.2:5
    % fh = plotVelocityHistogrtam(fh,vel,states,binEdges)
    if ishandle(vel)
        fh = vel;
        vel = states;
        states = binEdges;
        if ~isempty(varargin)
            binEdges = varargin{1};
        end
    else
        fh = figure();
    end
    setFigureProperties(fh);
    if isrow(vel)
        vel = vel';
    end
    if ~exist('binEdges','var')
        binEdges = 0:.2:max([ceil(vel);5]);
    end

    allStates = unique(states);
    hold on
    hh = nan(1,numel(allStates));
    for k = 1:numel(allStates)
        hh(k) = histogram(vel(strcmp(states,allStates{k})),binEdges);
    end
    legend(hh,allStates)
    xlabel('Velocity (cm/sec)')
    ylabel('# of episodes')

