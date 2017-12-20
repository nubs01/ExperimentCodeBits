function fh = plotAllStateSpectra(varargin)
    % fh = plotAllStatesSpectra(spectra,freq,states)
    % fh = plotAllStateSpectra(fh,spectra,freq,states)
    % fh = plotAllStateSpectra(...,sortOrder)

    if ishandle(varargin{1})
        fh = varargin{1};
        varargin = varargin(2:end);
    else
        fh = figure();
    end
    setFigureProperties(fh);

    spectra = varargin{1};
    freq = varargin{2};
    states = varargin{3};
    if numel(varargin)==4
        sortOrder = varargin{4};
    end

    if ~exist('sortOrder','var')
        [states,sortIdx] = sort(states);
    else
        sortIdx = [];
        for k=1:numel(sortOrder)
            sortIdx = [sortIdx find(strcmp(states,sortOrder{k}))];
            if ~isrow(sortIdx)
                sortIdx = sortIdx';
            end
        end
        states = states(sortIdx);
    end

    if iscell(spectra)
        spectra = spectra(sortdx);
        specDat = spectra;
        spectra = zeros(numel(specDat),numel(freqs));
        for k=1:numel(specDat)
            spectra(k,:) = specDat{k};
        end 
    else
        spectra = spectra(sortIdx,:);
    end
    allStates = unique(states);

    % Z-score spectra along each frequency
    %spectra = zscore(spectra);

    imagesc(freq,1:size(spectra,1),spectra)
    set(gca,'YDir','normal')
    cb = colorbar;
    colormap('jet');
    % yticks = get(gca,'YTick');
    % yticklabels = get(gca,'YTickLabel');
    yticks = [];
    yticklabels = {};
    for k=1:numel(allStates)
        sIdx = strcmp(states,allStates{k});
        if ~any(sIdx)
            disp(['No ' allStates{k} ' episodes found.'])
            continue;
        end
        if sIdx(end)
            breakPoint = find(sIdx,1,'first')-0.5;
        else
            breakPoint = find(sIdx,1,'last')+0.5;
        end
        hold on
        plot([freq(1) freq(end)],[breakPoint breakPoint],'k-','LineWidth',3);
        midPt = median(find(sIdx));
        yticks = [yticks midPt];
        yticklabels = [yticklabels allStates{k}];
    end
    xlabel('Frequency (Hz)')
    ylabel('Behavioral Episode')
    [yticks,sIdx] = sort(yticks);
    yticklabels = yticklabels(sIdx);
    set(gca,'YTick',yticks,'YTickLabel',yticklabels,'YTickLabelRotation',90)
    ylabel(cb,'Power')


  
