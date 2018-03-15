function [Nepisodes,specVar,specFigures] = getStateSpectralVariance(varargin)
    % [episodeN,specVar,specFigures] = getStateSpectralVariance(eegTraces,eegFS,velTraces,velFS)
    % eegTraces & velTraces should be matrices with each row being a trace from
    % a single animal and single epoch. velocity should be in cm/s and should
    % proably be resampled at a consistent sample rate if camera timing is not
    % exact
    % NAME-VALUE pairs:
    %   - immThresh     : velocity upper threshold for immobility (cm/s) [default=0.4]
    %   - movThresh     : velocity lower threshold for moving (cm/s) [default=0.4]
    %   - sleepDur      : min duration of immobility to classify as sleep (s)  [default=40]
    %   - movInterrupt  : amount of time during a movement episode that vel is allowed to fall under threshold (s) [default 2]
    %   - minEpDur      : min length of any behavioral episode (s) [default=10]
    %   - freqRanges    : cell vector of frequency ranges {[a b],...} to get specgtral variance of [default={[0.5 20],[40,80]}
    %   - winSize       : window size (in seconds) to use when generating spectra [default=2] (default step is winSize/2 for now)
    %   - plotSpec      : flag whether to plot the mean spectra for each state [default=1]
    %   - specFigures   : vector of figures (one per frequency range) on which to plot the spectra. [] (default) results in new figures being created
    %   - figRows       : number of rows of subplots the figure will have [default=1]. Figure will plot each state in a column.
    %   - currFigRow    : row of the figure subplots that the current plots should go on

    immThresh = 0.4;
    movThresh = 0.4;
    sleepDur = 40;
    movInterrupt = 2;
    minEpDur = 10;
    eegTraces = varargin{1};
    eegFs = varargin{2};
    velTraces = varargin{3};
    velFs = varargin{4};
    freqRanges = {[0.5 20],[40 80]};
    Nfreqs = 128;
    winSize = 2;
    winStep =1;
    plotSpec = 1;
    specFigures = [];
    figRows = 1;
    currFigRow = 1;


    assignVars(varargin{5:end})

    if movThresh<immThresh
        movThresh = immThresh;
    end

    behavioralStates = {'sleeping','resting','moving'};
    stateVelRanges = {[0 immThresh],[0 immThresh],[movThresh 80]}; % cm/sec - need to look at states found to check limits
    stateDurationRanges = {[sleepDur 1800],[minEpDur sleepDur],[minEpDur 1800]}; % seconds
    allowedInterrupts = {0,0,movInterrupt}; % seconds outside of velocity thresh allowed in an  episode (to account for breif stops while moving and lags in pos tracking)
    stateParams = struct('state',behavioralStates,'velRange',stateVelRanges,'durRange',stateDurationRanges,'allowedInterrupt',allowedInterrupts);

    eegTime = 0:1/eegFs:(size(eegTraces,2)-1)/eegFs;
    velTime = 0:1/velFs:(size(velTraces,2)-1)/velFs;
    Ntrials = size(eegTraces,1);
    Nstates = numel(behavioralStates);
    Nbands  = numel(freqRanges);
    specVar = zeros(numel(freqRanges),Nstates);

    % Setup new figures for each freqBand if figures are not provided in specFigures
    if plotSpec
        if isempty(specFigures)
            specFigures = gobjects(1,Nbands);
            for l=1:Nbands
                specFigures(l) = figure('Visible','off');
                setFigureProperties(specFigures(l),'Position',[1 1 1900 980])
                subplot(figRows,Nstates,Nstates*(currFigRow-1)+1)
            end
        end
        figI = Nstates*(currFigRow-1);
        if numel(specFigures) ~= Nbands
            error('Not enough figures for bands')
        end
    end

    Nepisodes = zeros(Ntrials,Nstates);
    stateMats = cell(Ntrials,1);
    for k=1:Ntrials
        stateMats{k} = getBehavioralEpisodes(velTime,velTraces(k,:),stateParams);
        Nepisodes(k,:) = arrayfun(@(x) sum(stateMats{k}(:,3)==x),1:Nstates);
    end

    for l=1:numel(freqRanges)
        freqRange = freqRanges{l};
        Spectra = cell(Ntrials,1);
        for k=1:size(eegTraces,1)
            stateMat = stateMats{k};
            [Spectra{k},freqs] = getStateSpectra(eegTraces(k,:),eegTime,eegFs,stateMat,freqRange,winSize); % winStep is automatically 50% of winSize
        end

        specMat = cell2mat(Spectra);
        states = cell2mat(stateMats);
        states = states(:,3);
        if numel(states) ~= size(specMat,1)
            error('Something went wrong, states don''t line up with spectra')
        end

        for k=1:Nstates
            idx = (states==k);
            stateSpec = specMat(idx,:);
            specVar(l,k) = mean(var(stateSpec));
            
            if plotSpec
                figure(specFigures(l))
                subplot(figRows,Nstates,figI+k)
                hold on
                [~,gh] = shadedErrorPlot(freqs,10*log10(mean(specMat)),10*log10(std(specMat)));
                if currFigRow==1
                    title(behavioralStates{k})
                end
                if currFigRow==figRows && k==ceil(Nstates/2)
                    xlabel('Frequency (Hz)')
                end
                if currFigRow==ceil(figRows/2) && k==1
                    ylabel('Power(dB/Hz)')
                end
                legend(gh,sprintf('N=%g',sum(Nepisodes(:,k))))
                xlim(freqRange)
                drawnow;
            end
        end
        if plotSpec
           % suptitle(sprintf('Mean LFP spectra: %g - %g Hz',freqRange))
        end
        Nepisodes = sum(Nepisodes);
    end
