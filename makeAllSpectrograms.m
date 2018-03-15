function makeAllSpectrograms(animID,dataDir,days,figDir,epochNames)
    diary([figDir filesep animID 'makeSpectrogram.log'])
    set(0,'defaultaxesfontsize',14)
    set(0,'defaultaxesfontname','Arial')
    figSize =[680 228 810 615]; 

    % varaibles
    maxTrace = 30*60; % time to truncate eeg traces to
    stftWin = 5; % sec to window for chronux spectrogram
    stftStep = 2.5; % sec
    stftFpass = [0 50];
    stftTapers = [3 4];

    cwtFrange = [0.5 .1 40];

    deltaBand = [0.5 5];
    thetaBand = [6 10];

    if dataDir(end)==filesep
        dataDir = dataDir(1:end-1);
    end
    if figDir(end)~=filesep
        figDir = [figDir filesep];
    end
    if ~isrow(days)
        days = days';
    end
    figBaseDir = figDir;
    tetinfo = load([dataDir filesep animID 'tetinfo.mat']);
    tetinfo = tetinfo.tetinfo;
    tetDescrip = cellfetch(tetinfo,'descrip');
    a = strcmp(tetDescrip.values,'riptet');
    ripIdx = tetDescrip.index(a,:);
    %ripIdx = ripIdx(arrayfun(@(x) any(days==x),ripIdx(:,1)),:);
    
    eegDir = [dataDir filesep 'EEG' filesep];
    for d=days
        tic
        fprintf('Analyzing Day %02i Data...\n',d);
        figDir = [figBaseDir sprintf('%02i_SpectralAnalysis',d) filesep];
        specDatDir = [dataDir filesep 'SpectralData' filesep];
        if ~exist(specDatDir,'dir')
            mkdir(specDatDir);
        end
        if ~exist(figDir,'dir')
            mkdir(figDir);
        end
        idx = ripIdx(ripIdx(:,1)==d,:);
        epochs = sort(unique(idx(:,2)));
        tets = sort(unique(idx(:,3)));
        eN = epochNames{d};
        sleepEpoch = find(strcmp(lower(eN),'sleep'));
        for tt=tets'
            fprintf('    Processing tetrode %02i...\n',tt)
            stftSpec = [];
            stftTime = [];
            timeBreak = [];
            timeTicks = [];
            timeTickLabels = {};
            stftFreq = [];

            cwtTimeBreak = [];
            cwtTicks = [];
            cwtTickLabels = {};
            cwtSpec = [];
            cwtTime = [];
            cwtFreq = [];
            cwtCOI = [];

            sleepEpochTime = [];

            saveFiles = {'stftSpectrogram','stftDeltaTrace','stftDeltaTraceZ','cwtSpectrogram','cwtDeltaTrace','stftDeltaTrace-ThetaNorm','cwtDeltaTraceZ','cwtDeltaTrace-ThetaNorm'}';
            saveFiles = strcat(animID,saveFiles);
            saveFiles = strcat(figDir,saveFiles);
            saveFiles = strcat(saveFiles,sprintf('%02i-%02i',d,tt));
            saveFilename = @(s) sprintf([figDir filesep animID s '%02i-%02i'],d,tt); 

            stftSpecFile =sprintf('%s%s%sstftSpecData%02i-%02i.mat',specDatDir,filesep,animID,d,tt);
            cwtSpecFile = sprintf('%s%s%scwtSpecData%02i-%02i.mat',specDatDir,filesep,animID,d,tt);  

            if ~exist(stftSpecFile,'file') || ~exist(cwtSpecFile,'file')
                for ep = epochs'
                    fprintf('        Epoch %02i...\n',ep)
                    eeg = load(sprintf('%s%seeg%02i-%02i-%02i.mat',eegDir,animID,d,ep,tt));
                    eeg = eeg.eeg;
                    eeg = eeg{d}{ep}{tt};
                    
                    fs = eeg.samprate;
                    t = eeg.starttime:1/fs:eeg.endtime;
                    e = eeg.data;
                    if ~isrow(e)
                        e = e';
                    end

                    % get rid of pinging noise by filtering with 200 Hz low-pass
                    fprintf('         - Low-pass filtering eeg @ 200Hz\n')
                    [b,a] = butter(10,200/(fs/2));
                    % interpolate infinite and nan points
                    infIdx = ~isfinite(e);
                    if any(infIdx)
                        fprintf('         - Interpolating over %i non-finite points in eeg...\n',sum(infIdx));
                        T = 1:numel(e);
                        e(infIdx) = interp1(T(~infIdx),e(~infIdx),T(infIdx));
                    end
                    ef = filtfilt(b,a,e);
                    clear infIdx b a
                    % truncate to maxTrace
                    if numel(e)>maxTrace*fs
                        fprintf('         - Truncate eeg to 30 min\n')
                        e = e(1:maxTrace*fs);
                        ef = ef(1:maxTrace*fs);
                        t = t(1:maxTrace*fs);
                    else
                        % cut down to nearest minute if < maxTrace
                        maxMin = fix((t(end)-t(1))/60);
                        fprintf('         - EEG data for d-e-t %02i-%02i-%02i is less than maxTrace. Truncating to %d minutes...\n',d,ep,tt,maxMin);
                        e = e(1:maxMin*60*fs);
                        t = t(1:maxMin*60*fs);
                        ef = ef(1:maxMin*60*fs);
                    end

                    if ep==sleepEpoch
                        sleepEpochTime = [sleepEpochTime;t(1) t(end)]/60;
                    end

                    % Short-time fourier transform analysis (chronux)
                    fprintf('         - Calculate Chronux STFT Spectrogram\n')
                    params = struct('Fs',fs,'err',[2 0.05],'fpass',stftFpass,'tapers',stftTapers);
                    [specS,specT,specF,specE] = mtspecgramc(ef,[stftWin stftStep],params);
                    specT = specT + t(1); % add time offset
                    specT = specT/60; % convert time from sec to minutes for plotting
                    if isempty(stftSpec)
                        stftSpec = specS;
                        stftTime = specT;
                        stftFreq = specF;
                        timeTicks = [specT(1)  mean(specT) specT(end)];
                        timeTickLabels = {sprintf('^{%0.0f}',specT(1)) eN{ep} sprintf('^{%0.0f}',specT(end))};
                    else
                        stftSpec = [stftSpec;specS];
                        timeBreak = [timeBreak;stftTime(end) specT(1)];
                        stftTime = [stftTime specT];
                        timeTicks = [timeTicks specT(1) mean(specT) specT(end)];
                        timeTickLabels = [timeTickLabels sprintf('^{%0.0f}',specT(1)) eN{ep} sprintf('^{%0.0f}',specT(end))];
                    end
                    clear specS specT specF specE params

                    % Continuous wavelet transform analysis
                    fprintf('         - Calculate Wavelet Spectrogram (cwt)\n')
                    [~,specS,specT,specF,specCOI] = waveletSpectrogram(ef,t/60,fs,cwtFrange,0,0);
                    if ~isempty(cwtSpec)
                        cwtSpec = [cwtSpec specS];
                        cwtTimeBreak = [cwtTimeBreak;cwtTime(end) specT(1)];
                        cwtTime = [cwtTime specT];
                        cwtCOI = [cwtCOI specCOI];
                        cwtTicks = [cwtTicks specT(1) mean(specT) specT(end)];
                        cwtTickLabels = [cwtTickLabels sprintf('^{%0.0f}',specT(1)) eN{ep} sprintf('^{%0.0f}',specT(end))];
                    else
                        cwtSpec = specS;
                        cwtTime = specT;
                        cwtCOI = specCOI;
                        cwtFreq = specF;
                        cwtTicks = [specT(1)  mean(specT) specT(end)];
                        cwtTickLabels = {sprintf('^{%0.0f}',specT(1)) eN{ep} sprintf('^{%0.0f}',specT(end))};
                    end

                    clear e t ef fs eeg
                end
                % save spectrogram data
                fprintf('    - Saving spectra for tetrode %02i\n',tt)
                save(stftSpecFile,'stftSpec','stftTime','stftFreq','timeBreak','timeTicks','timeTickLabels','sleepEpochTime')
                save(cwtSpecFile,'cwtSpec','cwtTime','cwtFreq','cwtCOI','cwtTimeBreak','cwtTicks','cwtTickLabels')

            else
                fprintf('    - Spectrogram data found for tetrode %02i. Loading...\n',tt)
                load(stftSpecFile)
                load(cwtSpecFile)
            end
            
            if ~exist('sleepEpochTime','var') || isempty(sleepEpochTime)
                fprintf('     - Cannot find sleepEpochTime. Setting to first epoch...\n')
                sleepEpochTime = [stftTime(1) timeBreak(1,1)];
                save(stftSpecFile,'sleepEpochTime','-append')
            end

            % Plot stft Spectrogram
            if ~exist([saveFiles{1} '.svg'],'file')
                fprintf('    - Plotting STFT spectrogram\n')
                T2 = stftTime(1):stftTime(2)-stftTime(1):stftTime(end);
                S2 = interp1(stftTime,stftSpec,T2);

                stftFig1 =  figure('Position',figSize,'Visible','off');
                freqTrim = (stftFreq>=0.5 & stftFreq<=40);
                plot_matrix(S2(:,freqTrim),T2,stftFreq(freqTrim))
                hold on
                for tb = 1:size(timeBreak,1)
                    wid = timeBreak(tb,2)-timeBreak(tb,1);
                    ylim = get(gca,'ylim');
                    hei = diff(ylim);
                    rectangle('Position',[timeBreak(tb,1) ylim(1) wid hei],'FaceColor',[1 1 1],'EdgeColor','none')
                end
                colormap('jet')
                title({[animID ' STFT Spectrogram'],sprintf('Day %02i, Tet %02i',d,tt)})
                xlabel('Time (min)')
                ylabel('Frequency (Hz)')
                cb = colorbar;
                ylabel(cb,'Power')
                set(gca,'XTick',timeTicks,'XTickLabels',timeTickLabels)
                saveas(stftFig1,saveFiles{1},'svg')
                close(stftFig1)
                fprintf('    - STFT Spectrogram saved as %s\n',saveFiles{1}(find(saveFiles{1}==filesep,1,'last')+1:end))
            end
            
            % Analyze stft delta power
            % So I'm trying to see if
            % (delta_ket-delta_sleep) > (delta_sal-delta_sleep) -- for both df1 & wt,
            % (delta_ket-delta_sal)_df1 > (delta_ket-delta_sal)_wt,
            % (theta_sleep)_df1 == (theta_sleep)_wt and
            % (delta_sleep/thea_sleep)_df1 > (delta_sleep/theta_sleep)_wt
            % Raw
            if ~exist([saveFilename('stftDeltaTrace-raw') '.svg'],'file')
                fprintf('    - Plotting STFT Traces\n')

                % Plot Raw delta trace
                [deltaPow,deltaSTD,dbMean,dbSTD] = getBandPower(stftSpec,stftTime,stftFreq,deltaBand,0,sleepEpochTime(1,:)); 
                deltaFig1 = makeTracePlot(stftTime,smoothdata(deltaPow,'gaussian',24),smoothdata(deltaSTD,'gaussian',24),...
                                          timeBreak,timeTicks,timeTickLabels,...
                                          {[animID ' STFT Delta Power'],sprintf('Day %02i, Tet %02i',d,tt)},...
                                          'Time (min)','Power','b'); 
                saveas(deltaFig1,saveFilename('stftDeltaTrace-raw'),'svg')
                close(deltaFig1)
            


                % Raw theta trace
                [thetaPow,thetaSTD,tbMean,tbSTD] = getBandPower(stftSpec,stftTime,stftFreq,thetaBand,0,sleepEpochTime(1,:));
                thetaFig1 = makeTracePlot(stftTime,smoothdata(thetaPow,'gaussian',24),smoothdata(thetaSTD,'gaussian',24),...
                                          timeBreak,timeTicks,timeTickLabels,...
                                          {[animID ' STFT Theta Power'],sprintf('Day %02i, Tet %02i',d,tt)},...
                                          'Time (min)','Power','b');
                saveas(thetaFig1,saveFilename('stftThetaTrace-raw'),'svg')
                close(thetaFig1)

                % Theta normalized delta trace
                tnormDeltaPow = deltaPow./tbMean;
                tnormDeltaSTD = abs(tnormDeltaPow).*sqrt((deltaSTD./deltaPow).^2+(tbSTD/tbMean)^2);
                deltaFig2 = makeTracePlot(stftTime,smoothdata(tnormDeltaPow,'gaussian',24),smoothdata(tnormDeltaSTD,'gaussian',24),...
                                          timeBreak,timeTicks,timeTickLabels,...
                                          {[animID ' STFT Delta Power'],sprintf('Day %02i, Tet %02i',d,tt),'^{^{normalized to baseline theta power}}'},...
                                          'Time (min)','Normalized Power','b');
                saveas(deltaFig2,saveFilename('stftDeltaTrace-ThetaNormalized'),'svg')
                close(deltaFig2)

                % Z-scored to sleep mean & std Delta trace
                [deltaZPow,deltaZSTD,dbMean,dbSTD] = getBandPower(stftSpec,stftTime,stftFreq,deltaBand,1,sleepEpochTime(1,:));
                deltaFig3 = makeTracePlot(stftTime,smoothdata(deltaZPow,'gaussian',24),smoothdata(deltaZSTD,'gaussian',24),...
                                          timeBreak,timeTicks,timeTickLabels,...
                                          {[animID ' STFT Delta Power'],sprintf('Day %02i, Tet %02i',d,tt),'^{^{Z-scored to baseline epoch}}'},...
                                          'Time (min)','Z-Scored Power','b');
                saveas(deltaFig3,saveFilename('stftDeltaTrace-Z-Scored'),'svg')
                close(deltaFig3)

            end
            clear stftTime T2 stftFreq stftSpec freqTrim timeBreak timeTicks timeTickLabels

            % Plot CWT spectrogram
            if ~exist([saveFiles{4} '.svg'],'file')
                fprintf('    - Plotting CWT spectrogram\n')
                T2 = cwtTime(1):cwtTime(2)-cwtTime(1):cwtTime(end);
                cwtSpec = interp1(cwtTime,cwtSpec',T2)';
                cwtFig1 = figure('Position',figSize,'Visible','off');
                imagesc(T2,cwtFreq,cwtSpec);
                set(gca,'ydir','normal')
                for tb = 1:size(cwtTimeBreak,1)
                    wid = cwtTimeBreak(tb,2)-cwtTimeBreak(tb,1);
                    ylim = get(gca,'ylim');
                    hei = diff(ylim);
                    rectangle('Position',[cwtTimeBreak(tb,1) ylim(1) wid hei],'FaceColor',[1 1 1],'EdgeColor','none')
                end
                cb = colorbar;
                colormap('jet')
                ylabel('Frequency (Hz)')
                ylabel(cb,'Power')
                title({[animID ' CWT Spectrogram'],sprintf('Day %02i, Tet %02i',d,tt)})
                mS = mean(mean(cwtSpec));
                stdS = std(std(cwtSpec));
                caxis([0 mS+3*stdS])
                set(gca,'XTick',cwtTicks,'XTickLabels',cwtTickLabels)
                xlabel('Time (min)')
                saveas(cwtFig1,saveFiles{4},'svg')
                close(cwtFig1)
                fprintf('    - CWT Spectrogram saved as %s\n',saveFiles{4}(find(saveFiles{4}==filesep,1,'last')+1:end))
            end
            clear cwtTicks cwtTickLabels cwtTimeBreak cwtSpec cwtTime cwtCOI cwtFreq cwtFig1 stftFig1
        end
        toc
        clear epochs tets idx eN
    end
    diary off
