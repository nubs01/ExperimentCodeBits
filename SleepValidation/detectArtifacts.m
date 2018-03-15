function artifactTimes = detectArtifacts(eegTime,eegDat,Fs,varargin)
    % artfactTimes = detectArtifacts(lfpTime,lfpDat,NAME,VALUE)
    % detects artifacts by:
    %       a) Bandpass filtering the lfp trace at 50-400Hz, and marking
    %          artifacts as those points >5*st.dev. from the mean

    artifactBand = [50 400];
    stdThresh = 5;
    if ~isempty(varargin)
        assignVars(varargin{:})
    end

    [b,a] = butter(10,artifactBand/(Fs/2),'bandpass');
    eegDat = interpInf(eegDat); 
    filtDat = filtfilt(b,a,eegDat);
    fMean = mean(filtDat);
    fSTD = std(filtDat);
    artThresh = stdThresh*fSTD;
    artifactTimes = eegTime(filtDat>(fMean+artThresh) | filtDat<(fMean-artThresh));
