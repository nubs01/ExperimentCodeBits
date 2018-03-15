function [spectra,freqs] = getStateSpectra(eegDat,eegTime,fs,stateMat,freqRange,winSize)
    % [spectra,freqs] = getStateSpectra(eegDat,eegTime,stateMat,states,freqRange,winSize)
    % takes the eegDat (raw lfp trace) and for each behavioral episode, breaks
    % it into winSize second epochs and averages fft spectra to get an average
    % spectrum for eac h behavioral episode. uses chronux function
    % mtspectrumsegc. Spectra are then interpolated to get 256 frequency points
    % covering freqRange. spectra is an episodes x frequency array with rows
    % corresponding to the spectrum for the corresponding row of stateMat. 
    

    nFreqs = 128;
    spectra = zeros(size(stateMat,1),nFreqs);
    freqs = linspace(freqRange(1),freqRange(2),nFreqs);
    for k=1:size(stateMat,1)
        t1 = stateMat(k,1);
        t2 = stateMat(k,2);
        epSeg = eegDat(eegTime>=t1 & eegTime<=t2);
        [S,freqs] = getPSD(epSeg,fs,'freqRange',freqRange,'winSize',winSize,'freqRes',mean(diff(freqs)),'noverlap',winSize/2,'psdType','stft'); 
        % TODO: can try to use C to whiten matrix
        % TODO: how to average confidence interval/errors, then I can use Serr
        spectra(k,:) = S;
        clear S f t1 t2 epSeg varS C Serr S1
    end
    

