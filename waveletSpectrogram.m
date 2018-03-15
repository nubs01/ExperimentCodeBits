function [fh,S,t,f,COI] = waveletSpectrogram(eeg,time,fs,freqs,zFlag,fh)
% makes a spectrogram of eeg data using the continuous wavelet transform. eeg
% should be an array containing an eeg trace
if ~isrow(time)
    time = time';
end
if ~isrow(eeg)
    eeg = eeg';
end

if ~exist('freqs','var')
    freqs = 0.5:.25:50;
elseif ~isrow(freqs)
    freqs = freqs';
end
if ~exist('zFlag','var')
    zFlag = 0;
end
if numel(freqs)==2 % if frequency is range 
    freqs = freqs(1):.25:freqs(2);
elseif numel(freqs)==3 % if freqs is range with step [start step end]
    freqs = freqs(1):freqs(2):freqs(3);
end

% divide eeg trace into equal segments less than 90000 points each and run cwt
maxPts = 2700000;
nSeg = ceil(numel(eeg)/maxPts);
segPts = fix(numel(eeg)/nSeg);

cwtIdx = @(i) (i-1)*segPts+1:i*segPts;
S = [];
f = [];
t = [];
COI = [];
for i=1:nSeg
    idx = cwtIdx(i);
    [WT,FF,C] = cwt(eeg(idx),fs);
    FF = FF(end:-1:1);
    WT = WT(end:-1:1,:);
    Si = interp1(FF,WT,freqs);
    Sic  = abs(Si.*conj(Si));
    if isempty(S)
        S = Sic;
        f = freqs;
        t = time(idx);
        COI = C';
    else
        S = [S Sic];
        t = [t time(idx)];
        COI = [COI C'];
    end
    clear WT FF C Si Sic idx
end
if zFlag
    S = zscore(S,0,2);
    cbLabel = 'Power Z-Score';
else
    cbLabel = 'Power';
end

% Downsample data to reduce memory usage and variable size
S = smoothdata(S,2,'gaussian',fs/4); % smooth with .25 sec gaussian window
S = downsample(S',fs/4)'; % Then downsample taking every fs/4th point
t = downsample(t,fs/4); % same with time


if ~exist('fh','var')
    fh = figure('Position',[680 338 810 615]);
elseif fh~=0
    figure(fh);
    hold on
end
if fh~=0
    imagesc(t,f,S);
    set(gca,'ydir','normal')
    cb = colorbar;
    colormap('jet')
    ylabel('Frequency (Hz)')
    xlabel('Time')
    ylabel(cb,cbLabel)
    title('Wavelet Spectrogram')
    mS = mean(mean(S));
    stdS = std(std(S));
    caxis([mS-3*stdS mS+3*stdS])
end
