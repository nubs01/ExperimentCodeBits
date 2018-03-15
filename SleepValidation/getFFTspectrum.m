function [S,f] =  getFFTspectrum(x,fs,varargin)
%OBSOLETE
    if isrow(x)
        x = x';
    end
    nfft = size(x,1);
    w = hann(nfft);
    x2 = x.*w;
    X = fft(x2,nfft);
    mx = abs(X).^2;
    mx = mx/(w'*w);
    npts = nfft/2+1;
    mx = mx(1:npts,:);
    mx(2:end-1,:) = mx(2:end-1,:).*2;
    S = mx./fs;
    f = (0:npts-1)*fs/nfft;
