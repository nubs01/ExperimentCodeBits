function outVar = geAverageVariance(varargin)
% outVar = getAverageVariance(spectra,stateMat) with spectra having a row
% for each spectrum, and their being one spectrum per state (row) in stateMat

spectra = varargin{1};
stateMat = varargin{2};
