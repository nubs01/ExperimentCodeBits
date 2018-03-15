% metadata: animIDs, free params, 'fixed' params, animDirs, epochs to test, figDirs
%
% get day, epoch and tetrode indices for each animal
% make eeg and velocity matrices (trim eeg, clean-up velocity)
% maintain an index array [animal day epoch tet]
%
% setup behavioralState parameters & parameters to vary
% run getStateSpectralVariance to generate plots & collect variances
% for now, save data
%
% TODO: plot surfaces and see if variances can be minimized on all states
%       without sacrificing Nepisodes greatly

animIDs = {'RW2','RW3'};
projDir = '~/Projects/rn_Schizophrenia_Project/';
expDirs = strcat(animIDs,'_WT_Experiment/');
dataDirs = strcat(projDir,strcat(expDirs,strcat(animIDs,'_direct/')));
figDir = [projDir 'RZ-RW_PrelimAnalysisFigs/ImmThresh-Test1/'];
mkdir(figDir)
testEpochs = {'sleep'};
Nanim = numel(animIDs);

minOfData = 26; % cutoff first and last 2 minutes of data
eegTraces = cell(Nanim,1);
velTraces = cell(Nanim,1);
datIdx = cell(Nanim,1);

for k=1:Nanim
    tetinfo = load(sprintf('%s%stetinfo.mat',dataDirs{k},animIDs{k}));
    rntask  = load(sprintf('%s%srntask.mat',dataDirs{k},animIDs{k}));
    tetinfo = tetinfo.tetinfo;
    rntask  = rntask.rntask;
    tetIdx = getCellStructIdx(tetinfo,'descrip','riptet','area','CA1');
    taskIdx = getCellStructIdx(rntask,'epoch','sleep');
    days = taskIdx(:,1);
    epochs = taskIdx(:,2);
    tmp1 = arrayfun(@(x,y) any(days==x) & any(epochs==y),tetIdx(:,1),tetIdx(:,2));
    tetIdx = tetIdx(tmp1,:);
    
    eegTraces{k} = getAnimEEG(dataDirs{k},animIDs{k},tetIdx,'truncateDat',28);
    datIdx{k} = [repmat(k,size(tetIdx,1),1) tetIdx];
    velTraces{k} = getAnimVel(dataDirs{k},animIDs{k},tetIdx,'truncateDat',28);
end


eegDat = cell2mat(eegTraces);
velDat = cell2mat(velTraces);
eegFs = 1500; % TODO: fix hardcoding
velFs = 25; % TODO: fix hardcoding
dataIdx = cell2mat(datIdx);

immThresh = linspace(0,2,10);
Nthresh = numel(immThresh);
specFigures = [];
freqRanges = {[0.5 20],[40 80]};

% output variables
Nepisodes = zeros(Nthresh,3); %column for each behavioral state: sleeping,resting,moving
specVar = cell(numel(freqRanges),1);
for l=1:numel(freqRanges)
    specVar{l} = zeros(Nthresh,3);
end
for k=1:Nthresh
    iTh = immThresh(k);
    mTh = immThresh(k); % later make this vary as well
    [Nepisodes(k,:),sV,specFigures] = getStateSpectralVariance(eegDat,eegFs,velDat,velFs,'freqRanges',freqRanges,'immThresh',iTh,'movThresh',mTh,'specFigures',specFigures,'figRows',Nthresh,'currFigRow',k);
    for l=1:numel(freqRanges)
        specVar{l}(k,:)=sV(l,:);
    end
end
save([figDir 'ImmobileThresh-SpectralVariance-Test1.mat'],'Nepisodes','specVar','immThresh','freqRanges','animIDs','projDir','dataDirs','figDir','epochs','dataIdx')

% save specFigures
fname = @(x) sprintf('%sStateSpectra-VariedImmThresh-freq%g',figDir,x);
for k=1:numel(specFigures)
    figure(specFigures(k))
    suptitle(sprintf('Mean LFP spectra: %g - %g Hz',freqRanges{k}))
    saveas(specFigures(k),fname(k),'svg')
    varianceFig = figure();
    setFigureProperties(varianceFig,'Position',[1 1 1440 400])
    subplot(1,2,1)
    plot(immThresh,specVar{k})
    legend('sleeping','resting','moving')
    ylabel('Mean Spectral Variance')
    xlabel('Immobility Threshold (cm/s)')
    title({'Spectral variance',sprintf('%g-%g Hz',freqRanges{k})})
    subplot(1,2,2)
    plot(immThresh,Nepisodes)
    legend('sleeping','resting','moving')
    ylabel('Total number of episodes')
    xlabel('Immobility Threshold (cm/s)') 
    title('Number of episodes per state')
    suptitle('Varying Velocity Threshold for Immobility')
    saveas(varianceFig,sprintf('%sVariancePlot-VariedImmThresh-freqRange%g',figDir,k),'svg')
    close(varianceFig)
    close(specFigures(k))
end


