function setupExperimentDirectories(animID,descriptor)

if ~exist('animID','var')
    animID = inputdlg('Enter animal ID:','Setup Experiment Directories');
    animID = animID{1};
end

if ~exist('descriptor','var')
    descriptor = inputdlg('Enter Experiment Descriptor: (1-3 words, no spaces, upperCamelCase)','Experiment Descriptor')
    descriptor = descriptor{1};
end
topDir = sprintf('%s_%s_Experiment%s',animID,descriptor,filesep)
dirsToMake = {topDir,...
              sprintf('%s%s',topDir,animID),...
              sprintf('%s%s_direct',topDir,animID),...
              sprintf('%s%s_analysis',topDir,animID),...
              sprintf('%s%s_histology',topDir,animID),...
              sprintf('%s%s_preprocessCode',topDir,animID),...
              sprintf('%s%s_configs',topDir,animID)};

for i=1:numel(dirsToMake),
    if ~exist(dirsToMake{i},'dir')
        mkdir(dirsToMake{i});
    end
end
