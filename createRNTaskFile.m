function createRNTaskFile(animID,expDir,dataDir,genotype,epochTypes,dayEpochs)
    rntask = cell(1,numel(dayEpochs));
    saveFile = [dataDir filesep animID 'rntask.mat'];
    if exist(saveFile,'file')
        return;
    end
    fprintf('Creating rntask file...\n')
    for k=1:numel(dayEpochs)
        rntask{k} = cell(1,numel(dayEpochs{k}));
        for l=1:numel(dayEpochs{k})
            rntask{k}{l} = struct('animal',animID,'genotype',genotype,'epoch',epochTypes{dayEpochs{k}(l)},'expDir',expDir);
        end
    end
    save(saveFile,'rntask')
