function tetIdx = getTetIdx(dataDir,animID,tetReq,taskReq)
    % tetIdx = getTetIdx(dataDir,animID,tetReq,taskReq) returns the [day epoch tet] coordinates that match the requirements
    % tetReq is a cell array structured {fieldname,requirement,fieldname2,req2,...} where the fieldnames are fieldnames in tetinfo
    % taskReq is formatted the same except the structure searched is rntask

    tetinfo = load(sprintf('%s%stetinfo.mat',dataDir,animID));
    rntask = load(sprintf('%s%srntask.mat',dataDir,animID));
    tetinfo = tetinfo.tetinfo;
    rntask = rntask.rntask;

    tetIdx = getCellStructIdx(tetinfo,tetReq{:});
    taskIdx = getCellStructIdx(rntask,taskReq{:});
    days = taskIdx(:,1);
    epochs = taskIdx(:,2);

    idx = arrayfun(@(x,y) any(x==days) & any(y==epochs),tetIdx(:,1),tetIdx(:,2));
    tetIdx = tetIdx(idx,:);

