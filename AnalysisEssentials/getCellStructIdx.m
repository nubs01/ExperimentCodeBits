function structIdx = getCellStructIdx(cellStruct,varargin)
% structIdx = getCellStructIdx(cellStruct,NAME,VALUE) 
% eg. structIdx = getCellStructIdx(cellStruct,'descrip','riptet','area','CA1')

if mod(numel(varargin),2) ~=0
    error('all search tags must be paired with a value')
end
tags = varargin(1:2:numel(varargin));
searchVals = varargin(2:2:numel(varargin));
structIdx = [];
for k=1:numel(tags)
    A = cellfetch(cellStruct,tags{k});
    sv = searchVals{k};
    if all(cellfun(@isempty,A.values))
        disp([tags{k} ' is not a valid field of the structure. Skipping...'])
        continue;
    end
    if ischar(sv)
        idx  = find(strcmp(A.values,sv));
    else
        idx = find(A.values==sv);
    end

    if isempty(structIdx)
        structIdx = A.index(idx,:);
    else
        structIdx = intersect(structIdx,A.index(idx,:),'rows');
    end
end
