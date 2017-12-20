function out = checkDiodeNumbers(expDir)
a = subdir([expDir '*.videoPositionTracking']);
out = {};
for k=1:numel(a)
    b = readTrodesExtractedDataFile(a(k).name);
    out = [out;{a(k).folder(end-8:end) 1+any(b.fields(4).data)}];
end
