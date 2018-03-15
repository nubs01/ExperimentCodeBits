function out = checkPosTrackingRange(expDir)
a = subdir([expDir filesep '*.videoPositionTracking']);
out = {};
for k=1:numel(a)
    b = readTrodesExtractedDataFile(a(k).name);
    out = [out;{a(k).folder(end-8:end) b.pixelscale}];
end
