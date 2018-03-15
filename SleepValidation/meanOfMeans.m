function [MoM,EoM] = meanOfMeans(dat,err)
    MoM = mean(dat);
    EoM = sqrt(sum(err.^2))/numel(dat);
