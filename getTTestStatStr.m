function outStr = getTTestStatStr(alpha,P,CI,Stats,groupNames)

    if P<alpha && all(CI>0)
        compSign = '>';
        hypothesis = sprintf('Null Rejected: %s is significantly different from %s',groupNames{1},groupNames{2});
    elseif P<alpha && all(CI<0)
        compSign = '<';
        hypothesis = sprintf('Null Rejected: %s is significantly different from %s',groupNames{1},groupNames{2});
    else
        compSign = '=';
        hypothesis = sprintf('Null Accepted: %s is not significantly different from %s',groupNames{1},groupNames{2});
    end

    outStr = sprintf('%s\n%s    %s  %s\n@ the %0.4g%% confidence level\nTwo-Tailed t-Test\np-value : %0.03g\n%0.4g%% CI : [%0.4g %0.4g]\nt-stat : %0.4g\ndf : %0.4g\nsd : %0.4g',...
                    hypothesis,groupNames{1},compSign,groupNames{2},1-alpha,P,1-alpha,CI(1),CI(2),Stats.tstat,Stats.df,Stats.sd);
    outStr = strjust(char(strsplit(outStr,newline)),'center');

