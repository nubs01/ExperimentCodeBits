function outStr = getTTestStatStr(alpha,P,CI,Stats,groupNames,N)

    if P<alpha && all(CI>0)
        compSign = '>';
        hypothesis = 'Null Rejected';
    elseif P<alpha && all(CI<0)
        compSign = '<';
        hypothesis = 'Null Rejected';
    else
        compSign = '=';
        hypothesis = 'Null Accepted';
    end

    outStr = sprintf('%s\n%s    %s  %s\n@ the %0.4g%% confidence level\nTwo-Tailed t-Test\np-value : %0.03g\n%0.4g%% CI : [%0.4g %0.4g]\nt-stat : %0.4g\ndf : %0.4g\nsd : %0.4g\nN: %i vs %i',...
                    hypothesis,groupNames{1},compSign,groupNames{2},(1-alpha)*100,P,(1-alpha)*100,CI(1),CI(2),Stats.tstat,Stats.df,Stats.sd,N(1),N(2));
    outStr = strjust(char(strsplit(outStr,newline)),'center');

