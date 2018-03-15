function bootStruct = bootstrapDat(X,nboot,varargin)
    % bootStruct = bootstrapDat(X,nboot) will bootstrap each column in X and
    % return a struct array with the mean, SEM and bootstats for each column
    % can also pass NAME-VALUE pairs to override:
    %   - bootfun  : function to apply to data (default = @mean)

    bootfun = @mean;
    assignVars(varargin)

    Ncol = size(X,2);
    meanDat = cell(1,Ncol);
    semDat = cell(1,Ncol);
    bootstats = cell(1,Ncol);
    for k=1:Ncol
        y = X(:,k);
        meanDat{k} = mean(y);
        bootstats{k} = bootstrp(nboot,bootfun,y);
        semDat{k} = std(bootstats{k});
    end
    bootStruct = struct('mean',meanDat,'SEM',semDat,'bootstats',bootstats);
