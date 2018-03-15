function permStats = permutationTest(X,Y,nboot,varargin)

    paired = 0;
    if isrow(X)
        X = X';
    end
    if isrow(Y)
        Y = Y';
    end
    bootfun = @(x,y) mean(x)-mean(y);
    assignVars(varargin)

    Nx = numel(X);
    Ny = numel(Y);
    if ~paired
        N = Nx+Ny;
        Z = [X;Y];
        permDiffs = zeros(nboot,1);
        testStat = bootfun(X,Y);
        for k=1:nboot
            idx = randperm(N);
            i1 = idx(1:Nx);
            i2 = idx(Nx+1:end);
            permDiffs(k) = bootfun(Z(i1),Z(i2));
        end
    else
        Z = [X Y];
        N = size(Z,1);
        permDiffs = zeros(nboot,1);
        testStat = bootfun(X,Y);
        for k=1:nboot
            idx = rand(N,1);
            Z(idx>.5,:) = fliplr(Z(idx>.5,:));
            permDiffs(k) = bootfun(Z(:,1),Z(:,2));
        end
    end
            

    if testStat > mean(permDiffs)
        Pone = sum(permDiffs>testStat)/nboot;
    else
        Pone = sum(permDiffs<testStat)/nboot;
    end
    Ptwo = sum(abs(permDiffs)>abs(testStat))/nboot;
    permStats = struct('P_one_tail',Pone,'P_two_tail',Ptwo,'permutation_stats',permDiffs);
