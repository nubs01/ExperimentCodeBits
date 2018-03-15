function out = interpInf(datVec)
    A = ~isfinite(datVec);
    T = 1:numel(datVec);
    out = datVec;
    if ~any(A)
        return;
    end
    out(A) = interp1(T(~A),datVec(~A),T(A));

