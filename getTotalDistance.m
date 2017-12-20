function out = getTotalDistance(x,y)
    % returns the total distance moved between sequential x-y points
    out = 0;
    for i=2:numel(x)
        out = out + sqrt((x(i)-x(i-1))^2+(y(i)-y(i-1))^2);
    end

