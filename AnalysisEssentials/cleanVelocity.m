function [outVel,outTime] = cleanVelocity(varargin)
    % [outVel,outTime] = cleanVelocity(velocity,time,NAME,VALUE) with take a velocity
    % vector and interpolate of nan and inf values, low-pass filter the
    % velocity and then resample to a consistent sampling rate using linear
    % interpolation. time should be in seconds and velocity in cm/s or px/s. 
    % Name-Value pairs:
    %   - fs: initial sampling rate in Hz (default 30Hz)
    %   - lpf: low-pass filter, structure with fields b and a containing filter
    %          parameters (default 5Hz 5th-order butterworth low-pass filter)
    %   - newFs: sampling rate to resample to in Hz (default 25 Hz)

    vel = varargin{1};
    T = varargin{2};
    fs = 30;
    [b,a] = butter(5,5/(fs/2),'low');
    lpf.b = b;
    lpf.a = a;
    clear a b
    newFs = 25;
    args = varargin(3:end);
    if ~isempty(args)
        assignVars(args{:});
    end

    outTime = T(1):1/newFs:T(end);
    infIdx = ~isfinite(vel);
    outVel = vel;
    outVel(infIdx) = interp1(T(~infIdx),vel(~infIdx),T(infIdx));
    outVel = filtfilt(lpf.b,lpf.a,outVel);
    [T,idx] = unique(T);
    outVel = interp1(T,outVel(idx),outTime,'linear','extrap');






