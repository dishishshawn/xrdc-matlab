function y = smoothCounts(counts, windowSize, method)
%SMOOTHCOUNTS  Smooth a counts trace.
%   y = xrdc.signal.smoothCounts(counts, windowSize)
%   y = xrdc.signal.smoothCounts(counts, windowSize, method)
%
%   Defaults to 'movmean' (matching the legacy XRDC NoiseSuppression).
%   method can be any option accepted by smoothdata():
%     'movmean' (default), 'movmedian', 'gaussian', 'lowess', 'loess',
%     'rlowess', 'rloess', 'sgolay'.
%
%   Inputs
%     counts     : counts vector
%     windowSize : integer window size (points)
%     method     : char/string, optional
%
%   Output
%     y          : smoothed counts, same size as input.

    arguments
        counts     (:,1) double
        windowSize (1,1) double {mustBePositive, mustBeInteger}
        method     (1,:) char = 'movmean'
    end

    y = smoothdata(counts, method, windowSize);
end
