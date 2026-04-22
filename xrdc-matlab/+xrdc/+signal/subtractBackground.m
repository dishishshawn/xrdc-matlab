function [y, baseline] = subtractBackground(counts, windowSize, method)
%SUBTRACTBACKGROUND  Subtract a moving-window baseline from counts.
%   [y, baseline] = xrdc.signal.subtractBackground(counts, windowSize)
%   [y, baseline] = xrdc.signal.subtractBackground(counts, windowSize, method)
%
%   Default method ('movmean') matches the legacy XRDC SubstractBackground
%   in xrdc1.pas:1457 — a large-window moving average that is subtracted
%   pointwise, with negative values clipped to zero.
%
%   Inputs
%     counts     : counts vector
%     windowSize : integer window size (points)
%     method     : 'movmean' (default) | 'movmin' | 'rollingPercentile'
%
%   Output
%     y        : counts with baseline removed, negative values clipped to 0
%     baseline : the baseline that was subtracted

    arguments
        counts     (:,1) double
        windowSize (1,1) double {mustBePositive, mustBeInteger}
        method     (1,:) char = 'movmean'
    end

    switch lower(method)
        case 'movmean'
            baseline = movmean(counts, windowSize, 'Endpoints', 'shrink');
        case 'movmin'
            baseline = movmin(counts, windowSize, 'Endpoints', 'shrink');
        case 'rollingpercentile'
            % 10th percentile in a rolling window — more robust to sharp peaks
            baseline = movmedian(counts, windowSize, 'Endpoints', 'shrink');
            % a simple rolling-min as a second-pass refinement:
            baseline = movmin(baseline, windowSize, 'Endpoints', 'shrink');
        otherwise
            error('xrdc:signal:unknownMethod', ...
                'Unknown background method: %s', method);
    end

    y = counts - baseline;
    y(y < 0) = 0;
end
