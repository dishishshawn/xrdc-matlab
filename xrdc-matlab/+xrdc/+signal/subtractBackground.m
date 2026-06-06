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
%                  | 'spline' — Curve Fitting Toolbox csaps smoothing
%                    spline fit to a moving-min envelope. Best for XRR
%                    backgrounds with curvature. Falls back to 'movmin'
%                    when the Curve Fitting Toolbox is missing.
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
        case 'spline'
            if isempty(which('csaps'))
                warning('xrdc:signal:noCurveFitting', ...
                    'csaps not available (Curve Fitting Toolbox missing). Falling back to movmin.');
                baseline = movmin(counts, windowSize, 'Endpoints', 'shrink');
            else
                % Fit a smoothing spline to the rolling-min envelope so
                % the resulting baseline is C² and free of step-like
                % artefacts at peak edges.
                env = movmin(counts, windowSize, 'Endpoints', 'shrink');
                idx = (1:numel(counts)).';
                h = double(windowSize);
                p = 1 / (1 + h^3 / 6);
                baseline = csaps(idx, env, p, idx);
                baseline = baseline(:);
            end
        otherwise
            error('xrdc:signal:unknownMethod', ...
                'Unknown background method: %s', method);
    end

    y = counts - baseline;
    y(y < 0) = 0;
end
