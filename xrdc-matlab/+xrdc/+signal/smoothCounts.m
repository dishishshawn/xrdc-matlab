function y = smoothCounts(counts, windowSize, method)
%SMOOTHCOUNTS  Smooth a counts trace.
%   y = xrdc.signal.smoothCounts(counts, windowSize)
%   y = xrdc.signal.smoothCounts(counts, windowSize, method)
%
%   Defaults to 'movmean' (matching the legacy XRDC NoiseSuppression).
%   method can be any option accepted by smoothdata():
%     'movmean' (default), 'movmedian', 'gaussian', 'lowess', 'loess',
%     'rlowess', 'rloess', 'sgolay'
%   or the Curve Fitting Toolbox smoothing-spline path:
%     'csaps'  — Reinsch smoothing spline; windowSize sets the smoothing
%                length scale (in samples). Falls back to 'gaussian'
%                smoothing if the Curve Fitting Toolbox is unavailable.
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

    if strcmpi(method, 'csaps')
        if isempty(which('csaps'))
            warning('xrdc:signal:noCurveFitting', ...
                'csaps not available (Curve Fitting Toolbox missing). Falling back to gaussian smoothing.');
            y = smoothdata(counts, 'gaussian', windowSize);
            return
        end
        % Map windowSize → smoothing parameter p ∈ (0,1].
        % csaps' p is a balance between exact interpolation (p=1) and
        % a straight-line fit (p=0). A useful rule of thumb is
        %   p = 1 / (1 + h^3 / 6)
        % where h is the characteristic smoothing length in samples.
        % windowSize is taken as that length.
        h = double(windowSize);
        p = 1 / (1 + h^3 / 6);
        idx = (1:numel(counts)).';
        y = csaps(idx, counts, p, idx);
        y = y(:);
        return
    end

    y = smoothdata(counts, method, windowSize);
end
