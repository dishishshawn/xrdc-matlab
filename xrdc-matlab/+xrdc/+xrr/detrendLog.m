function yDet = detrendLog(twoTheta, yLog, options)
%DETRENDLOG  Polynomial detrend of log10(counts) for XRR fringe analysis.
%   yDet = xrdc.xrr.detrendLog(twoTheta, yLog)
%   yDet = xrdc.xrr.detrendLog(twoTheta, yLog, Order=4)
%
%   The slow XRR decay below the critical edge follows roughly
%   I ∝ (θ_c/θ)^4 with an additional absorption falloff; in log space
%   that's smooth and well captured by a polynomial of moderate order.
%   A windowed sgolayfilt envelope (the previous approach) absorbed
%   fringes for thin films whenever the envelope span approached the
%   fringe period, causing ~3× thickness errors. A polynomial subtraction
%   sidesteps the period-vs-window trap.
%
%   Inputs
%     twoTheta : 2θ axis (degrees)
%     yLog     : log10(counts) on the same grid (already post-edge)
%
%   Name/Value
%     'Order'  (default 4) — polynomial order
%
%   Output
%     yDet     : yLog minus polyval(p, twoTheta), same size as inputs

    arguments
        twoTheta        (:,1) double
        yLog            (:,1) double
        options.Order   (1,1) double {mustBeInteger, mustBePositive} = 4
    end

    if numel(twoTheta) ~= numel(yLog)
        error('xrdc:xrr:sizeMismatch', ...
            'twoTheta and yLog must have the same length.');
    end
    if numel(twoTheta) <= options.Order + 1
        error('xrdc:xrr:tooFewPoints', ...
            'Need more than Order+1 points (got %d, Order=%d).', ...
            numel(twoTheta), options.Order);
    end

    p    = polyfit(twoTheta, yLog, options.Order);
    yDet = yLog - polyval(p, twoTheta);
end
