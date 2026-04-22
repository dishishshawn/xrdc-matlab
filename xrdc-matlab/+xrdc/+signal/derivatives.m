function [slope, slope2] = derivatives(twoTheta, counts, frameSize, polyOrder)
%DERIVATIVES  First and second derivatives of a counts trace.
%   [slope, slope2] = xrdc.signal.derivatives(twoTheta, counts)
%   [slope, slope2] = xrdc.signal.derivatives(twoTheta, counts, frameSize)
%   [slope, slope2] = xrdc.signal.derivatives(twoTheta, counts, frameSize, polyOrder)
%
%   Modern replacement for the sliding-window OLS fit used in the Delphi
%   XRDC (CalcSlopes / Geradenanpassung in xrdc1.pas:1353+). Uses
%   Savitzky-Golay, which produces the same result as a sliding OLS fit
%   of the given polynomial order but is better-conditioned.
%
%   Defaults: frameSize = 11 points, polyOrder = 3 — appropriate for
%   typical XRD step sizes (0.01°–0.02°).
%
%   Inputs
%     twoTheta   : 2θ vector (degrees); assumed uniformly spaced.
%     counts     : counts vector, same length as twoTheta
%     frameSize  : odd integer > polyOrder (window width)
%     polyOrder  : polynomial order (typically 2 or 3)
%
%   Outputs
%     slope   : dCounts/d(2θ)   (counts per degree)
%     slope2  : d²Counts/d(2θ)² (counts per degree²)
%
%   Note on edges: the first and last (frameSize-1)/2 points of slope and
%   slope2 use MATLAB's default sgolayfilt end behaviour, which is a
%   lower-order polynomial fit to the available points. For peak
%   detection you probably want to ignore derivative values near the
%   very edges of the scan anyway.

    arguments
        twoTheta   (:,1) double
        counts     (:,1) double
        frameSize  (1,1) double {mustBeInteger, mustBePositive} = 11
        polyOrder  (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    end

    if numel(twoTheta) ~= numel(counts)
        error('xrdc:signal:sizeMismatch', ...
            'twoTheta and counts must be the same length.');
    end

    if mod(frameSize, 2) == 0
        frameSize = frameSize + 1;
    end
    if frameSize <= polyOrder
        error('xrdc:signal:badFrame', ...
            'frameSize (%d) must exceed polyOrder (%d).', frameSize, polyOrder);
    end

    step = mean(diff(twoTheta));
    if any(abs(diff(twoTheta) - step) > 1e-6 * max(abs(step), 1))
        warning('xrdc:signal:nonUniform', ...
            'twoTheta is not uniformly spaced; derivatives will be approximate.');
    end

    % Use Signal Processing Toolbox sgolay if available; otherwise the
    % pure-MATLAB fallback. Same coefficient convention either way.
    if isempty(which('sgolay'))
        [~, g] = xrdc.signal.sgolay_fallback(polyOrder, frameSize);
    else
        [~, g] = sgolay(polyOrder, frameSize);
    end
    % g(:, k+1) gives Savitzky-Golay coefficients for the k-th derivative
    % scaled so that y^{(k)} ≈ factorial(k) * (g(:,k+1).' * window) / step^k.
    slope  = conv(counts, factorial(1) * flipud(g(:, 2)), 'same') / step;
    slope2 = conv(counts, factorial(2) * flipud(g(:, 3)), 'same') / step^2;
end
