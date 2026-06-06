function envelope = logEnvelope(yLog, spanPts, polyOrder)
%LOGENVELOPE  Smooth envelope of a log-counts trace for fringe detrending.
%   envelope = xrdc.signal.logEnvelope(yLog, spanPts)
%   envelope = xrdc.signal.logEnvelope(yLog, spanPts, polyOrder)
%
%   Removes the slow decay underneath log10(counts) in XRR / Kiessig data
%   so that each fringe contributes roughly the same ripple amplitude
%   regardless of absolute intensity. Used by demoXRR and the GUI's XRR
%   analysis path.
%
%   Uses sgolayfilt (Signal Processing Toolbox) when available — it
%   follows steep, curved decays much better than movmean and doesn't
%   distort the ripple. Falls back to a centred moving mean otherwise.
%
%   Inputs
%     yLog      : log10(counts) vector (or any 1-D signal to detrend)
%     spanPts   : window length in samples (forced odd internally)
%     polyOrder : Savitzky-Golay polynomial order (default 3)
%
%   Output
%     envelope  : the slow component, same size as yLog

    arguments
        yLog      (:,1) double
        spanPts   (1,1) double {mustBePositive, mustBeInteger}
        polyOrder (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    end

    n = numel(yLog);
    if n == 0
        envelope = yLog;
        return
    end

    win = spanPts;
    if mod(win, 2) == 0, win = win + 1; end
    if win > n
        win = max(3, n - 1 + mod(n, 2));   % largest odd ≤ n
    end

    if win <= polyOrder || isempty(which('sgolayfilt'))
        % Fallback: centred moving mean (Endpoints='shrink' keeps the
        % trace inside the data on both ends).
        envelope = movmean(yLog, spanPts, 'Endpoints', 'shrink');
        return
    end

    envelope = sgolayfilt(yLog, polyOrder, win);
end
