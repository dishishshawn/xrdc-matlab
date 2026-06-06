function info = dominantPeriod(twoTheta, yDetrended, options)
%DOMINANTPERIOD  FFT-based dominant fringe period of a detrended XRR signal.
%   info = xrdc.xrr.dominantPeriod(twoTheta, yDetrended)
%   info = xrdc.xrr.dominantPeriod(twoTheta, yDetrended, Name=Value)
%
%   Picks the strongest peak in |FFT(yDetrended)| inside a physically
%   reasonable band. Used to (a) report a robust thickness estimate that
%   doesn't depend on individual-peak detection, and (b) set an informed
%   MinSeparation for findPeaks so the peak detector doesn't pick up
%   spurious high-frequency residuals from imperfect detrending.
%
%   Inputs
%     twoTheta    : 2θ axis (degrees), uniformly spaced
%     yDetrended  : background-subtracted log-counts (or any detrended
%                   trace) on the same grid
%
%   Name/Value options
%     'MinFreq' (default  0.5 cycles/°)
%         Excludes DC and the lowest band where polynomial-detrend
%         residuals tend to leak. 0.5/° corresponds to a fringe period
%         of 2° → film thicknesses well below ~5 nm map here.
%     'MaxFreq' (default 30 cycles/°)
%         Upper bound — past this band, Rigaku step-noise dominates.
%         30/° ≈ films above ~300 nm; rarely relevant for thin-film XRR.
%
%   Output (struct)
%     .freqPerDeg     dominant frequency (cycles per degree of 2θ)
%     .periodDeg      1 / freqPerDeg (degrees per fringe)
%     .amplitude      |FFT| at the peak
%     .signalToNoise  peak amplitude / median(|FFT|) in band

    arguments
        twoTheta            (:,1) double
        yDetrended          (:,1) double
        options.MinFreq     (1,1) double {mustBePositive} = 0.5
        options.MaxFreq     (1,1) double {mustBePositive} = 30
    end

    if numel(twoTheta) ~= numel(yDetrended)
        error('xrdc:xrr:sizeMismatch', ...
            'twoTheta and yDetrended must have the same length.');
    end
    if numel(twoTheta) < 16
        error('xrdc:xrr:tooFewPoints', ...
            'Need at least 16 points for an FFT-based period estimate.');
    end

    step = median(diff(twoTheta));
    if step <= 0 || ~isfinite(step)
        error('xrdc:xrr:badStep', 'twoTheta must be strictly increasing.');
    end

    n   = numel(yDetrended);
    Y   = abs(fft(yDetrended - mean(yDetrended)));
    freq = (0:n-1).' / (step * n);

    half = floor(n / 2);
    Y    = Y(2:half);          % drop DC
    freq = freq(2:half);

    band = freq >= options.MinFreq & freq <= options.MaxFreq;
    if ~any(band)
        error('xrdc:xrr:noBand', ...
            'No FFT bins fall in [%.2f, %.2f] cycles/° — adjust MinFreq/MaxFreq.', ...
            options.MinFreq, options.MaxFreq);
    end
    YB = Y(band);
    fB = freq(band);

    [pkAmp, iPk] = max(YB);
    fPeak = fB(iPk);
    snr   = pkAmp / max(median(YB), eps);

    info = struct( ...
        'freqPerDeg',    fPeak, ...
        'periodDeg',     1 / fPeak, ...
        'amplitude',     pkAmp, ...
        'signalToNoise', snr);
end
