function peaks = findPeaks(scan, options)
%FINDPEAKS  Detect diffraction peaks using MATLAB's findpeaks.
%   peaks = xrdc.peaks.findPeaks(scan)
%   peaks = xrdc.peaks.findPeaks(scan, Name, Value, ...)
%
%   Default analysis path. Wraps MATLAB's `findpeaks` with XRD-tuned
%   defaults — prominence-based detection is far more robust than the
%   slope-threshold detectors in the Delphi original and exposes FWHM
%   directly. For bit-for-bit reproduction of old Delphi analyses, use
%   `xrdc.peaks.findPeaksLegacy` instead.
%
%   Input
%     scan  : scan struct from any xrdc.io reader
%             (requires .twoTheta, .counts).
%
%   Name/Value options
%     'MinProminence'  (default: automatic, log-domain)
%         Required prominence above neighbouring troughs, in counts.
%         When omitted, an automatic scale-invariant criterion is used:
%         a peak must rise >= 0.3 decades above its neighbouring troughs
%         on the log10(counts) curve AND clear a Poisson significance test
%         (>= 5 sigma) against the local (~1 degree) median background. The
%         significance test is computed in photon units (intensity / count
%         quantum), so it is correct for both raw photon counts and Rigaku
%         counts-per-second exports, where one photon is ~5-15 cps and the
%         baseline is mostly exact zeros — single-photon cps blips are
%         rejected as ~1 sigma rather than mistaken for peaks. One default
%         covers a 1e6-count substrate peak and a 1e3-count film peak in
%         the same scan. Pass a numeric value (counts) for the classic
%         fixed linear threshold.
%     'MinHeight'      (default: -Inf)
%         Absolute count threshold (roughly equivalent to SpinEditMinPeak
%         in the Delphi UI).
%     'MinSeparation'  (default: 0.05°)
%         Minimum 2θ separation between accepted peaks. Peaks closer than
%         this are merged (the taller one wins). Matches FloatEdit13 in
%         the Delphi UI.
%     'MinWidth'       (default: 0° — no lower bound)
%         Minimum FWHM in 2θ to accept a candidate.
%     'MaxWidth'       (default: Inf)
%         Upper bound on FWHM — useful to reject broad background humps.
%     'TwoThetaRange'  (default: full scan)
%         [min, max] in degrees. Peaks outside are ignored.
%     'WidthReference' (default: "halfheight")
%         Passed through to `findpeaks`; either "halfheight" (FWHM) or
%         "halfprom" (prominence half-width).
%
%   Returns (struct array, one entry per detected peak)
%     .twoTheta        Peak 2θ (degrees)
%     .counts          Intensity at peak (same units as scan.counts)
%     .prominence      findpeaks prominence
%     .fwhm            Full width at half maximum (degrees)
%     .leftHalf        Left 2θ where profile crosses half-max
%     .rightHalf       Right 2θ where profile crosses half-max
%     .index           Index into scan.twoTheta (nearest; sub-sample
%                      positions are recorded in .twoTheta)
%
%   Requires the Signal Processing Toolbox (for `findpeaks`).
%
%   See also: xrdc.peaks.findPeaksLegacy, xrdc.peaks.adjustPeaks,
%             xrdc.peaks.fitPeak.

    arguments
        scan                    (1,1) struct
        options.MinProminence   (1,1) double = NaN   % NaN → auto
        options.MinHeight       (1,1) double = -Inf
        options.MinSeparation   (1,1) double {mustBeNonnegative} = 0.05
        options.MinWidth        (1,1) double {mustBeNonnegative} = 0
        options.MaxWidth        (1,1) double {mustBePositive}    = Inf
        options.TwoThetaRange   (1,2) double = [-Inf, Inf]
        options.WidthReference  (1,1) string {mustBeMember( ...
            options.WidthReference, ["halfheight","halfprom"])} = "halfheight"
    end
    if ~isnan(options.MinProminence) && options.MinProminence < 0
        error('xrdc:peaks:badProminence', ...
            'MinProminence must be non-negative.');
    end

    if ~isfield(scan, 'twoTheta') || ~isfield(scan, 'counts')
        error('xrdc:peaks:badScan', ...
            'scan must have .twoTheta and .counts fields.');
    end

    x = double(scan.twoTheta(:));
    y = double(scan.counts(:));
    if numel(x) ~= numel(y)
        error('xrdc:peaks:sizeMismatch', ...
            'twoTheta and counts must have the same length.');
    end
    if numel(x) < 3
        peaks = emptyPeakArray();
        return
    end

    % Restrict to requested 2θ window before calling findpeaks so the
    % MinPeakDistance constraint is interpreted in-window.
    inWin = x >= options.TwoThetaRange(1) & x <= options.TwoThetaRange(2);
    if nnz(inWin) < 3
        % findpeaks needs >= 3 samples; a too-narrow window means no peaks.
        peaks = emptyPeakArray();
        return
    end
    xw = x(inWin);
    yw = y(inWin);
    idxInScan = find(inWin);

    % MinPeakDistance is in samples (findpeaks semantics); convert from
    % 2θ using the median step. For irregular spacing this is still a
    % reasonable default.
    step = median(diff(xw));
    if step <= 0 || ~isfinite(step)
        minDistSamples = 1;
    else
        minDistSamples = max(1, round(options.MinSeparation / step));
    end

    if isnan(options.MinProminence)
        % Automatic log-domain criterion (see help text above).
        [pks, locs, widths, proms] = autoDetect(yw, step, minDistSamples, options);
    else
        fpArgs = { ...
            'MinPeakProminence', options.MinProminence, ...
            'MinPeakHeight',     options.MinHeight, ...
            'MinPeakDistance',   minDistSamples, ...
            'WidthReference',    char(options.WidthReference)};
        if options.MinWidth > 0
            fpArgs = [fpArgs, {'MinPeakWidth', options.MinWidth / max(step,eps)}];
        end
        if isfinite(options.MaxWidth)
            fpArgs = [fpArgs, {'MaxPeakWidth', options.MaxWidth / max(step,eps)}];
        end
        [pks, locs, widths, proms] = callFindpeaks(yw, fpArgs{:});
    end

    % findpeaks returns `locs` as indices into yw (the cropped vector).
    if isempty(pks)
        peaks = emptyPeakArray();
        return
    end

    % Map back to 2θ; findpeaks returns integer indices so we can resolve
    % the sub-sample peak position later with adjustPeaks or fitPeak.
    twoThetaPeak = xw(locs);
    fwhmDeg      = widths * step;
    idxGlobal    = idxInScan(locs);

    peaks = repmat(blankPeak(), numel(pks), 1);
    for i = 1:numel(pks)
        peaks(i).twoTheta   = twoThetaPeak(i);
        peaks(i).counts     = pks(i);
        peaks(i).prominence = proms(i);
        peaks(i).fwhm       = fwhmDeg(i);
        peaks(i).leftHalf   = twoThetaPeak(i) - fwhmDeg(i) / 2;
        peaks(i).rightHalf  = twoThetaPeak(i) + fwhmDeg(i) / 2;
        peaks(i).index      = idxGlobal(i);
    end

    % Sort ascending in 2θ (findpeaks returns sorted by sample index which
    % for ascending x is already sorted; guard against descending scans).
    [~, order] = sort([peaks.twoTheta]);
    peaks = peaks(order);
end

% -------------------------------------------------------------------------

function [pks, locs, widths, proms] = autoDetect(yw, step, minDistSamples, options)
%AUTODETECT  Log-domain auto prominence (MinProminence omitted).
%   A peak must (a) rise PROM_DECADES above its neighbouring troughs on
%   the log10 intensity curve — scale-invariant across the ~6-decade
%   dynamic range of a θ-2θ scan — and (b) clear a Poisson significance
%   test against the local median background. The significance test is
%   computed in *photon* units, not raw intensity units, so it is correct
%   regardless of whether the data are raw counts or counts-per-second
%   (cps). See estimateQuantum below. Reported metrics stay linear-domain:
%   log10 (clamped at 1) is monotone non-decreasing, so every accepted
%   log-domain maximum is a linear local maximum at the same index.
    PROM_DECADES  = 0.3;   % ≈2x above the surrounding troughs (log shape test)
    NOISE_SIGMAS  = 5;     % min Poisson significance (z) of the photon excess
    BG_WINDOW_DEG = 1.0;   % local-background median window, in 2θ

    pks = []; locs = []; widths = []; proms = [];

    yLog = log10(max(yw, 1));
    [~, locsLog] = callFindpeaks(yLog, ...
        'MinPeakProminence', PROM_DECADES, ...
        'MinPeakDistance',   minDistSamples);
    if isempty(locsLog), return, end

    if step > 0 && isfinite(step)
        bgWin = max(3, round(BG_WINDOW_DEG / step));
    else
        bgWin = min(numel(yw), 51);
    end
    bg  = movmedian(yw, bgWin);

    % Unit-aware Poisson guard. Estimate the count quantum q (intensity per
    % detected photon: 1 for raw counts, ~5–15 for cps), convert peak and
    % background to photon counts N = y/q, Nbg = bg/q, and require the
    % excess to be significant as a difference of two Poisson variables:
    %   z = (N - Nbg) / sqrt(N + Nbg + 1)  >=  NOISE_SIGMAS .
    % Using N (the peak's own counts) in the variance — rather than only the
    % background — is what rejects single-photon cps blips on a near-zero
    % baseline: one photon is z ≈ 1, far below threshold, whereas the old
    % sqrt(max(bg,1)) guard floored sigma at 1 *count* and let a 5–15 cps
    % blip clear a flat 5σ line. The +1 keeps z finite at zero background.
    q   = estimateQuantum(yw);
    N   = yw(locsLog) / q;
    Nbg = bg(locsLog) / q;
    z   = (N - Nbg) ./ sqrt(N + Nbg + 1);
    sig = z >= NOISE_SIGMAS;
    locsKeep = locsLog(sig);
    if isempty(locsKeep), return, end

    % Linear metrics: enumerate all linear local maxima with widths and
    % prominences, then keep the accepted subset by index.
    [pksAll, locsAll, widthsAll, promsAll] = callFindpeaks(yw, ...
        'MinPeakHeight',  options.MinHeight, ...
        'WidthReference', char(options.WidthReference));
    [locs, ia] = intersect(locsAll, locsKeep);
    pks = pksAll(ia); widths = widthsAll(ia); proms = promsAll(ia);

    keep = true(size(locs));
    if options.MinWidth > 0
        keep = keep & widths >= options.MinWidth / max(step, eps);
    end
    if isfinite(options.MaxWidth)
        keep = keep & widths <= options.MaxWidth / max(step, eps);
    end
    pks = pks(keep); locs = locs(keep); widths = widths(keep); proms = proms(keep);
end

function q = estimateQuantum(y)
%ESTIMATEQUANTUM  Intensity per detected photon (the count quantum).
%   Raw photon-count data has quantum 1. Rigaku cps exports divide raw
%   counts by the per-point integration time, so one photon shows up as a
%   fixed small intensity (~4.5–15 cps here) and the baseline is dominated
%   by exact zeros. The Poisson significance guard needs photon counts
%   N = y/q to be unit-correct, so estimate q here.
%
%   Estimator: when a substantial fraction of samples are *exactly* zero
%   (the signature of quantised cps data — continuous/noisy data almost
%   never lands on exact 0), the single-photon level is the small recurring
%   floor among the positive values; the 1st percentile of positive values
%   is a robust, outlier-resistant estimate of it. Otherwise (continuous or
%   raw-count data) fall back to q = 1, which reduces the guard to the
%   ordinary raw-count Poisson test. q is floored at 1 so the guard can
%   never become *more* permissive than raw-count Poisson.
%
%   The percentile is computed by sorting (base MATLAB) rather than
%   prctile, so this never pulls in the Statistics Toolbox.
    y   = y(:);
    pos = sort(y(y > 0));
    if numel(pos) < 20
        q = 1;
        return
    end
    if mean(y == 0) >= 0.05
        % 1st percentile by nearest-rank on the sorted positives.
        qc = pos(max(1, ceil(0.01 * numel(pos))));
        if qc > 1
            q = qc;
            return
        end
    end
    q = 1;
end

function [pks, locs, widths, proms] = callFindpeaks(y, varargin)
%CALLFINDPEAKS  Signal Processing Toolbox findpeaks, or the fallback.
    if isempty(which('findpeaks'))
        [pks, locs, widths, proms] = xrdc.peaks.findpeaks_fallback(y, varargin{:});
    else
        [pks, locs, widths, proms] = findpeaks(y, varargin{:});
    end
end

% -------------------------------------------------------------------------

function s = blankPeak()
    % Scalar peak struct with the canonical fields, all NaN.
    s = struct( ...
        'twoTheta',   NaN, ...
        'counts',     NaN, ...
        'prominence', NaN, ...
        'fwhm',       NaN, ...
        'leftHalf',   NaN, ...
        'rightHalf',  NaN, ...
        'index',      NaN);
end

function s = emptyPeakArray()
    % 0x1 struct array with the canonical fields.
    s = repmat(blankPeak(), 0, 1);
end
