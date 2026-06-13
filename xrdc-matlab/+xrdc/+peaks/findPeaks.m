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
%         on the log10(counts) curve AND stand >= 5 Poisson sigma above
%         the local (~1 degree) median background. One default covers a
%         1e6-count substrate peak and a 1e3-count film peak in the same
%         scan. Counts below 1 are clamped before the log, so sub-1-cps
%         baselines yield no peaks. Pass a numeric value (counts) for the
%         classic fixed linear threshold.
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
    if ~any(inWin)
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
%   test against the local median background, which rejects the large
%   fake log-prominence of quantisation jitter on near-zero-count
%   baselines. Reported metrics stay linear-domain: log10 (clamped at 1)
%   is monotone non-decreasing, so every accepted log-domain maximum is a
%   linear local maximum at the same index.
    PROM_DECADES  = 0.3;   % ≈2x above the surrounding troughs
    NOISE_SIGMAS  = 5;     % Poisson significance vs local background
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
    sig = (yw(locsLog) - bg(locsLog)) >= NOISE_SIGMAS * sqrt(max(bg(locsLog), 1));
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
