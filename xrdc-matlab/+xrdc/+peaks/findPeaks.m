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
%     'MinProminence'  (default: 1.5% of the scan's peak-to-trough range)
%         Required prominence above neighbouring troughs, in counts.
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

    % Default prominence = 1.5% of peak-to-trough range; gives reasonable
    % sensitivity without catching quantisation noise in low-count scans.
    minProm = options.MinProminence;
    if isnan(minProm)
        minProm = max(1, 0.015 * (max(y) - min(y)));
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

    fpArgs = { ...
        'MinPeakProminence', minProm, ...
        'MinPeakHeight',     options.MinHeight, ...
        'MinPeakDistance',   minDistSamples, ...
        'WidthReference',    char(options.WidthReference)};
    if options.MinWidth > 0
        fpArgs = [fpArgs, {'MinPeakWidth', options.MinWidth / max(step,eps)}];
    end
    if isfinite(options.MaxWidth)
        fpArgs = [fpArgs, {'MaxPeakWidth', options.MaxWidth / max(step,eps)}];
    end

    [pks, locs, widths, proms] = findpeaks(yw, fpArgs{:});

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
