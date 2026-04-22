function peaks = findPeaksLegacy(scan, options)
%FINDPEAKSLEGACY  Slope/2nd-derivative peak detectors ported from Delphi.
%   peaks = xrdc.peaks.findPeaksLegacy(scan)
%   peaks = xrdc.peaks.findPeaksLegacy(scan, 'Algorithm','slope'|'slope2', ...)
%
%   Ports the Delphi ScanPeak1/ScanPeak2 detectors from xrdc1.pas:1482 and
%   1542. Preserved for bit-for-bit reproduction of old XRDC analyses —
%   prefer `xrdc.peaks.findPeaks` for new work (prominence-based detection
%   is more robust and returns FWHM directly).
%
%   Input
%     scan : scan struct. Requires .twoTheta, .counts, .slope, .slope2.
%            If .slope / .slope2 are missing, they are computed on the
%            fly via `xrdc.signal.derivatives` with defaults.
%
%   Name/Value options
%     'Algorithm'       ("slope" | "slope2"), default "slope".
%                       "slope"  = ScanPeak1 (xrdc1.pas:1482).
%                       "slope2" = ScanPeak2 (xrdc1.pas:1542).
%     'SlopeThreshold'  Default 50. Minimum positive slope magnitude to
%                       accept a peak candidate (ScanPeak1: maxSlope;
%                       ScanPeak2: |slope2|). Matches SpinEdit1/SpinEdit6.
%     'MinHeight'       Default 0. Absolute count threshold applied after
%                       the detector runs (Delphi SpinEditMinPeak).
%     'MinSeparation'   Default 0.05°. Peaks closer than this are merged
%                       (counts-weighted average) — matches FloatEdit13
%                       and the Delphi merge loop in `FindPeaks`.
%     'TwoThetaRange'   Default [-Inf, Inf].
%
%   Returns the same struct shape as xrdc.peaks.findPeaks, minus .fwhm
%   and the half-height positions (the slope detectors don't compute
%   them — run `xrdc.peaks.adjustPeaks` to refine position and FWHM).
%
%   See also: xrdc.peaks.findPeaks, xrdc.peaks.adjustPeaks.

    arguments
        scan                       (1,1) struct
        options.Algorithm          (1,1) string {mustBeMember( ...
            options.Algorithm, ["slope","slope2"])} = "slope"
        options.SlopeThreshold     (1,1) double {mustBePositive} = 50
        options.MinHeight          (1,1) double = 0
        options.MinSeparation      (1,1) double {mustBeNonnegative} = 0.05
        options.TwoThetaRange      (1,2) double = [-Inf, Inf]
    end

    x = double(scan.twoTheta(:));
    y = double(scan.counts(:));
    if numel(x) < 5
        peaks = emptyPeakArray();
        return
    end

    % Derivatives: use what's on the scan if present, otherwise compute
    % them with default smoothing parameters.
    if isfield(scan, 'slope') && isfield(scan, 'slope2') ...
            && numel(scan.slope) == numel(x)
        slope  = double(scan.slope(:));
        slope2 = double(scan.slope2(:));
    else
        [slope, slope2] = xrdc.signal.derivatives(x, y);
    end

    % Run the requested detector over the whole trace, recursively
    % restarting after each hit (this matches the Delphi loop in
    % FindPeaks at xrdc1.pas:1589-1605).
    n = numel(x);
    j = 1;
    hits = zeros(0, 1);
    while true
        switch options.Algorithm
            case "slope"
                i = scanPeak1(j, n, slope, options.SlopeThreshold);
            case "slope2"
                i = scanPeak2(j, n, slope, slope2, options.SlopeThreshold);
        end
        if i < 1, break; end
        hits(end+1, 1) = i;   %#ok<AGROW>
        j = i + 1;
        if j > n, break; end
    end

    if isempty(hits)
        peaks = emptyPeakArray();
        return
    end

    % Apply MinHeight and TwoThetaRange filters
    tt = x(hits);
    cc = y(hits);
    keep = cc >= options.MinHeight & ...
           tt >= options.TwoThetaRange(1) & ...
           tt <= options.TwoThetaRange(2);
    hits = hits(keep);
    tt   = tt(keep);
    cc   = cc(keep);

    if isempty(hits)
        peaks = emptyPeakArray();
        return
    end

    % Sort by 2θ (the detector already walks left-to-right, but guard
    % against descending scans just in case).
    [tt, order] = sort(tt);
    cc          = cc(order);
    hits        = hits(order);

    % Merge close pairs — matches the Delphi "Peaks zu dicht beisammen"
    % branch at xrdc1.pas:1656. The merged position is the simple mean
    % (not counts-weighted) to exactly reproduce Delphi behaviour; the
    % intensity is the average of the pair.
    merged = false(size(tt));
    i = 1;
    while i < numel(tt)
        if ~merged(i) && abs(tt(i+1) - tt(i)) < options.MinSeparation
            tt(i)   = (tt(i) + tt(i+1)) / 2;
            cc(i)   = (cc(i) + cc(i+1)) / 2;
            hits(i) = round((hits(i) + hits(i+1)) / 2);
            merged(i+1) = true;
        end
        i = i + 1;
    end
    tt   = tt(~merged);
    cc   = cc(~merged);
    hits = hits(~merged);

    peaks = repmat(blankPeak(), numel(tt), 1);
    for k = 1:numel(tt)
        peaks(k).twoTheta = tt(k);
        peaks(k).counts   = cc(k);
        peaks(k).index    = hits(k);
        % FWHM / half-positions are populated by adjustPeaks or fitPeak.
    end
end

% -------------------------------------------------------------------------

function idx = scanPeak1(startIdx, endIdx, slope, threshold)
% Direct port of TForm1.ScanPeak1 (xrdc1.pas:1482).
% 1-based index convention (Delphi uses 0-based, but the algorithm is
% the same).
    idx = -1;
    if startIdx >= endIdx, return; end

    % 1. find the slope maximum (must exceed +threshold)
    maxSlope = -Inf;
    maxPos   = -1;
    zeroCross = -1;
    for i = startIdx:endIdx
        if slope(i) > threshold
            if slope(i) > maxSlope
                maxSlope = slope(i);
                maxPos   = i;
            end
        elseif maxSlope > 0 && slope(i) < 0
            zeroCross = i;
            break;
        end
    end
    if maxPos < 1 || zeroCross < 1, return; end

    % 2. slope minimum beyond zero-crossing; min must be < -maxSlope/2
    se = zeroCross + 2 * (zeroCross - maxPos);
    if se > endIdx, se = endIdx; end
    minSlope = +Inf;
    minPos   = -1;
    for i = zeroCross:se
        if slope(i) < -maxSlope / 2
            if slope(i) < minSlope
                minSlope = slope(i);
                minPos   = i;
            end
        end
    end
    if minPos < 1
        % nothing found — recurse from zeroCross
        if zeroCross < endIdx
            idx = scanPeak1(zeroCross, endIdx, slope, threshold);
        end
        return
    end

    idx = floor((maxPos + minPos) / 2);
    if idx < startIdx, idx = startIdx; end
    if idx > endIdx,   idx = endIdx;   end
end

function idx = scanPeak2(startIdx, endIdx, slope, slope2, threshold)
% Direct port of TForm1.ScanPeak2 (xrdc1.pas:1542).
    idx = -1;
    if startIdx >= endIdx, return; end
    last = 0;
    for i = startIdx:endIdx
        signChange = (last < 0 && slope(i) > 0) || ...
                     (last > 0 && slope(i) < 0) || ...
                     (slope(i) == 0);
        if signChange && slope2(i) < 0 && -slope2(i) >= threshold
            idx = i;
            return
        end
        last = slope(i);
    end
end

% -------------------------------------------------------------------------

function s = blankPeak()
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
    s = repmat(blankPeak(), 0, 1);
end
