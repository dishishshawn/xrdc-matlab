function peaksOut = adjustPeaks(scan, peaks, options)
%ADJUSTPEAKS  Refine peak positions and FWHMs by interpolated area-bisection.
%   refined = xrdc.peaks.adjustPeaks(scan, peaks)
%   refined = xrdc.peaks.adjustPeaks(scan, peaks, Name, Value, ...)
%
%   Ports the FWHM-based refinement from xrdc1.pas:1715 (`AdjustPeaks`).
%   For each peak the trace is interpolated to a 100× finer grid, the
%   half-height window is located, and the peak position is replaced by
%   the **area bisector** of that window (not the centroid, not the
%   parabola vertex — this is the Delphi convention, see
%   ALGORITHM_SPEC §4.4).
%
%   Input
%     scan       : scan struct with .twoTheta, .counts.
%     peaks      : struct array from findPeaks / findPeaksLegacy.
%
%   Name/Value options
%     'HeightFraction' (default 0.5)
%         Fraction of peak height used to define the refinement window.
%         0.5 → FWHM window. Matches SpinEdit3 / 100 in the Delphi UI.
%     'Oversample'     (default 100)
%         Interpolation factor. 100 matches the Delphi default.
%     'MergeTolerance' (default 0.03°)
%         After refinement, peaks closer than this in 2θ are collapsed to
%         a single peak (the higher-intensity one wins). Matches the
%         0.03° de-dup step at xrdc1.pas:1774.
%
%   Output
%     Same struct shape as the input; .twoTheta, .counts, .fwhm,
%     .leftHalf, .rightHalf are updated. .index is re-mapped to the
%     nearest original-grid sample.
%
%   See also: xrdc.peaks.findPeaks, xrdc.peaks.fitPeak.

    arguments
        scan                      (1,1) struct
        peaks                           struct
        options.HeightFraction    (1,1) double {mustBePositive, ...
                                         mustBeLessThanOrEqual(options.HeightFraction, 1)} = 0.5
        options.Oversample        (1,1) double {mustBeInteger, ...
                                         mustBePositive} = 100
        options.MergeTolerance    (1,1) double {mustBeNonnegative} = 0.03
    end

    if isempty(peaks)
        peaksOut = peaks;
        return
    end

    x = double(scan.twoTheta(:));
    y = double(scan.counts(:));
    n = numel(x);
    if n < 3 || numel(x) ~= numel(y)
        error('xrdc:peaks:badScan', ...
            'scan.twoTheta / scan.counts must be matching vectors of length >= 3.');
    end

    % Build the over-sampled grid. Unlike the Delphi code, we use the
    % actual (possibly non-uniform) 2θ axis rather than index-based
    % interpolation, so this also works on scans with variable step.
    os = options.Oversample;
    mFine = (n - 1) * os + 1;
    idxBase     = (0:mFine-1).' / os;          % fractional index into x
    lIdx        = floor(idxBase);
    frac        = idxBase - lIdx;
    lIdx(end)   = n - 2;                       % clamp the final point
    frac(end)   = 1;
    xFine = x(lIdx + 1) .* (1 - frac) + x(lIdx + 2) .* frac;
    yFine = y(lIdx + 1) .* (1 - frac) + y(lIdx + 2) .* frac;

    peaksOut = peaks;
    for p = 1:numel(peaksOut)
        peakX = peaksOut(p).twoTheta;
        peakH = peaksOut(p).counts;
        if ~isfinite(peakX) || ~isfinite(peakH) || peakH <= 0
            continue
        end

        % Locate the peak on the fine grid (nearest sample).
        [~, ci] = min(abs(xFine - peakX));
        threshold = peakH * options.HeightFraction;

        % Walk left/right from the peak until the counts drop below
        % threshold — matches the Delphi while-loops at lines 1743, 1748.
        li = ci;
        while li > 1 && yFine(li) > threshold
            li = li - 1;
        end
        ri = ci;
        while ri < mFine && yFine(ri) > threshold
            ri = ri + 1;
        end
        if ri <= li
            continue   % degenerate — leave this peak alone
        end

        % Area bisector. The Delphi code does
        %     b := 0; i := l;
        %     while b < a/2 do { i++; b += myCounts[i]; }
        % (i.e. accumulate from l+1 until half the area is reached). We
        % reproduce that behaviour exactly with cumsum — using counts,
        % not counts * Δx, so the result matches the Delphi integer-sum
        % behaviour on uniform grids.
        window  = yFine(li:ri);
        totalA  = sum(window);
        cs      = cumsum(window);
        bisIdx  = find(cs >= totalA / 2, 1, 'first');
        if isempty(bisIdx)
            continue
        end
        % Translate local window index back to fine-grid index; the
        % Delphi loop starts at l and increments before adding, so the
        % index that satisfies the sum lies at l + bisIdx. We clamp to
        % ri to be safe.
        bisFine = min(li + bisIdx, ri);
        newX    = xFine(bisFine);
        newH    = yFine(bisFine);

        peaksOut(p).twoTheta  = newX;
        peaksOut(p).counts    = newH;
        peaksOut(p).fwhm      = xFine(ri) - xFine(li);
        peaksOut(p).leftHalf  = xFine(li);
        peaksOut(p).rightHalf = xFine(ri);
        % Re-map .index to the nearest original-grid sample
        [~, peaksOut(p).index] = min(abs(x - newX));
    end

    % Post-refinement de-dup — xrdc1.pas:1770-1784. Repeat until stable.
    tol = options.MergeTolerance;
    if tol > 0
        changed = true;
        while changed && numel(peaksOut) > 1
            changed = false;
            tt = [peaksOut.twoTheta];
            for i = 1:numel(peaksOut)-1
                if abs(tt(i+1) - tt(i)) < tol
                    % Keep the taller of the two
                    if peaksOut(i).counts >= peaksOut(i+1).counts
                        peaksOut(i+1) = [];
                    else
                        peaksOut(i) = [];
                    end
                    changed = true;
                    break
                end
            end
        end
    end
end
