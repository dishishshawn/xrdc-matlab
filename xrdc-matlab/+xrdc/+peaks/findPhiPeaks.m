function [peaks, info] = findPhiPeaks(scan, options)
%FINDPHIPEAKS  Detect poles in a phi (azimuthal) scan, robustly.
%   [peaks, info] = xrdc.peaks.findPhiPeaks(scan)
%   [peaks, info] = xrdc.peaks.findPhiPeaks(scan, Name, Value, ...)
%
%   A phi scan of an off-axis reflection shows one sharp pole per
%   symmetry-equivalent in-plane direction (e.g. 4 poles 90 deg apart for
%   4-fold epitaxy). Two things make plain prominence detection fail here,
%   both seen on real lab data:
%
%     1. Weak films. When the tallest pole is only ~20 counts, a
%        "prominence = 10% of max" rule sits near ~2 counts, BELOW the
%        Poisson noise, so single-count noise spikes get counted as poles
%        (a 4-pole scan reported 9). This routine sets the threshold from
%        the NOISE level instead: height >= background + NoiseSigmas*sigma,
%        with a robust median/MAD background and sigma, so it adapts to
%        strong and weak scans alike.
%
%     2. The 360-deg wrap. phi is periodic, so a pole can sit at the scan
%        boundary or appear at both ends. findpeaks cannot flag a maximum
%        at the first/last sample, so the trace is circularly extended
%        before detection, and poles one Period apart (the same crystal
%        direction) are de-duplicated for the symmetry count.
%
%   Input
%     scan : scan struct with .twoTheta (= phi, degrees) and .counts.
%
%   Name/Value
%     'NoiseSigmas'         (default 6)   height threshold = bg + k*sigma
%     'MinProminenceSigmas' (default 4)   prominence threshold = k*sigma
%     'MinSeparation'       (default 30)  min degrees between poles
%     'Period'              (default 360) azimuthal period for wrap/dedup
%     'Margin'              (default 60)  degrees of circular pad each side
%
%   Output
%     peaks : struct array (xrdc.peaks.findPeaks shape) for the in-range
%             poles, sorted by phi.
%     info  : struct with
%               .background, .noiseSigma   robust bg / noise estimate
%               .heightThreshold           bg + NoiseSigmas*sigma
%               .nPoles                    numel(peaks)
%               .nUnique                   poles after 360-deg de-dup
%               .spacings                  deg between unique poles (sorted)
%               .fold                      round(Period/median(spacing)),
%                                          NaN with < 2 unique poles
%
%   See also XRDC.PEAKS.FINDPEAKS.

    arguments
        scan                            (1,1) struct
        options.NoiseSigmas             (1,1) double {mustBePositive}    = 6
        options.MinProminenceSigmas     (1,1) double {mustBeNonnegative} = 4
        options.MinSeparation           (1,1) double {mustBePositive}    = 30
        options.Period                  (1,1) double {mustBePositive}    = 360
        options.Margin                  (1,1) double {mustBeNonnegative} = 60
    end

    if ~isfield(scan, 'twoTheta') || ~isfield(scan, 'counts')
        error('xrdc:peaks:badScan', 'scan must have .twoTheta and .counts.');
    end

    [x, ord] = sort(double(scan.twoTheta(:)));
    y = double(scan.counts(ord));

    % Robust background + noise (median / MAD). Fall back to a Poisson-ish
    % estimate when the trace is flat enough that the MAD collapses to 0.
    bg    = median(y);
    sigma = 1.4826 * median(abs(y - bg));
    if sigma < eps
        sigma = max(sqrt(max(bg, 1)), 1);
    end
    H = bg + options.NoiseSigmas * sigma;
    P = options.MinProminenceSigmas * sigma;

    % Circular extension so a pole at the scan boundary is detectable.
    m   = options.Margin;
    iL  = x > x(end) - m;
    iH  = x < x(1)   + m;
    extScan          = scan;
    extScan.twoTheta = [x(iL) - options.Period; x; x(iH) + options.Period];
    extScan.counts   = [y(iL);                  y; y(iH)];

    pk = xrdc.peaks.findPeaks(extScan, ...
        'MinHeight',     H, ...
        'MinProminence', P, ...
        'MinSeparation', options.MinSeparation);

    % Keep only poles whose centre lies within the original phi range.
    if ~isempty(pk)
        inRange = [pk.twoTheta] >= x(1) - 1e-6 & [pk.twoTheta] <= x(end) + 1e-6;
        pk = pk(inRange);
    end
    if ~isempty(pk)
        [~, si] = sort([pk.twoTheta]);
        pk = pk(si);
    end
    peaks = pk;

    info = struct('background', bg, 'noiseSigma', sigma, ...
        'heightThreshold', H, 'nPoles', numel(pk), ...
        'nUnique', numel(pk), 'spacings', [], 'fold', NaN);
    if isempty(pk)
        return
    end

    % De-duplicate poles one Period apart (same crystal direction) so the
    % symmetry count reflects unique poles, not boundary repeats.
    canon = mod([pk.twoTheta], options.Period);
    uniq  = [];
    for v = canon
        if isempty(uniq)
            uniq(end+1) = v; %#ok<AGROW>
        else
            d = abs(v - uniq);
            d = min(d, options.Period - d);   % circular distance
            if all(d > options.MinSeparation)
                uniq(end+1) = v; %#ok<AGROW>
            end
        end
    end
    uniq = sort(uniq);
    info.nUnique = numel(uniq);
    if numel(uniq) >= 2
        info.spacings = diff(uniq);
        info.fold     = round(options.Period / median(info.spacings));
    end
end
