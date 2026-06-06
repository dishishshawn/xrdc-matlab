function result = analyzeFringes(scan, options)
%ANALYZEFRINGES  End-to-end Kiessig-fringe analysis for an XRR scan.
%   result = xrdc.xrr.analyzeFringes(scan)
%   result = xrdc.xrr.analyzeFringes(scan, Name=Value, ...)
%
%   Pipeline (drop-in replacement for the inline code in demoXRR.m and
%   xrdcApp's runXRR):
%
%     1. Critical edge: direct-beam dip → TER plateau peak
%        (xrdc.xrr.findCriticalEdge). The earlier "steepest descent of
%        log(counts) below 1.5°" heuristic locked onto the knife-edge
%        slope at ~0.04°, not the actual XRR edge, which gave ~3× thickness
%        errors on real Rigaku SmartLab data.
%
%     2. Polynomial detrend (xrdc.xrr.detrendLog). Order-4 polynomial fit
%        to log10(counts) on [θ_c+EdgeBuffer, ScanEnd]. Replaces the
%        window-based sgolayfilt envelope, which absorbed fringes whenever
%        the envelope span approached the fringe period.
%
%     3. Dominant fringe period from FFT (xrdc.xrr.dominantPeriod).
%        Robust against individual peak quality; provides a frequency-domain
%        thickness as the primary estimate.
%
%     4. Peak detection seeded by the FFT period: MinSeparation =
%        SeparationFactor × period, so the peak detector cannot pick up
%        sub-fringe residuals.
%
%     5. Linear fit of sin(θ) vs fringe index → thicknessFromFringes
%        (Curve Fitting path with confidence interval when available).
%
%   Input
%     scan : scan struct with .twoTheta, .counts, .lambda
%
%   Name/Value options
%     'LowerBound'       (default NaN)
%         Manual override for the lower 2θ analysis bound. NaN → auto
%         (θ_c + EdgeBuffer).
%     'UpperBound'       (default 5.0)
%         Upper 2θ bound of the analysis window (degrees).
%     'EdgeBuffer'       (default 0.05)
%         Buffer past the critical edge before the analysis window starts.
%     'PolyOrder'        (default 4)
%         Polynomial order for log-counts detrend.
%     'MinProminence'    (default 0.02)
%         Prominence threshold (in log10-counts units) for findPeaks.
%     'SeparationFactor' (default 0.7)
%         MinSeparation = factor × FFT-period (degrees). 0.7 rejects
%         sub-fringe artefacts without missing real fringes.
%     'EdgeOptions'      (default struct())
%         Forwarded to xrdc.xrr.findCriticalEdge as Name/Value.
%
%   Output (struct)
%     .twoThetaC          Critical edge (degrees)
%     .lowerBound         Actual analysis lower bound used
%     .upperBound         Actual analysis upper bound used
%     .thicknessNm        Best-estimate thickness — quadratic Kiessig
%                          when ≥3 fringes detected, else FFT, else NaN
%     .thicknessFftNm     Thickness from FFT dominant period (nm)
%     .thicknessFitNm     Thickness from linear sin θ fit (nm); NaN
%                          if fewer than 2 fringes were detected
%     .thicknessQuadNm    Thickness from the full quadratic Kiessig
%                          relation sin²θ = sin²θ_c + (nλ/2t)²;
%                          NaN with fewer than 3 fringes
%     .thicknessFitSeNm   1-σ SE on the linear fit thickness
%     .thicknessCi95Nm    [lo hi] 95% CI on the quadratic thickness
%                          (NaN when polyfit fallback)
%     .twoThetaCRecovered Critical edge recovered by the quadratic fit
%                          (degrees) — usually a few tenths of a degree
%                          below the TER plateau peak
%     .fringePeriodDeg    Dominant fringe period (2θ degrees)
%     .fringeSnr          FFT peak / median(band)
%     .nFringesDetected   Number of fringes returned by findPeaks
%     .peaks              Fringe struct array (xrdc.peaks.findPeaks shape)
%                          with .counts re-attached to the raw curve
%     .xDetrend, .yDetrend  Detrended trace inside the analysis window
%     .edgeInfo           Full struct from xrdc.xrr.findCriticalEdge
%     .periodInfo         Full struct from xrdc.xrr.dominantPeriod
%     .thicknessInfo      Full struct from xrdc.lattice.thicknessFromFringes

    arguments
        scan                            (1,1) struct
        options.LowerBound              (1,1) double = NaN
        options.UpperBound              (1,1) double {mustBePositive} = 5.0
        options.EdgeBuffer              (1,1) double {mustBeNonnegative} = 0.05
        options.PolyOrder               (1,1) double {mustBeInteger, mustBePositive} = 4
        options.MinProminence           (1,1) double {mustBeNonnegative} = 0.02
        options.SeparationFactor        (1,1) double {mustBePositive} = 0.7
        options.EdgeOptions             (1,1) struct = struct()
    end

    if ~isfield(scan, 'twoTheta') || ~isfield(scan, 'counts') || ~isfield(scan, 'lambda')
        error('xrdc:xrr:badScan', ...
            'scan must have .twoTheta, .counts and .lambda.');
    end

    tt  = double(scan.twoTheta(:));
    cnt = double(scan.counts(:));
    if numel(tt) ~= numel(cnt)
        error('xrdc:xrr:sizeMismatch', ...
            'scan.twoTheta and scan.counts must have the same length.');
    end

    % --- 1. Critical edge ---
    edgeArgs = optionsToNVPairs(options.EdgeOptions);
    edgeInfo = xrdc.xrr.findCriticalEdge(tt, cnt, edgeArgs{:});
    twoThetaC = edgeInfo.twoThetaC;

    % --- 2. Analysis window + detrend ---
    if isnan(options.LowerBound)
        lo = twoThetaC + options.EdgeBuffer;
    else
        lo = options.LowerBound;
    end
    hi = min(options.UpperBound, tt(end));
    mask = tt >= lo & tt <= hi;
    if nnz(mask) < max(options.PolyOrder + 5, 20)
        result = buildEmptyResult(twoThetaC, lo, hi, edgeInfo);
        return
    end
    xs   = tt(mask);
    ys   = log10(max(cnt(mask), 1));
    yDet = xrdc.xrr.detrendLog(xs, ys, 'Order', options.PolyOrder);

    % --- 3. Dominant FFT period ---
    periodInfo  = xrdc.xrr.dominantPeriod(xs, yDet);
    periodDeg   = periodInfo.periodDeg;
    lambdaNm    = scan.lambda / 10;
    thetaCenter = mean(xs) / 2 * pi / 180;   % radians, half of mid-2θ
    sinPeriod   = (periodDeg * pi / 180) / 2 * cos(thetaCenter);
    tFftNm      = lambdaNm / (2 * sinPeriod);

    % --- 4. Peak detection seeded by FFT period ---
    subScan = scan;
    subScan.twoTheta = xs;
    subScan.counts   = yDet;
    minSep  = options.SeparationFactor * periodDeg;
    pk = xrdc.peaks.findPeaks(subScan, ...
        'MinProminence', options.MinProminence, ...
        'MinSeparation', minSep);

    % Re-attach raw counts to fringe markers (so plots land on the real curve).
    if ~isempty(pk)
        for k = 1:numel(pk)
            [~, j] = min(abs(tt - pk(k).twoTheta));
            pk(k).counts = cnt(j);
            pk(k).index  = j;
        end
    end

    % --- 5. Linear and quadratic Kiessig fits ---
    if numel(pk) >= 2
        thickInfo = xrdc.lattice.thicknessFromFringes( ...
            [pk.twoTheta]', scan.lambda);
        tFitNm   = thickInfo.thicknessFitNm;
        tFitSe   = thickInfo.thicknessFitSeNm;
    else
        thickInfo = struct();
        tFitNm = NaN; tFitSe = NaN;
    end

    if numel(pk) >= 3
        quadInfo = xrdc.xrr.thicknessQuadratic( ...
            [pk.twoTheta]', scan.lambda);
        tQuadNm   = quadInfo.thicknessNm;
        tQuadSe   = quadInfo.thicknessSeNm;
        tQuadCi   = quadInfo.thicknessCi95Nm;
        twoThetaCRecovered = quadInfo.twoThetaC;
    else
        quadInfo = struct();
        tQuadNm = NaN; tQuadSe = NaN; tQuadCi = [NaN NaN];
        twoThetaCRecovered = NaN;
    end

    % Best-estimate thickness: quadratic when available, FFT otherwise.
    if ~isnan(tQuadNm)
        tNm = tQuadNm;
    else
        tNm = tFftNm;
    end

    result = struct( ...
        'twoThetaC',           twoThetaC, ...
        'twoThetaCRecovered',  twoThetaCRecovered, ...
        'lowerBound',          lo, ...
        'upperBound',          hi, ...
        'thicknessNm',         tNm, ...
        'thicknessFftNm',      tFftNm, ...
        'thicknessFitNm',      tFitNm, ...
        'thicknessQuadNm',     tQuadNm, ...
        'thicknessFitSeNm',    tFitSe, ...
        'thicknessQuadSeNm',   tQuadSe, ...
        'thicknessCi95Nm',     tQuadCi, ...
        'fringePeriodDeg',     periodDeg, ...
        'fringeSnr',           periodInfo.signalToNoise, ...
        'nFringesDetected',    numel(pk), ...
        'peaks',               pk, ...
        'xDetrend',            xs, ...
        'yDetrend',            yDet, ...
        'edgeInfo',            edgeInfo, ...
        'periodInfo',          periodInfo, ...
        'thicknessInfo',       thickInfo, ...
        'quadInfo',            quadInfo);
end

% =========================================================================

function r = buildEmptyResult(twoThetaC, lo, hi, edgeInfo)
    r = struct( ...
        'twoThetaC',           twoThetaC, ...
        'twoThetaCRecovered',  NaN, ...
        'lowerBound',          lo, ...
        'upperBound',          hi, ...
        'thicknessNm',         NaN, ...
        'thicknessFftNm',      NaN, ...
        'thicknessFitNm',      NaN, ...
        'thicknessQuadNm',     NaN, ...
        'thicknessFitSeNm',    NaN, ...
        'thicknessQuadSeNm',   NaN, ...
        'thicknessCi95Nm',     [NaN NaN], ...
        'fringePeriodDeg',     NaN, ...
        'fringeSnr',           NaN, ...
        'nFringesDetected',    0, ...
        'peaks',               struct([]), ...
        'xDetrend',            [], ...
        'yDetrend',            [], ...
        'edgeInfo',            edgeInfo, ...
        'periodInfo',          struct(), ...
        'thicknessInfo',       struct(), ...
        'quadInfo',            struct());
end

function nv = optionsToNVPairs(s)
    f = fieldnames(s);
    nv = cell(1, 2*numel(f));
    for i = 1:numel(f)
        nv{2*i-1} = f{i};
        nv{2*i}   = s.(f{i});
    end
end
