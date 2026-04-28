%% Demo: XRR (specular reflectivity) scan with Kiessig-fringe thickness.
%  Fringe-index → polyfit of sin(θ_n) vs n → slope → film thickness.

addpath(fileparts(fileparts(mfilename('fullpath'))));

dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '..', '..', 'rexdrctomatlabport_rigakudatasets');
defaultFname = 'TR_S10_PTO_STO(100)_500c_150mT_20000sh_5hz_XRR_04162026.txt';
% Use workspace fname only if it exists AND looks like an XRR scan (by filename
% keyword). Prevents stale fnames from previous demos being reused on the wrong
% scan type.
if ~exist('fname', 'var') || isempty(fname) ...
        || ~isfile(fullfile(dataDir, fname)) ...
        || ~contains(lower(string(fname)), "xrr")
    fname = defaultFname;
end
scan = xrdc.io.readScan(fullfile(dataDir, fname));
fprintf('Loaded %s: %d points, 2θ ∈ [%.3f, %.3f]°\n', ...
    scan.identifier, numel(scan.twoTheta), scan.twoTheta(1), scan.twoTheta(end));

%% Find Kiessig fringes (small 2θ, after critical edge)
% XRR intensity drops 5–6 decades over the fringe region, so any prominence
% threshold tied to max(counts) only ever sees the first 1–2 fringes near
% the critical edge. Detect on log10(counts) with the slow envelope removed,
% so each fringe contributes the same ~few-percent ripple regardless of
% absolute height. Critical-edge position is auto-found from the steepest
% descent of log(counts) below 1.5° — removes the need to hand-tune the
% lower bound per sample.
tt    = scan.twoTheta(:);
ylog  = log10(max(double(scan.counts(:)), 1));

dy = gradient(ylog, tt);
edgeRegion = tt < 1.5;
dyMasked = dy;
dyMasked(~edgeRegion) = +Inf;     % force argmin into the edge region
[minSlope, iEdge] = min(dyMasked);
if ~isfinite(minSlope) || minSlope >= 0
    theta_c = 0.4;                % fallback if no clear edge below 1.5°
else
    theta_c = tt(iEdge);
end

winLo = theta_c + 0.05;
winHi = min(5.0, tt(end));
mask  = tt > winLo & tt < winHi;
xseg  = tt(mask);
yseg  = ylog(mask);

% Detrend with a moving-average envelope. Span ~0.4° is wide enough to be
% smooth across several fringes for typical 5–50 nm films, narrow enough to
% follow the steep post-edge decay.
step    = median(diff(xseg));
spanPts = max(5, round(0.4 / max(step, eps)));
envelope = movmean(yseg, spanPts);
ydet     = yseg - envelope;

subScan = scan;
subScan.twoTheta = xseg;
subScan.counts   = ydet;     % detrended log signal, dimensionless

pk = xrdc.peaks.findPeaks(subScan, ...
    'MinProminence', 0.015, ...   % ~3.5% intensity ripple in log units
    'MinSeparation', 0.05);

% Re-attach real counts at each fringe so plotScan markers land on the
% physical curve rather than on the detrended residual.
if ~isempty(pk)
    for k = 1:numel(pk)
        [~, j] = min(abs(scan.twoTheta - pk(k).twoTheta));
        pk(k).counts = scan.counts(j);
        pk(k).index  = j;
    end
end
fprintf('Critical edge ≈ %.3f°. Found %d Kiessig fringes in [%.2f°, %.2f°].\n', ...
    theta_c, numel(pk), winLo, winHi);

%% Thickness from fringe periodicity
if numel(pk) < 2
    warning('Too few fringes (need ≥2) — skipping thickness estimate.');
    t_nm = NaN;
    thick = struct();
else
    if numel(pk) < 3
        warning('Only %d fringe(s) — thickness estimate will be noisy.', numel(pk));
    end
    ttPk   = [pk.twoTheta];     % fringe 2θ in degrees (thicknessFromFringes wants 2θ)
    thick  = xrdc.lattice.thicknessFromFringes(ttPk(:), scan.lambda);
    t_nm   = thick.thicknessFitNm;
    fprintf('Film thickness  d = %.2f ± %.2f nm (fit of %d fringes, λ=%.4f Å)\n', ...
        t_nm, thick.thicknessFitSeNm, numel(pk), scan.lambda);
end

%% Plot
scan.peaks = pk;
if isnan(t_nm)
    title_str = "XRR — insufficient fringe visibility";
else
    title_str = sprintf("XRR — d = %.1f nm", t_nm);
end
h = xrdc.plot.plotScan(scan, ...
    'Title',     title_str, ...
    'ShowPeaks', true);
xlim(h.ax, [0, min(5, scan.twoTheta(end))]);

[~, stem, ~] = fileparts(fname);
outPath = fullfile(pwd, sprintf('xrr_%s.png', stem));
exportgraphics(h.figure, outPath, 'Resolution', 600);
fprintf('Saved: %s\n', outPath);
