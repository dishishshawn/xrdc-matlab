%% Demo: XRR (specular reflectivity) scan with Kiessig-fringe thickness.
%  Pipeline lives in +xrdc/+xrr/analyzeFringes.m. See that file for the
%  edge-detection / detrend / FFT / fringe-fit details.

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

%% Run the full XRR fringe pipeline
res = xrdc.xrr.analyzeFringes(scan);

fprintf('Critical edge ≈ %.3f°. Found %d Kiessig fringes in [%.2f°, %.2f°].\n', ...
    res.twoThetaC, res.nFringesDetected, res.lowerBound, res.upperBound);

if isnan(res.thicknessFftNm)
    warning('Could not estimate a fringe period — check the scan range.');
elseif res.nFringesDetected < 3
    warning('Only %d fringe(s) detected — reporting FFT thickness only.', ...
        res.nFringesDetected);
    fprintf('FFT thickness  d ≈ %.2f nm (period = %.3f° in 2θ)\n', ...
        res.thicknessFftNm, res.fringePeriodDeg);
else
    fprintf(['Film thickness  d = %.2f ± %.2f nm  (quadratic Kiessig, ' ...
             '%d fringes, λ=%.4f Å)\n'], ...
        res.thicknessQuadNm, res.thicknessQuadSeNm, res.nFringesDetected, scan.lambda);
    fprintf(['  cross-checks:  linear sinθ = %.2f nm,  FFT = %.2f nm  ' ...
             '(fringe period %.3f°, SNR %.1f)\n'], ...
        res.thicknessFitNm, res.thicknessFftNm, res.fringePeriodDeg, res.fringeSnr);
    fprintf('  recovered θ_c   = %.3f° (TER plateau peak was at %.3f°)\n', ...
        res.twoThetaCRecovered, res.twoThetaC);
end

%% Plot
scan.peaks = res.peaks;
if isnan(res.thicknessNm)
    title_str = "XRR — insufficient fringe visibility";
elseif isnan(res.thicknessQuadNm)
    title_str = sprintf("XRR — d ≈ %.1f nm (FFT only)", res.thicknessFftNm);
else
    title_str = sprintf("XRR — d = %.1f nm", res.thicknessQuadNm);
end
h = xrdc.plot.plotScan(scan, ...
    'Title',     title_str, ...
    'ShowPeaks', true);
xlim(h.ax, [0, min(5, scan.twoTheta(end))]);

[~, stem, ~] = fileparts(fname);
outPath = fullfile(pwd, sprintf('xrr_%s.png', stem));
exportgraphics(h.figure, outPath, 'Resolution', 600);
fprintf('Saved: %s\n', outPath);
