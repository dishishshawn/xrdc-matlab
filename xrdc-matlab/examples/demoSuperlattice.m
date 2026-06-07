%% Demo: superlattice period Λ from θ-2θ satellite peaks.
%  A superlattice [A/B]×N shows satellite reflections evenly spaced in sin θ
%  around the average Bragg peak. Their spacing gives the bilayer period via
%  xrdc.lattice.superlatticePeriod (same math as Kiessig/Laue fringes).
%
%  This demo SYNTHESIZES a pattern with a known period so it runs without a
%  data file. Point it at a real θ-2θ superlattice scan by replacing the
%  "Synthesize" block with xrdc.io.readScan(...) and choosing the window
%  around the main peak.

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% --- Synthesize a superlattice θ-2θ pattern with a known period ---
lambda     = 1.5406;                 % Cu Kα, Å
truePeriod = 15;                     % nm  (the answer we should recover)
lambda_nm  = lambda / 10;
slope      = lambda_nm / (2 * truePeriod);    % Δ(sin θ) per satellite order

sinTh0 = sind(46/2);                 % main (0-order) peak near 2θ = 46°
orders = (-3:3).';                   % satellite orders -3 … +3
sinSat = sinTh0 + orders * slope;
twoThetaSatTrue = 2 * asind(sinSat);
amps   = 1e6 * 0.45 .^ abs(orders);  % intensity falls off with |order|

tt = (40 : 0.004 : 52).';
counts = 30 * ones(size(tt));        % flat background
sigma  = 0.03;                       % peak σ in 2θ (FWHM ≈ 0.07°)
for k = 1:numel(orders)
    counts = counts + amps(k) * exp(-(tt - twoThetaSatTrue(k)).^2 / (2*sigma^2));
end

scan = xrdc.io.emptyScan();
scan.twoTheta   = tt;
scan.counts     = counts;
scan.lambda     = lambda;
scan.scanType   = "twoThetaOmega";
scan.identifier = "synthetic SL [A/B] (Λ=15 nm)";

%% --- Detect the satellites, compute the period ---
pk = xrdc.peaks.findPeaks(scan, ...
    'MinProminence', max(scan.counts) * 0.01, ...
    'MinSeparation', 0.2, ...
    'TwoThetaRange', [42, 50]);
assert(numel(pk) >= 2, 'Need ≥2 satellites; found %d.', numel(pk));

twoThetaSat = sort([pk.twoTheta].');
res = xrdc.lattice.superlatticePeriod(twoThetaSat, scan.lambda);

fprintf('Detected %d satellites in [%.2f, %.2f]°\n', ...
    res.nSatellites, twoThetaSat(1), twoThetaSat(end));
fprintf('Superlattice period  Λ = %.2f nm  (fit %.2f ± %.2f nm, R²=%.4f)\n', ...
    res.periodNm, res.periodFitNm, res.periodFitSeNm, res.rSquared);
fprintf('True period was %.2f nm → recovery error %.2f%%\n', ...
    truePeriod, 100*(res.periodFitNm - truePeriod)/truePeriod);

%% --- Plot ---
scan.peaks = pk;
h = xrdc.plot.plotScan(scan, ...
    'Title',     sprintf("Superlattice — Λ = %.1f nm (%d satellites)", ...
                         res.periodFitNm, res.nSatellites), ...
    'LogY',      true, ...
    'ShowPeaks', true);
xlim(h.ax, [42, 50]);

outPath = fullfile(pwd, 'superlattice_demo.png');
exportgraphics(h.figure, outPath, 'Resolution', 600);
fprintf('Saved: %s\n', outPath);
