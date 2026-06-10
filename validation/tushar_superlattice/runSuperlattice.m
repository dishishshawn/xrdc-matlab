%% runSuperlattice.m — test xrdc.lattice.superlatticePeriod on Tushar's SP79 SL.
%  Sample: SP79 WO3/LAO superlattice on LAO(100), nominally 3 uc x 26 stacks,
%          ~60 nm total (700 C, 50 mT). Cu K-alpha, theta-2theta scan.
%
%  Two parts:
%    A) Known-answer sanity check on a synthetic SL pattern (period = 15 nm).
%    B) The real SP79 scan: detect the film-peak satellite/fringe series around
%       the (001) and (002) reflections and run superlatticePeriod on them.
%
%  Run:  matlab -batch "run('validation/tushar_superlattice/runSuperlattice.m')"

here    = fileparts(mfilename('fullpath'));
pkgRoot = fullfile(here, '..', '..', 'xrdc-matlab');
dataDir = fullfile(here, '..', '..', 'data', 'tushar superlattice data');
addpath(pkgRoot);

lambda = 1.5406;   % Cu K-alpha, Angstrom

%% ===== A. Known-answer sanity check (synthetic, period = 15 nm) =====
fprintf('\n===== A. Synthetic known-answer check =====\n');
truePeriod = 15;                       % nm
lambda_nm  = lambda / 10;
slope      = lambda_nm / (2*truePeriod);     % d(sin theta) per order
sinTh0     = sind(46/2);
orders     = (-3:3).';
twoThetaSatTrue = 2 * asind(sinTh0 + orders*slope);
sA = xrdc.lattice.superlatticePeriod(sort(twoThetaSatTrue), lambda);
fprintf('  true = %.2f nm | recovered fit = %.3f +/- %.3f nm, R^2 = %.6f\n', ...
    truePeriod, sA.periodFitNm, sA.periodFitSeNm, sA.rSquared);
fprintf('  recovery error = %.3f%%\n', 100*(sA.periodFitNm - truePeriod)/truePeriod);

%% ===== B. Real SP79 superlattice scan =====
fprintf('\n===== B. SP79 WO3/LAO superlattice (real scan) =====\n');
fSL = fullfile(dataDir, ...
    'SP79_WO3_LAO SL_on_LAO_100_3_uc_26_stacks_700C_50mT_60nm_2theta_omega_04232026.txt');
assert(isfile(fSL), 'Missing data file: %s', fSL);

scan = xrdc.io.readScan(fSL);
fprintf('Loaded "%s"\n  scanType=%s, %d pts, 2theta in [%.2f, %.2f] deg\n', ...
    scan.identifier, scan.scanType, numel(scan.twoTheta), ...
    scan.twoTheta(1), scan.twoTheta(end));

% Analyse both accessible orders. The window stops short of the substrate
% Bragg spike so findPeaks locks onto film satellites, not the substrate.
windows = struct( ...
    'name',   {'LAO(001)',     'LAO(002)'}, ...
    'range',  {[21.0, 23.45],  [45.6, 47.90]}, ...
    'subPeak',{23.46,          47.96});

for w = windows
    fprintf('\n--- %s film region  %.2f-%.2f deg ---\n', w.name, w.range(1), w.range(2));

    % Detect satellite maxima. Prominence on the log-intensity envelope is the
    % right knob for the gentle SL humps; use a generous separation so the fine
    % total-thickness Laue fringes are not each picked as a "satellite".
    pk = xrdc.peaks.findPeaks(scan, ...
        'TwoThetaRange',  w.range, ...
        'MinProminence',  max(scan.counts)*2e-4, ...
        'MinSeparation',  0.12);

    if numel(pk) < 2
        fprintf('  Only %d peak(s) detected -> cannot form a satellite series.\n', numel(pk));
        continue;
    end

    twoThetaSat = sort([pk.twoTheta].');
    fprintf('  %d candidate satellites: %s\n', numel(twoThetaSat), ...
        strjoin(compose('%.3f', twoThetaSat), ', '));

    res = xrdc.lattice.superlatticePeriod(twoThetaSat, lambda);
    fprintf('  Period Lambda = %.2f nm (endpoint) | fit %.2f +/- %.2f nm | R^2 = %.4f | %s\n', ...
        res.periodNm, res.periodFitNm, res.periodFitSeNm, res.rSquared, res.fitBackend);
    fprintf('  max |residual| = %.4f deg over %d satellites\n', ...
        max(abs(res.residualsDeg)), res.nSatellites);
end

%% ===== Plot the scan with the (001) window annotated =====
h = xrdc.plot.plotScan(scan, 'Title', ...
    "SP79 WO3/LAO SL on LAO(100) — \theta-2\theta", 'LogY', true);
xlim(h.ax, [20.5, 24.0]);
outPng = fullfile(here, 'SP79_superlattice_001.png');
exportgraphics(h.figure, outPng, 'Resolution', 300);
fprintf('\nSaved %s\n', outPng);

fprintf('\nDone.\n');
