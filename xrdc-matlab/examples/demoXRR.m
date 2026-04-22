%% Demo: XRR (specular reflectivity) scan with Kiessig-fringe thickness.
%  Fringe-index → polyfit of sin(θ_n) vs n → slope → film thickness.

addpath(fileparts(fileparts(mfilename('fullpath'))));

dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '..', '..', 'rexdrctomatlabport_rigakudatasets');
fname = 'TR_S10_PTO_STO(100)_500c_150mT_20000sh_5hz_XRR_04162026.txt';
scan = xrdc.io.readScan(fullfile(dataDir, fname));
fprintf('Loaded %s: %d points, 2θ ∈ [%.3f, %.3f]°\n', ...
    scan.identifier, numel(scan.twoTheta), scan.twoTheta(1), scan.twoTheta(end));

%% Find Kiessig fringes (small 2θ, after critical edge)
% Work on a sub-scan in the fringe region
subScan = scan;
mask    = scan.twoTheta > 0.5 & scan.twoTheta < 3.0;
subScan.twoTheta = scan.twoTheta(mask);
subScan.counts   = scan.counts(mask);

pk = xrdc.peaks.findPeaks(subScan, ...
    'MinProminence',  max(subScan.counts) * 0.02, ...
    'MinSeparation',  0.05);
fprintf('Found %d Kiessig fringes between 0.5° and 3°.\n', numel(pk));

if numel(pk) < 3
    warning('Too few fringes — thickness estimate will be noisy.');
end

%% Thickness from fringe periodicity
ttPk   = [pk.twoTheta];     % fringe 2θ in degrees (thicknessFromFringes wants 2θ)
thick  = xrdc.lattice.thicknessFromFringes(ttPk(:), scan.lambda);
t_nm   = thick.thicknessFitNm;
fprintf('Film thickness  d = %.2f ± %.2f nm (fit of %d fringes, λ=%.4f Å)\n', ...
    t_nm, thick.thicknessFitSeNm, numel(pk), scan.lambda);

%% Plot
scan.peaks = pk;
h = xrdc.plot.plotScan(scan, ...
    'Title',     sprintf("XRR — d = %.1f nm", t_nm), ...
    'ShowPeaks', true);
xlim(h.ax, [0, min(5, scan.twoTheta(end))]);

outPath = fullfile(pwd, 'demoXRR.png');
exportgraphics(h.figure, outPath, 'Resolution', 600);
fprintf('Saved: %s\n', outPath);
