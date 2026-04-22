%% Demo: XRR (specular reflectivity) scan with Kiessig-fringe thickness.
%  Fringe-index → polyfit of sin(θ_n) vs n → slope → film thickness.

addpath(fileparts(fileparts(mfilename('fullpath'))));

dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '..', '..', 'rexdrctomatlabport_rigakudatasets');
if ~exist('fname', 'var') || isempty(fname)
    fname = 'TR_S10_PTO_STO(100)_500c_150mT_20000sh_5hz_XRR_04162026.txt';
end
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
