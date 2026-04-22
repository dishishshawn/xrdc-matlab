%DEMORSMKTAO3  Reciprocal-space map in the style of Schwaigert et al. JVST A 41, 022703 (2023).
%
%   The reference paper's Fig 2(e) shows an RSM on an 18 nm KTaO3 / GdScO3(110)o
%   heterostructure around the (103)/(332) reflections, plotted as a filled
%   log-intensity contour with a decade-tick colorbar (1 … 10^5 counts).
%
%   Dataset: the 112 area scan captured on a PtO2/TiO2(001) sample lives in
%   ../rexdrctomatlabport/ as a single multi-scan XRDML file. The reader
%   xrdc.io.readXrdmlArea expands the ~1500 embedded ω-slices automatically.
%
%   Run from xrdc-matlab/ root.

addpath(fileparts(fileparts(mfilename('fullpath'))));

dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '..', '..', 'rexdrctomatlabport');
if ~exist('fname', 'var') || isempty(fname)
    fname = 'HP PtO2 on TiO2 001 112 RSM_C_HP PtO2 on TiO2 001 112 RSM_C.xrdml';
end

scans = xrdc.rsm.loadAreaScan({fullfile(dataDir, fname)});
fprintf('Loaded %d ω-slices from %s\n', numel(scans), fname);
fprintf('  ω ∈ [%.3f°, %.3f°]\n', ...
    min([scans.secondAxis]), max([scans.secondAxis]));
fprintf('  2θ ∈ [%.3f°, %.3f°], %d points per slice\n', ...
    scans(1).twoTheta(1), scans(1).twoTheta(end), numel(scans(1).twoTheta));

% Fig 2(e) style: contourf, log decade colorbar, Painters renderer.
[~, stem, ~] = fileparts(fname);
outPath = fullfile(pwd, sprintf('rsm_%s.png', stem));

h = xrdc.plot.plotRsm(scans, ...
    'Mode',       "contourf", ...
    'NContours',  40, ...
    'Imin',       1, ...
    'Imax',       1e5, ...
    'Colormap',   "turbo", ...
    'ExportPath', outPath);

title(h.ax, 'PtO_2 (112) RSM — Fig 2(e) style');
fprintf('Saved: %s (600 dpi)\n', outPath);
