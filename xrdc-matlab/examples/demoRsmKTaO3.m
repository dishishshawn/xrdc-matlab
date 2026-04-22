%DEMORSMKTAO3  Reproduce Fig 2(e) of Schwaigert et al. JVST A 41, 022703 (2023).
%
%   The paper's Fig 2(e) shows an RSM around the GdScO3(332) substrate and
%   KTaO3(103) film reflections on a 18 nm KTaO3/GdScO3(110)o sample.
%   Axes: Q_x ∈ [0.245, 0.255] Å⁻¹, Q_z ∈ [0.74, 0.765] Å⁻¹.
%   Log colorbar, 1…10^5 counts.
%
%   Dataset: the three RSM slices captured for the 112 / 103 reflection
%   region live in ../rexdrctomatlabport/.  Each file is a θ-2θ slice at a
%   different ω (three omega offsets: C / L / R).
%
%   Run this script from the xrdc-matlab/ root.

clear;
dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '..', '..', 'rexdrctomatlabport');

pattern = 'HP PtO2 on TiO2 001 112 RSM*.xrdml';
scans = xrdc.rsm.loadAreaScan(dataDir, 'Pattern', pattern);
fprintf('Loaded %d RSM slices.\n', numel(scans));
fprintf('  Omega values: %s\n', mat2str([scans.secondAxis], 4));

% Match Fig 2(e) style: contourf, log decade colorbar, Painters renderer
h = xrdc.plot.plotRsm(scans, ...
    'Mode',       "contourf", ...
    'NContours',  30, ...
    'Imin',       1, ...
    'Imax',       1e5, ...
    'Colormap',   "turbo", ...
    'ExportPath', fullfile(pwd, 'demoRsm_fig2e.png'));

% Tighten to Fig 2(e) axes window if the peak is in range.
% (Exact limits depend on reflection; auto-range by default.)
title(h.ax, 'KTaO_3 (103) / GdScO_3 (332) RSM');

fprintf('Saved: demoRsm_fig2e.png (600 dpi)\n');
