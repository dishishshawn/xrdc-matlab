%% Demo: θ-2θ scan from the Paik lab Rigaku SmartLab (PbTiO3 / SrTiO3).
%  Reproduces the style of Fig 2(a) in Schwaigert et al. JVST A 41, 022703 (2023):
%  log-Y, Arial 18 pt ticks, substrate ticks overlaid.
%
%  Run from the xrdc-matlab root.

addpath(fileparts(fileparts(mfilename('fullpath'))));

dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '..', '..', 'rexdrctomatlabport_rigakudatasets');
defaultFname = 'TR_S04_PTO_STO(100)_750c_200mT_1000sh_3hz_2theta omega_04072026.txt';
if ~exist('fname', 'var') || isempty(fname) ...
        || ~isfile(fullfile(dataDir, fname)) ...
        || ~contains(lower(string(fname)), ["2theta", "th2th", "2th"])
    fname = defaultFname;
end
scan = xrdc.io.readScan(fullfile(dataDir, fname));

fprintf('Loaded: %s\n', scan.identifier);
fprintf('  format = %s,  scanType = %s\n', scan.sourceFormat, scan.scanType);
fprintf('  2θ range = [%.2f, %.2f]°,  %d points\n', ...
    scan.twoTheta(1), scan.twoTheta(end), numel(scan.twoTheta));

%% Detect peaks
pk = xrdc.peaks.findPeaks(scan, ...
    'MinProminence',  100, ...
    'MinHeight',      50);
scan.peaks = pk;
fprintf('\nFound %d peaks:\n', numel(pk));
for i = 1:numel(pk)
    fprintf('  2θ = %7.3f°   I = %9.1f cps   FWHM = %.3f°\n', ...
        pk(i).twoTheta, pk(i).counts, pk(i).fwhm);
end

%% Overlay substrate (STO 00L family) predictions
lat  = struct('system', 'cubic', 'a', 3.905);      % SrTiO3
T    = xrdc.lattice.simulatePattern(lat, [0 0 0 0 1 4], scan.lambda, ...
    'TwoThetaRange', [scan.twoTheta(1), scan.twoTheta(end)]);

%% Plot
h = xrdc.plot.plotScan(scan, 'Title', "PbTiO_3 / SrTiO_3(100) — θ-2θ");
hold(h.ax, 'on');
for i = 1:numel(T.twoTheta)
    xline(h.ax, T.twoTheta(i), '--', char(T.label(i)), ...
        'Color', [0.6 0 0], 'LineWidth', 1.0, 'FontSize', 12, ...
        'LabelVerticalAlignment', 'top');
end

[~, stem, ~] = fileparts(fname);
outPath = fullfile(pwd, sprintf('th2th_%s.png', stem));
exportgraphics(h.figure, outPath, 'Resolution', 600);
fprintf('\nSaved: %s\n', outPath);
