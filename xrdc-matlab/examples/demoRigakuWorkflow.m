%% Demo: Complete Rigaku workflow — load, find peak, fit, export paper figure.
%  Takes a single .txt export from the Rigaku SmartLab machine and
%  produces a publication-ready PNG. Designed as the starting template
%  when a new sample comes off the machine.

addpath(fileparts(fileparts(mfilename('fullpath'))));

dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '..', '..', 'rexdrctomatlabport_rigakudatasets');
fname = 'TR_S05_PTO_STO(100)_600c_200mT_1000sh_2hz_film RC_04092026.txt';
scan  = xrdc.io.readScan(fullfile(dataDir, fname));

fprintf('Loaded %s\n', scan.identifier);
fprintf('  scanType = %s, %d points, lambda = %.4f A\n', ...
    scan.scanType, numel(scan.twoTheta), scan.lambda);

%% Detect the tallest peak
pk = xrdc.peaks.findPeaks(scan, 'MinProminence', max(scan.counts) * 0.05);
assert(~isempty(pk), 'No peaks detected — check threshold.');
[~, idxMax] = max([pk.counts]);
pkMain = pk(idxMax);

%% Fit a Lorentzian around it (+/-0.5 deg window)
window = [pkMain.twoTheta - 0.5, pkMain.twoTheta + 0.5];
fit    = xrdc.peaks.fitPeak(scan, window, 'Shape', "lorentz");

fprintf('Peak centre : %.4f deg\n',      fit.twoTheta);
fprintf('FWHM        : %.4f deg (%.1f arcsec)\n', fit.fwhm, fit.fwhm*3600);
fprintf('R^2         : %.4f\n',          fit.rSquared);

%% Plot: scan + fit overlay in publication style
h = xrdc.plot.plotScan(scan, ...
    'Title',     sprintf("PTO/STO film RC — FWHM = %.1f arcsec", fit.fwhm*3600), ...
    'LogY',      true, ...
    'ShowPeaks', false);
xlabel(h.ax, '\omega (\circ)');     % RC x-axis is really omega

hold(h.ax, 'on');
plot(h.ax, fit.xFit, fit.yFit, '--', ...
    'Color', [0.85 0.2 0.2], 'LineWidth', 1.5, 'DisplayName', 'Lorentzian fit');
plot(h.ax, pkMain.twoTheta, pkMain.counts, 'o', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.8 0.2], ...
    'MarkerSize', 10, 'DisplayName', 'Peak');
legend(h.ax, 'Location', 'best');
hold(h.ax, 'off');

outPath = fullfile(pwd, 'demoRigakuWorkflow.png');
exportgraphics(h.figure, outPath, 'Resolution', 600);
fprintf('Saved: %s\n', outPath);
