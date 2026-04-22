%% Demo: Rocking curve — Lorentzian fit, report FWHM in arcsec.
%  Fig 2(c)/(d) style in Schwaigert et al. JVST A 41, 022703 (2023).
%  Rigaku RC files label the x-axis "2θ" but the data are actually ω
%  (detuning from the Bragg peak); readRigakuTxt sets scanType = "omega"
%  from the filename.

addpath(fileparts(fileparts(mfilename('fullpath'))));

dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '..', '..', 'rexdrctomatlabport_rigakudatasets');
defaultFname = 'TR_S11_PTO_STO(100)_580c_150mT_20000sh_5hz_Film RC_04162026.txt';
% Use workspace fname only if it exists AND looks like a rocking curve. Keywords
% cover Rigaku (_RC_, "RC") and PANalytical ("rocking") naming conventions.
if ~exist('fname', 'var') || isempty(fname) ...
        || ~isfile(fullfile(dataDir, fname)) ...
        || ~contains(lower(string(fname)), ["rc", "rocking"])
    fname = defaultFname;
end
scan = xrdc.io.readScan(fullfile(dataDir, fname));
fprintf('Loaded %s  (scanType = %s)\n', scan.identifier, scan.scanType);

%% Find the main peak (largest prominence)
pk = xrdc.peaks.findPeaks(scan, ...
    'MinProminence',  max(scan.counts) * 0.1);
assert(~isempty(pk), 'No peak found in the RC.');
% Pick the tallest
[~, idxMax] = max([pk.counts]);
pkMain = pk(idxMax);

%% Fit Lorentzian in a ±0.5° window
w      = 0.5;
window = [pkMain.twoTheta - w, pkMain.twoTheta + w];
fit    = xrdc.peaks.fitPeak(scan, window, 'Shape', "lorentz");

fwhm_deg    = fit.fwhm;
fwhm_arcsec = fwhm_deg * 3600;
fprintf('Peak centre ω₀ = %.4f°\n', fit.twoTheta);
fprintf('FWHM         = %.4f°  =  %.1f arcsec\n', fwhm_deg, fwhm_arcsec);

%% Plot
h = xrdc.plot.plotScan(scan, ...
    'Title',     sprintf("Film RC — FWHM = %.1f arcsec", fwhm_arcsec), ...
    'LogY',      true, ...
    'ShowPeaks', false);
xlabel(h.ax, '\omega (\circ)');    % correct label — x-axis is omega

% Overlay the Lorentzian fit (provided in fit.xFit, fit.yFit by fitPeak)
hold(h.ax, 'on');
plot(h.ax, fit.xFit, fit.yFit, '--', ...
    'Color', [0.85 0.2 0.2], 'LineWidth', 1.5);

[~, stem, ~] = fileparts(fname);
outPath = fullfile(pwd, sprintf('RC_%s.png', stem));
exportgraphics(h.figure, outPath, 'Resolution', 600);
fprintf('Saved: %s\n', outPath);
