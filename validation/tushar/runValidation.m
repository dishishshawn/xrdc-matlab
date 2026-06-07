%% runValidation.m — reproduce Tushar's S25 figures from the app.
%  Generates three PNGs in output/ to compare against reference/:
%    1. S25_XRR_xrdc.png         (Kiessig-fringe XRR)
%    2. S25_subfilmRC_xrdc.png   (sub + film rocking curves, normalised)
%  The theta-2theta figure (S25_2theta_omega_xrdc.png) was exported earlier.
%
%  Run from anywhere:  run('validation/tushar/runValidation.m')

here    = fileparts(mfilename('fullpath'));
pkgRoot = fullfile(here, '..', '..', 'xrdc-matlab');
inDir   = fullfile(here, 'input');
outDir  = fullfile(here, 'output');
addpath(pkgRoot);

fXRR  = fullfile(inDir, 'TR_S25_SRO_STO(100)_700c_100mT_10500sh_5hz_XRR_05062026_11.txt');
fSub  = fullfile(inDir, 'TR_S25_SRO_STO(100)_700c_100mT_10500sh_5hz_sub RC_05062026.txt');
fFilm = fullfile(inDir, 'TR_S25_SRO_STO(100)_700c_100mT_10500sh_5hz_film RC_05062026.txt');

%% ---------- 1. XRR ----------
fprintf('\n===== XRR =====\n');
scan = xrdc.io.readScan(fXRR);
fprintf('Loaded %s: %d pts, 2theta in [%.3f, %.3f]\n', ...
    scan.identifier, numel(scan.twoTheta), scan.twoTheta(1), scan.twoTheta(end));
res = xrdc.xrr.analyzeFringes(scan);
fprintf('Critical edge ~ %.3f deg. %d fringes in [%.2f, %.2f].\n', ...
    res.twoThetaC, res.nFringesDetected, res.lowerBound, res.upperBound);
if ~isnan(res.thicknessQuadNm)
    fprintf('Thickness d = %.2f +/- %.2f nm (quadratic Kiessig, %d fringes)\n', ...
        res.thicknessQuadNm, res.thicknessQuadSeNm, res.nFringesDetected);
    fprintf('  cross-checks: linear=%.2f nm, FFT=%.2f nm (period %.3f deg, SNR %.1f)\n', ...
        res.thicknessFitNm, res.thicknessFftNm, res.fringePeriodDeg, res.fringeSnr);
    ttl = sprintf('S25 XRR -- d = %.1f nm', res.thicknessQuadNm);
elseif ~isnan(res.thicknessFftNm)
    fprintf('FFT-only thickness d = %.2f nm (period %.3f deg)\n', ...
        res.thicknessFftNm, res.fringePeriodDeg);
    ttl = sprintf('S25 XRR -- d = %.1f nm (FFT)', res.thicknessFftNm);
else
    ttl = 'S25 XRR -- insufficient fringe visibility';
end
scan.peaks = res.peaks;
h = xrdc.plot.plotScan(scan, 'Title', ttl, 'LogY', true, 'ShowPeaks', true);
xlim(h.ax, [0, min(6, scan.twoTheta(end))]);
exportgraphics(h.figure, fullfile(outDir, 'S25_XRR_xrdc.png'), 'Resolution', 600);
fprintf('Saved S25_XRR_xrdc.png\n');

%% ---------- 2. Rocking curves (sub + film, combined, normalised) ----------
fprintf('\n===== Rocking curves =====\n');
[subW, subC, subG, subL] = fitRC(fSub,  'sub');
[flmW, flmC, flmG, flmL] = fitRC(fFilm, 'film');

fig = figure('Color', 'w', 'Position', [100 100 760 520]);
ax  = axes(fig); hold(ax, 'on');
% black = substrate, red = film, normalised to peak = 1, x = detuning
semilogy(ax, subW, subC, 'k-',  'LineWidth', 1.0);
semilogy(ax, flmW, flmC, 'r-',  'LineWidth', 1.0);
set(ax, 'YScale', 'log');
xlim(ax, [-2 2]);
xlabel(ax, '\Delta\omega (\circ)');
ylabel(ax, 'Intensity (a.u.)');
title(ax, 'S25 SRO/STO(100) rocking curves');
legend(ax, { sprintf('Sub  FWHM = %.4f\\circ (gauss)',  subG), ...
             sprintf('Film FWHM = %.4f\\circ (gauss)', flmG) }, ...
       'Location', 'northeast');
grid(ax, 'off'); box(ax, 'on');
exportgraphics(fig, fullfile(outDir, 'S25_subfilmRC_xrdc.png'), 'Resolution', 600);
fprintf('Saved S25_subfilmRC_xrdc.png\n');

fprintf('\n--- FWHM summary (deg) ---\n');
fprintf('sub : gauss %.4f | lorentz %.4f   (Tushar 0.0308)\n', subG, subL);
fprintf('film: gauss %.4f | lorentz %.4f   (Tushar 0.0920)\n', flmG, flmL);

%% ---------- helper ----------
function [dw, cn, fwhmG, fwhmL] = fitRC(fname, label)
    scan = xrdc.io.readScan(fname);
    fprintf('Loaded %s (scanType=%s, %d pts)\n', label, scan.scanType, numel(scan.twoTheta));
    pk = xrdc.peaks.findPeaks(scan, 'MinProminence', max(scan.counts)*0.1);
    assert(~isempty(pk), 'No peak found in %s RC', label);
    [~, iMax] = max([pk.counts]);
    w   = 0.5;
    win = [pk(iMax).twoTheta - w, pk(iMax).twoTheta + w];
    fg  = xrdc.peaks.fitPeak(scan, win, 'Shape', "gauss");
    fl  = xrdc.peaks.fitPeak(scan, win, 'Shape', "lorentz");
    fwhmG = fg.fwhm; fwhmL = fl.fwhm;
    omega0 = fg.twoTheta;
    dw = scan.twoTheta(:) - omega0;
    cn = double(scan.counts(:)) / max(scan.counts);
    cn(cn <= 0) = NaN;   % keep log axis clean
    fprintf('  %s: omega0=%.4f, FWHM gauss=%.4f lorentz=%.4f\n', ...
        label, omega0, fwhmG, fwhmL);
end
