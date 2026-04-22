%DEMORIG AKUWORKFLOW  Complete workflow: load, analyze, plot Rigaku data
%   Shows the typical steps when using xrdc-matlab with files from
%   the Rigaku SmartLab machine.
%
%   This demo loads a rocking curve (RC), detects the peak, overlays
%   it on the plot, and computes derivatives for analysis.

clear; close all; addpath(fullfile(fileparts(mfilename), '..'));

%% Step 1: Locate your Rigaku files
% When you export from the Rigaku SmartLab machine, you get a .txt file.
% The filename tells you what kind of scan it is (look for keywords):
%   "RC"  or "rocking" -> omega/theta rocking curve  (scanType = "omega")
%   "phi" -> phi scan                                 (scanType = "phi")
%   "2theta omega" -> reciprocal space map slice     (scanType = "twoThetaOmega")
%   "XRR" -> X-ray reflectivity                      (scanType = "twoThetaOmega")

rigaku_dir = fullfile(fileparts(mfilename), '..', '..', '..', 'rexdrctomatlabport_rigakudatasets');
rc_file = fullfile(rigaku_dir, 'TR_S05_PTO_STO(100)_600c_200mT_1000sh_2hz_film RC_04092026.txt');

if ~isfile(rc_file)
    error('Sample Rigaku file not found. Check the rigaku_dir path.');
end

%% Step 2: Load the scan
% xrdc.io.readScan automatically detects the format and calls the right parser.
scan = xrdc.io.readScan(rc_file);

fprintf('Loaded scan: %s\n', scan.identifier);
fprintf('  Type: %s\n', scan.scanType);
fprintf('  Points: %d\n', numel(scan.twoTheta));
fprintf('  Wavelength: %.4f Å (Cu Kα1)\n', scan.lambda);
fprintf('  Count range: [%d, %d]\n', min(scan.counts), max(scan.counts));

%% Step 3: Plot the raw scan
figure('Name', 'Raw Rocking Curve'); clf;
semilogy(scan.twoTheta, max(scan.counts, 1), 'b-', 'LineWidth', 1.5);
xlabel('2\theta (degrees)');
ylabel('Intensity (counts)');
title('Rocking Curve (Film Peak)');
grid on;

%% Step 4: Compute derivatives for peak detection
% The first derivative tells us where the peak slope is steepest.
% The second derivative is zero at the peak (inflection point).
[slope, slope2] = xrdc.signal.derivatives(scan.twoTheta, scan.counts, 11, 3);

%% Step 5: Find peaks using the second derivative
% Peaks occur where slope2 crosses zero (going negative).
% Simple peak detection: find local maxima without Signal Processing Toolbox
peaks = find(diff(sign(diff(scan.counts))) < 0) + 1;  % Local maxima
peaks(scan.counts(peaks) < max(scan.counts) * 0.1) = [];  % Filter noise

if ~isempty(peaks)
    fprintf('\nDetected %d peak(s):\n', numel(peaks));
    for p = peaks.'
        fprintf('  Index %d: θ = %.4f°, counts = %d\n', ...
            p, scan.twoTheta(p), scan.counts(p));
    end
end

%% Step 6: Overlay peak locations on plot
hold on;
semilogy(scan.twoTheta(peaks), scan.counts(peaks), 'r*', 'MarkerSize', 12, 'DisplayName', 'Peaks');
legend('Location', 'best');
hold off;

%% Step 7: Plot derivatives to visualize peak detection
fig2 = figure('Name', 'Derivatives');
ax1 = subplot(2, 1, 1);
plot(ax1, scan.twoTheta, slope, 'b', 'LineWidth', 1.5);
xlabel(ax1, '2\theta (degrees)');
ylabel(ax1, 'dI/d\theta (counts/°)');
title(ax1, 'First Derivative (Slope)');
grid(ax1, 'on');

ax2 = subplot(2, 1, 2);
plot(ax2, scan.twoTheta, slope2, 'r', 'LineWidth', 1.5);
hold(ax2, 'on');
yline(ax2, 0, 'k--', 'LineWidth', 1);
plot(ax2, scan.twoTheta(peaks), slope2(peaks), 'go', 'MarkerSize', 8, 'DisplayName', 'Zero Crossings');
xlabel(ax2, '2\theta (degrees)');
ylabel(ax2, 'd²I/d\theta² (counts/°²)');
title(ax2, 'Second Derivative (Curvature)');
legend(ax2, 'Location', 'best');
grid(ax2, 'on');
hold(ax2, 'off');

%% Step 8 (Optional): Load multiple scans from a folder
% If you have several related scans (e.g., different substrate and film peaks),
% you can load them all at once.
%
% Example: load all "RC" scans from a folder
% rc_files = dir(fullfile(rigaku_dir, '*RC*'));
% scans = [];
% for f = {rc_files.name}
%     try
%         scans = [scans, xrdc.io.readScan(fullfile(rigaku_dir, f{1}))];
%     catch
%         % Skip files that fail to parse
%     end
% end

fprintf('\n✓ Demo complete. Check figures for results.\n');
fprintf('  Next: Export peak positions, fit peak widths, or compare multiple scans.\n');
