%% Demo: φ scan — linear Y, identify 4-fold substrate symmetry peaks.
%  A well-aligned cubic/tetragonal (100) substrate shows 4 peaks per 360°
%  in the (101) family, spaced ~90° apart.

addpath(fileparts(fileparts(mfilename('fullpath'))));

% Phi scan PANalytical copy lives in rexdrctomatlabport; Rigaku archive
% doesn't have an equivalent phi export.
dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '..', '..', 'rexdrctomatlabport');
defaultFname = 'HP_TiO2 101 phi scan 3pixcel 1_8 slit.xrdml';
if ~exist('fname', 'var') || isempty(fname) || ~isfile(fullfile(dataDir, fname))
    fname = defaultFname;
end
scan    = xrdc.io.readScan(fullfile(dataDir, fname));
fprintf('Loaded %s  (scanType = %s)\n', scan.identifier, scan.scanType);

%% Find 4-fold peaks
pk = xrdc.peaks.findPeaks(scan, ...
    'MinProminence',  max(scan.counts) * 0.1, ...
    'MinSeparation',  30);       % ≥30° apart for a 4-fold symmetry
fprintf('Found %d peaks.  Expected 4 for a (100)-cut cubic substrate.\n', numel(pk));
if numel(pk) == 4
    spacings = diff([pk.twoTheta]);
    fprintf('Peak-to-peak spacings: %s\n', ...
        strjoin(arrayfun(@(x) sprintf('%.2f°', x), spacings, 'UniformOutput', false), ', '));
end

%% Plot (linear Y — φ scans don't benefit from log)
scan.peaks = pk;
h = xrdc.plot.plotScan(scan, ...
    'Title',     "φ scan — (101) family, 4-fold symmetry check", ...
    'LogY',      false, ...
    'ShowPeaks', true);
xlabel(h.ax, '\phi (\circ)');

[~, stem, ~] = fileparts(fname);
outPath = fullfile(pwd, sprintf('phi_%s.png', stem));
exportgraphics(h.figure, outPath, 'Resolution', 600);
fprintf('Saved: %s\n', outPath);
