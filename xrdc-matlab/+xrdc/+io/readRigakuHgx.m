function scan = readRigakuHgx(path)
%READRIGAKUHGX  Read a Rigaku SmartLab Studio II GlobalFit project (.hgx).
%   scan = xrdc.io.readRigakuHgx(path)
%
%   .hgx files are HDF5 containers (magic bytes "\x89HDF\r\n\x1a\n") written
%   by Rigaku's SmartLab Studio II GlobalFit / reflectivity analysis. The
%   measured curve lives under a fixed path:
%
%     /current/data/datasets/0/x_raw   — angle (°), same quantity the Rigaku
%                                         ASCII export puts in column 1
%     /current/data/datasets/0/y_raw   — intensity (cps)
%     /current/data/datasets/0/y_sim   — GlobalFit simulated curve (optional)
%
%   x_raw/y_raw were verified bit-identical to the matching .txt export, so
%   this reader yields the same scan as readRigakuTxt would on the .txt twin.
%   The .hgx additionally carries the layer-model fit (y_sim) and parameters,
%   which we surface in metadata for downstream use.
%
%   Scan type is inferred from the filename (Rigaku stores reflectivity,
%   θ-2θ, RC, etc. in the same container shape) — see readRigakuTxt.
%
%   Input
%     path : path to a .hgx file
%
%   Output
%     scan : struct in xrdc.io.emptyScan() shape

    arguments
        path (1,1) string
    end

    if ~isfile(path)
        error('xrdc:io:notFound', 'File not found: %s', path);
    end

    base = '/current/data/datasets/0';
    try
        x = h5read(char(path), [base '/x_raw']);
        y = h5read(char(path), [base '/y_raw']);
    catch ME
        error('xrdc:io:hgxStructure', ...
            'Unexpected .hgx layout (could not read %s/x_raw,y_raw): %s', ...
            base, ME.message);
    end

    x = double(x(:));
    y = double(y(:));
    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);
    if numel(x) < 2
        error('xrdc:io:noData', 'No usable curve in %s.', path);
    end

    scan = xrdc.io.emptyScan();
    scan.sourcePath   = path;
    scan.sourceFormat = "rigakuHgx";
    [~, stem, ~]      = fileparts(path);
    scan.identifier   = string(stem);
    scan.twoTheta     = x;
    scan.counts       = y;
    scan.scanType     = inferScanType(stem);
    scan.lambda       = 1.5406;   % Cu Kα1, Rigaku SmartLab default

    % Optionally surface the GlobalFit simulated curve for cross-checks.
    ySim = [];
    try
        ySim = double(h5read(char(path), [base '/y_sim']));
        ySim = ySim(:);
        ySim = ySim(keep);
    catch
        ySim = [];
    end

    step = x(2) - x(1);
    scan.metadata = struct( ...
        'nPoints',     numel(x), ...
        'step',        step, ...
        'hasSimCurve', ~isempty(ySim) && numel(ySim) == numel(x));
    if scan.metadata.hasSimCurve
        scan.metadata.ySim = ySim;
    end
end

% =====================================================================
function st = inferScanType(stem)
    % Mirror readRigakuTxt: the container shape is identical for all scan
    % types, so the filename is the only cue.
    s = lower(string(stem));
    if contains(s, " rc") || contains(s, "_rc") || contains(s, " rocking") ...
            || endsWith(s, "rc")
        st = "omega";
    elseif contains(s, "phi")
        st = "phi";
    elseif contains(s, "psi")
        st = "psi";
    else
        st = "twoThetaOmega";   % XRR, 2theta omega, th2th, etc.
    end
end
