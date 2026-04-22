function scan = readRigakuTxt(path)
%READRIGAKUTXT  Read Rigaku ASCII export (.txt, tab-delimited).
%   scan = xrdc.io.readRigakuTxt(path)
%
%   Parses the custom Rigaku SmartLab ASCII export. Two layout variants
%   are accepted:
%
%     Headered (most files):
%       <BOM><SampleName>\n
%       2θ, °\tIntensity, cps\n
%       <x>\t<y>\n
%       ...
%
%     Headerless:
%       <x>\t<y>\n
%       ...
%
%   Detection: parses line 1 as two tab-separated floats. Success → headerless.
%   Failure → treat first 2 lines as header, start data at line 3.
%
%   The axis label in the header reads "2θ, °" for every Rigaku export,
%   even for rocking curves where the x-axis is ω. Scan type is therefore
%   inferred from the *filename* (case-insensitive):
%     - contains "RC"               → "omega"   (rocking curve)
%     - contains "phi"              → "phi"
%     - contains "XRR"              → "twoThetaOmega"  (specular reflectivity)
%     - contains "2theta" or "2 theta" → "twoThetaOmega"
%     - otherwise                   → "twoThetaOmega"
%
%   Input
%     path : path to file
%
%   Output
%     scan : struct in xrdc.io.emptyScan() shape

    arguments
        path (1,1) string
    end

    if ~isfile(path)
        error('xrdc:io:notFound', 'File not found: %s', path);
    end

    scan = xrdc.io.emptyScan();
    scan.sourcePath   = path;
    scan.sourceFormat = "rigakuTxt";
    [~, stem, ~]      = fileparts(path);
    scan.identifier   = string(stem);

    % Read all lines, strip UTF-8 BOM on first line if present
    raw = fileread(char(path));
    if ~isempty(raw) && double(raw(1)) == 65279   % UTF-8 BOM as char
        raw(1) = [];
    end
    % Also handle the 3-byte BOM sequence if read as bytes
    if numel(raw) >= 3 && all(double(raw(1:3)) == [239 187 191])
        raw(1:3) = [];
    end
    lines = splitlines(string(raw));
    lines(strlength(strtrim(lines)) == 0) = [];   % drop blank lines

    if numel(lines) < 2
        error('xrdc:io:tooShort', 'File has fewer than 2 lines: %s', path);
    end

    % Decide headered vs headerless by trying to parse line 1 as [x y].
    [headered, dataStart, sampleName, colLabels] = detectHeader(lines);

    % Parse remaining lines as tab (or whitespace) separated floats
    dataLines = lines(dataStart:end);
    n = numel(dataLines);
    twoTheta = zeros(n, 1);
    counts   = zeros(n, 1);
    keep     = true(n, 1);
    for i = 1:n
        vals = sscanf(char(dataLines(i)), '%f');
        if numel(vals) < 2
            keep(i) = false;
            continue;
        end
        twoTheta(i) = vals(1);
        counts(i)   = vals(2);
    end
    twoTheta = twoTheta(keep);
    counts   = counts(keep);

    if numel(twoTheta) < 2
        error('xrdc:io:noData', 'Could not parse any data rows in %s.', path);
    end

    % Uniform-step check (warn only; Delphi does the same — xrdc1.pas:1235)
    step  = twoTheta(2) - twoTheta(1);
    diffs = diff(twoTheta);
    if any(abs(diffs - step) > 1e-3)
        warning('xrdc:io:nonUniformStep', ...
            '2θ step is non-uniform in %s (max dev %.4g°).', path, max(abs(diffs - step)));
    end

    scan.twoTheta = twoTheta;
    scan.counts   = counts;
    scan.scanType = inferScanType(stem);
    scan.lambda   = 1.5406;   % Cu Kα1, Rigaku SmartLab default (override if needed)

    scan.metadata = struct( ...
        'headered',   headered, ...
        'sampleName', sampleName, ...
        'colLabels',  colLabels, ...
        'nPoints',    numel(twoTheta), ...
        'step',       step);
end

% =====================================================================
function [headered, dataStart, sampleName, colLabels] = detectHeader(lines)
    % Try to parse line 1 as two floats
    vals = sscanf(char(lines(1)), '%f');
    if numel(vals) >= 2
        headered   = false;
        dataStart  = 1;
        sampleName = "";
        colLabels  = "";
        return;
    end

    % Headered: line 1 = sample name, line 2 = column labels, line 3+ = data
    headered   = true;
    sampleName = strtrim(lines(1));
    if numel(lines) >= 2
        colLabels = strtrim(lines(2));
        dataStart = 3;
    else
        colLabels = "";
        dataStart = 2;
    end
    % Edge case: some files have only a 1-line header.  If line 2 parses as data,
    % treat line 2 as the start of data.
    if numel(lines) >= 2
        vals2 = sscanf(char(lines(2)), '%f');
        if numel(vals2) >= 2
            colLabels = "";
            dataStart = 2;
        end
    end
end

function st = inferScanType(stem)
    s = lower(string(stem));
    if contains(s, " rc") || contains(s, "_rc") || contains(s, " rocking") ...
            || endsWith(s, "rc")
        st = "omega";
    elseif contains(s, "phi")
        st = "phi";
    elseif contains(s, "psi")
        st = "psi";
    else
        % Everything else from Rigaku SmartLab exports is a coupled θ-2θ:
        %   XRR, 2theta omega, 2 theta omega, th2th, etc.
        st = "twoThetaOmega";
    end
end
