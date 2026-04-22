function scan = readTextScan(path)
%READTEXTSCAN  Read a two-column (2θ, counts) plain-text scan.
%   scan = xrdc.io.readTextScan(path)
%
%   Matches ReadFileTypeText in xrdc1.pas:1235. Behaviour:
%     - Accepts period OR comma as decimal separator (legacy German locale).
%     - Skips header/comment lines that don't parse as two numbers.
%     - If any count exceeds 2e9, divides the entire column by 10
%       (undoes a fixed-point × 10 encoding some vendor exports use).
%     - Warns if the step width is not uniform to 0.001°.
%
%   Input
%     path : path to file (char or string)
%
%   Output
%     scan : struct with the shape defined by xrdc.io.emptyScan()

    arguments
        path (1,1) string
    end

    scan = xrdc.io.emptyScan();
    scan.sourcePath   = path;
    scan.sourceFormat = "text";
    [~, stem, ~]      = fileparts(path);
    scan.identifier   = string(stem);

    lines = splitlines(string(fileread(path)));
    twoTheta = [];
    counts   = [];
    for i = 1:numel(lines)
        line = strtrim(lines(i));
        if line == "" || startsWith(line, "#") || startsWith(line, "%")
            continue
        end
        line = replace(line, ",", ".");   % German-locale tolerance
        tok = regexp(line, '[-+]?\d*\.?\d+([eE][-+]?\d+)?', 'match');
        if numel(tok) < 2
            continue    % skip non-numeric headers
        end
        twoTheta(end+1, 1) = str2double(tok{1});  %#ok<AGROW>
        counts(end+1, 1)   = str2double(tok{2});  %#ok<AGROW>
    end

    if isempty(twoTheta)
        error('xrdc:io:emptyScan', 'No numeric data rows found in %s.', path);
    end

    if any(counts > 2e9)
        counts = counts / 10;
    end

    step = diff(twoTheta);
    if ~isempty(step) && (max(step) - min(step)) > 0.001
        warning('xrdc:io:nonUniformStep', ...
            'Non-uniform 2θ step in %s (range %.4f°–%.4f°).', path, min(step), max(step));
    end

    scan.twoTheta = twoTheta;
    scan.counts   = counts;
    scan.scanType = "twoThetaOmega";   % best guess for plain two-column data
end
