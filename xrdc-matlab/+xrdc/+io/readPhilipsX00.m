function scan = readPhilipsX00(path)
%READPHILIPSX00  Read a Philips HR-XRDSCAN (.x00) text scan.
%   scan = xrdc.io.readPhilipsX00(path)
%
%   Matches ReadFileTypePhilips in xrdc1.pas:991. The format is a
%   line-based key/value header (FIRSTANGLE, STEPWIDTH, NROFDATA,
%   LABDAALPHA1, etc.) followed by a SCANDATA marker and integer counts
%   one per line.
%
%   Input
%     path : path to .x00 file
%
%   Output
%     scan : struct (xrdc.io.emptyScan shape)

    arguments
        path (1,1) string
    end

    scan = xrdc.io.emptyScan();
    scan.sourcePath   = path;
    scan.sourceFormat = "philipsX00";
    [~, stem, ~]      = fileparts(path);
    scan.identifier   = string(stem);

    raw   = string(fileread(path));
    lines = splitlines(raw);

    md = struct();
    dataStart = -1;
    for i = 1:numel(lines)
        line = strtrim(lines(i));
        if startsWith(upper(line), "SCANDATA")
            dataStart = i + 1;
            break
        end
        if contains(line, "=")
            parts = split(line, "=");
            if numel(parts) >= 2
                key = matlab.lang.makeValidName(strtrim(parts{1}));
                val = strtrim(strjoin(parts(2:end), "="));
                val = replace(val, ",", ".");
                valNum = str2double(val);
                if isfinite(valNum)
                    md.(key) = valNum;
                else
                    md.(key) = string(val);
                end
            end
        end
    end

    if dataStart < 0
        error('xrdc:io:noData', 'No SCANDATA block found in %s.', path);
    end

    % parse counts (one integer per remaining non-empty line)
    counts = [];
    for i = dataStart:numel(lines)
        line = strtrim(lines(i));
        if line == "" || startsWith(line, "!")
            continue
        end
        line = replace(line, ",", ".");
        tok = regexp(line, '[-+]?\d*\.?\d+', 'match');
        if ~isempty(tok)
            counts(end+1, 1) = str2double(tok{1});  %#ok<AGROW>
        end
    end

    if isempty(counts)
        error('xrdc:io:emptyScan', 'No counts parsed from %s.', path);
    end

    firstAngle = getNum(md, 'FIRSTANGLE', NaN);
    stepWidth  = getNum(md, 'STEPWIDTH',  NaN);
    nPts       = getNum(md, 'NROFDATA',   numel(counts));

    if isnan(firstAngle) || isnan(stepWidth)
        error('xrdc:io:missingHeader', ...
            'FIRSTANGLE / STEPWIDTH missing from %s.', path);
    end

    if nPts ~= numel(counts)
        warning('xrdc:io:countMismatch', ...
            'Header NROFDATA=%d but %d counts parsed from %s.', ...
            nPts, numel(counts), path);
    end

    twoTheta = firstAngle + (0:numel(counts)-1).' * stepWidth;
    scan.twoTheta = twoTheta;
    scan.counts   = counts;
    scan.scanType = "twoThetaOmega";
    scan.metadata = md;
    lambdaKey = {'LABDAALPHA1', 'LAMBDAALPHA1', 'WAVELENGTH'};
    for kk = 1:numel(lambdaKey)
        if isfield(md, lambdaKey{kk})
            v = md.(lambdaKey{kk});
            if isnumeric(v) && isfinite(v)
                scan.lambda = v;
                break
            end
        end
    end
end

function v = getNum(s, field, fallback)
    if isfield(s, field) && isnumeric(s.(field))
        v = s.(field);
    else
        v = fallback;
    end
end
