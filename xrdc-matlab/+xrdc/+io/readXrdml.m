function scan = readXrdml(path)
%READXRDML  Read a PANalytical XRDML (.xrdml) scan.
%   scan = xrdc.io.readXrdml(path)
%
%   Matches ReadFileTypeXML in xrdc1.pas:494. Extracts:
%     - scanAxis attribute of <scan> → scanType
%     - startPosition / endPosition → 2θ grid
%     - intensities → counts (whitespace-separated)
%     - usedWavelength → lambda (Å)
%
%   For θ-2θ scans the companion axis (usually omega) is captured as
%   secondAxis / secondAxisName.
%
%   Input
%     path : path to .xrdml file
%
%   Output
%     scan : struct (xrdc.io.emptyScan shape)

    arguments
        path (1,1) string
    end

    scan = xrdc.io.emptyScan();
    scan.sourcePath   = path;
    scan.sourceFormat = "xrdml";
    [~, stem, ~]      = fileparts(path);
    scan.identifier   = string(stem);

    doc = xmlread(char(path));
    root = doc.getDocumentElement();

    % handle both <xrdMeasurement> and <xrdMeasurements>/<xrdMeasurement>
    scanNode = findFirst(root, 'scan');
    if isempty(scanNode)
        error('xrdc:io:noScanNode', 'No <scan> node in %s.', path);
    end

    % scanAxis → scanType
    scanAxis = "";
    if scanNode.hasAttribute('scanAxis')
        scanAxis = string(char(scanNode.getAttribute('scanAxis')));
    end
    scan.scanType = normaliseScanType(scanAxis);

    dataPointsNode = findFirst(scanNode, 'dataPoints');
    if isempty(dataPointsNode)
        error('xrdc:io:noDataPoints', 'No <dataPoints> in %s.', path);
    end

    % wavelength (<usedWavelength>/<kAlpha1> etc.)
    wlNode = findFirst(root, 'usedWavelength');
    if ~isempty(wlNode)
        kA1 = findFirst(wlNode, 'kAlpha1');
        if ~isempty(kA1)
            scan.lambda = str2double(char(kA1.getTextContent()));
        end
    end

    % positions: one or more <positions axis="..."> elements inside
    % <dataPoints>. The "primary" axis matches scanAxis; the others are
    % constants (secondary).
    posNodes = findAll(dataPointsNode, 'positions');
    primaryAxis = primaryAxisName(scan.scanType);
    primaryStart = NaN; primaryEnd = NaN;
    secondaryAxis = ""; secondaryValue = NaN;

    for i = 1:numel(posNodes)
        posNode = posNodes{i};
        axisName = "";
        if posNode.hasAttribute('axis')
            axisName = string(char(posNode.getAttribute('axis')));
        end
        startNode = findFirst(posNode, 'startPosition');
        endNode   = findFirst(posNode, 'endPosition');
        commonNode = findFirst(posNode, 'commonPosition');

        if ~isempty(startNode) && ~isempty(endNode)
            startVal = str2double(char(startNode.getTextContent()));
            endVal   = str2double(char(endNode.getTextContent()));
            if strcmpi(axisName, primaryAxis) || isnan(primaryStart)
                primaryStart = startVal;
                primaryEnd   = endVal;
            else
                secondaryAxis  = axisName;
                secondaryValue = (startVal + endVal) / 2;
            end
        elseif ~isempty(commonNode)
            secondaryAxis  = axisName;
            secondaryValue = str2double(char(commonNode.getTextContent()));
        end
    end

    % intensities
    intensityNode = findFirst(dataPointsNode, 'intensities');
    if isempty(intensityNode)
        intensityNode = findFirst(dataPointsNode, 'counts');
    end
    if isempty(intensityNode)
        error('xrdc:io:noIntensities', 'No <intensities>/<counts> in %s.', path);
    end
    txt = strtrim(string(char(intensityNode.getTextContent())));
    counts = sscanf(txt, '%f')';
    counts = counts(:);

    n = numel(counts);
    if n < 2
        error('xrdc:io:tooShort', 'Too few intensity values in %s (%d).', path, n);
    end

    step = (primaryEnd - primaryStart) / (n - 1);
    twoTheta = primaryStart + (0:n-1).' * step;

    scan.twoTheta       = twoTheta;
    scan.counts         = counts;
    scan.secondAxis     = secondaryValue;
    scan.secondAxisName = string(secondaryAxis);
    scan.metadata = struct( ...
        'scanAxis',    scanAxis, ...
        'startPosition', primaryStart, ...
        'endPosition',   primaryEnd, ...
        'nPoints',       n);
end

% ---- XML helpers ----

function node = findFirst(parent, tagName)
    nodes = parent.getElementsByTagName(tagName);
    if nodes.getLength() == 0
        node = [];
    else
        node = nodes.item(0);
    end
end

function out = findAll(parent, tagName)
    nodes = parent.getElementsByTagName(tagName);
    out = cell(nodes.getLength(), 1);
    for i = 1:nodes.getLength()
        out{i} = nodes.item(i - 1);
    end
end

function st = normaliseScanType(scanAxisStr)
    s = lower(string(scanAxisStr));
    if contains(s, "2theta") && contains(s, "omega")
        st = "twoThetaOmega";
    elseif s == "omega"
        st = "omega";
    elseif s == "phi"
        st = "phi";
    elseif s == "psi"
        st = "psi";
    elseif contains(s, "reciprocal")
        st = "area";
    else
        st = "unknown";
    end
end

function ax = primaryAxisName(scanType)
    switch scanType
        case "twoThetaOmega", ax = "2Theta-Omega";
        case "omega",         ax = "Omega";
        case "phi",           ax = "Phi";
        case "psi",           ax = "Psi";
        otherwise,            ax = "";
    end
end
