function scans = readXrdmlArea(path)
%READXRDMLAREA  Read a multi-scan PANalytical XRDML reciprocal-space map.
%   scans = xrdc.io.readXrdmlArea(path)
%
%   PANalytical area scans embed many <scan> elements inside one .xrdml
%   file — one 2θ slice per ω step. The single-scan reader
%   (xrdc.io.readXrdml) returns only the first slice, which is correct
%   for standalone scans but collapses an RSM to a line. This reader
%   returns a struct array with one element per slice, suitable for
%   xrdc.rsm.toReciprocalSpace and xrdc.plot.plotRsm.
%
%   Input
%     path : path to .xrdml file (must be a multi-scan area measurement).
%
%   Output
%     scans : (1×N) struct array in xrdc.io.emptyScan() shape.
%             Each element has scanType = "area" and secondAxis set to
%             the ω value for that slice.

    arguments
        path (1,1) string
    end

    doc  = xmlread(char(path));
    root = doc.getDocumentElement();

    % Top-level <usedWavelength> (shared across all slices).
    lambda = NaN;
    wlNode = findFirst(root, 'usedWavelength');
    if ~isempty(wlNode)
        kA1 = findFirst(wlNode, 'kAlpha1');
        if ~isempty(kA1)
            lambda = str2double(char(kA1.getTextContent()));
        end
    end

    % Enumerate all <scan> elements — each is one 2θ slice at a fixed ω.
    scanNodes = root.getElementsByTagName('scan');
    nSlices   = scanNodes.getLength();
    if nSlices == 0
        error('xrdc:io:noScanNode', 'No <scan> elements in %s.', path);
    end

    [~, stem, ~] = fileparts(path);
    identifier   = string(stem);

    scans = repmat(xrdc.io.emptyScan(), 1, nSlices);

    for i = 0:(nSlices - 1)
        scanNode = scanNodes.item(i);
        dataPointsNode = findFirst(scanNode, 'dataPoints');
        if isempty(dataPointsNode)
            error('xrdc:io:noDataPoints', ...
                'Slice %d in %s has no <dataPoints>.', i + 1, path);
        end

        % Find 2θ start/end and the ω commonPosition for this slice.
        [tt0, tt1, omegaVal] = readSlicePositions(dataPointsNode);

        intensityNode = findFirst(dataPointsNode, 'intensities');
        if isempty(intensityNode)
            intensityNode = findFirst(dataPointsNode, 'counts');
        end
        if isempty(intensityNode)
            error('xrdc:io:noIntensities', ...
                'Slice %d in %s has no <intensities>.', i + 1, path);
        end
        txt    = strtrim(string(char(intensityNode.getTextContent())));
        counts = sscanf(txt, '%f');
        n      = numel(counts);
        if n < 2
            error('xrdc:io:tooShort', ...
                'Slice %d in %s has only %d points.', i + 1, path, n);
        end
        step    = (tt1 - tt0) / (n - 1);
        twoThet = tt0 + (0:n-1).' * step;

        s = xrdc.io.emptyScan();
        s.sourcePath     = path;
        s.sourceFormat   = "xrdml";
        s.identifier     = sprintf("%s#%d", identifier, i + 1);
        s.scanType       = "area";
        s.twoTheta       = twoThet;
        s.counts         = counts;
        s.lambda         = lambda;
        s.secondAxis     = omegaVal;
        s.secondAxisName = "Omega";
        s.metadata = struct( ...
            'scanAxis',      '2Theta', ...
            'startPosition', tt0, ...
            'endPosition',   tt1, ...
            'nPoints',       n, ...
            'sliceIndex',    i + 1, ...
            'nSlices',       nSlices);
        scans(i + 1) = s;
    end

    % Sort by omega ascending so plotRsm sees monotonic slices.
    [~, order] = sort([scans.secondAxis]);
    scans = scans(order);
end

% --- helpers -----------------------------------------------------------

function [tt0, tt1, omegaVal] = readSlicePositions(dataPointsNode)
%READSLICEPOSITIONS  Pull 2θ start/end and ω common position from one slice.
    posNodes = dataPointsNode.getElementsByTagName('positions');
    tt0 = NaN; tt1 = NaN; omegaVal = NaN;
    for k = 0:(posNodes.getLength() - 1)
        node = posNodes.item(k);
        if node.getParentNode() ~= dataPointsNode
            continue  % skip nested <positions> from deeper children
        end
        axisName = "";
        if node.hasAttribute('axis')
            axisName = string(char(node.getAttribute('axis')));
        end
        startNode  = findFirst(node, 'startPosition');
        endNode    = findFirst(node, 'endPosition');
        commonNode = findFirst(node, 'commonPosition');

        if strcmpi(axisName, "2Theta") && ~isempty(startNode) && ~isempty(endNode)
            tt0 = str2double(char(startNode.getTextContent()));
            tt1 = str2double(char(endNode.getTextContent()));
        elseif strcmpi(axisName, "Omega") && ~isempty(commonNode)
            omegaVal = str2double(char(commonNode.getTextContent()));
        end
    end
    if isnan(tt0) || isnan(tt1)
        error('xrdc:io:no2Theta', 'Slice missing 2θ start/end positions.');
    end
    if isnan(omegaVal)
        error('xrdc:io:noOmega', 'Slice missing Omega commonPosition.');
    end
end

function node = findFirst(parent, tagName)
    nodes = parent.getElementsByTagName(tagName);
    if nodes.getLength() == 0
        node = [];
    else
        node = nodes.item(0);
    end
end
