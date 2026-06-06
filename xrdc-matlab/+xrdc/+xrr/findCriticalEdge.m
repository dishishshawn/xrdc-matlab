function info = findCriticalEdge(twoTheta, counts, options)
%FINDCRITICALEDGE  Locate the XRR critical edge robustly.
%   info = xrdc.xrr.findCriticalEdge(twoTheta, counts)
%   info = xrdc.xrr.findCriticalEdge(twoTheta, counts, Name=Value)
%
%   Rigaku SmartLab XRR scans usually contain a direct-beam spike at
%   2θ≈0 that decays through a knife-edge "dip" into the total-external-
%   reflection (TER) plateau, then drops sharply at the critical angle θ_c.
%   The old heuristic — "steepest descent of log10(counts) below 1.5°" —
%   locks onto the knife-edge slope (~0.04°), not the actual XRR edge,
%   which propagates into a 3× thickness error downstream.
%
%   This routine instead:
%     1. Finds the direct-beam dip   = argmin(counts) for 2θ < DirectBeamLimit
%     2. Finds the TER plateau peak  = argmax(counts) in
%                                      [dip + EdgeBuffer, EdgeSearchEnd]
%     3. Treats that plateau peak as 2θ_c (the XRR critical edge).
%
%   Inputs
%     twoTheta : 2θ axis (degrees, column vector)
%     counts   : counts vector, same length
%
%   Name/Value options
%     'DirectBeamLimit' (default 0.50°)
%         Upper bound of the search for the direct-beam dip.
%     'EdgeBuffer'      (default 0.05°)
%         Minimum gap between the direct-beam dip and the TER plateau peak.
%     'EdgeSearchEnd'   (default 1.50°)
%         Upper bound of the TER plateau peak search.
%
%   Output (struct)
%     .twoThetaC        Critical edge in 2θ (degrees)
%     .twoThetaDbDip    2θ at the direct-beam dip (degrees)
%     .countsAtCritical Counts at the critical edge (TER plateau height)
%     .indexCritical    Index into twoTheta of the critical edge
%     .indexDbDip       Index into twoTheta of the direct-beam dip
%
%   See also: xrdc.xrr.analyzeFringes.

    arguments
        twoTheta            (:,1) double
        counts              (:,1) double
        options.DirectBeamLimit (1,1) double {mustBePositive} = 0.50
        options.EdgeBuffer      (1,1) double {mustBeNonnegative} = 0.05
        options.EdgeSearchEnd   (1,1) double {mustBePositive} = 1.50
    end

    if numel(twoTheta) ~= numel(counts)
        error('xrdc:xrr:sizeMismatch', ...
            'twoTheta and counts must have the same length.');
    end
    if numel(twoTheta) < 5
        error('xrdc:xrr:tooFewPoints', ...
            'Need at least 5 points to detect a critical edge.');
    end

    % 1. Direct-beam dip
    preDB = twoTheta < options.DirectBeamLimit;
    if ~any(preDB)
        error('xrdc:xrr:noDirectBeamWindow', ...
            'No points below DirectBeamLimit=%.3f° for direct-beam dip search.', ...
            options.DirectBeamLimit);
    end
    [~, relDip] = min(counts(preDB));
    iDbDip = relDip;

    % 2. TER plateau peak after the dip
    searchMask = twoTheta > twoTheta(iDbDip) + options.EdgeBuffer & ...
                 twoTheta < options.EdgeSearchEnd;
    if ~any(searchMask)
        error('xrdc:xrr:noEdgeWindow', ...
            ['No points in [%.3f°, %.3f°] for the TER-plateau search. ' ...
             'Increase EdgeSearchEnd or check the scan range.'], ...
            twoTheta(iDbDip) + options.EdgeBuffer, options.EdgeSearchEnd);
    end
    ws = find(searchMask, 1, 'first');
    we = find(searchMask, 1, 'last');
    [~, relTer] = max(counts(ws:we));
    iCrit = ws + relTer - 1;

    info = struct( ...
        'twoThetaC',        twoTheta(iCrit), ...
        'twoThetaDbDip',    twoTheta(iDbDip), ...
        'countsAtCritical', counts(iCrit), ...
        'indexCritical',    iCrit, ...
        'indexDbDip',       iDbDip);
end
