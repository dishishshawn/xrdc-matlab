function [keep, ghosts] = filterGhostPeaks(twoTheta, counts, lambda, opts)
%FILTERGHOSTPEAKS  Flag Cu-Kbeta / W-Lalpha ghost peaks of strong reflections.
%   [keep, ghosts] = xrdc.peaks.filterGhostPeaks(twoTheta, counts, lambda)
%
%   A reflection with spacing d (from the main wavelength) also diffracts
%   the Cu Kbeta (1.3922 A) and W Lalpha1 (1.4763 A) tube lines, leaving
%   small satellite peaks at 2*asin(lambdaGhost/2d), below the parent.
%   Any peak within PositionTol of a predicted ghost position whose
%   intensity is at most MaxRatio x the parent is flagged.
%
%   Caveat: a genuine weak peak that happens to fall within PositionTol of
%   a predicted ghost position will also be flagged as a ghost.  Callers
%   should inspect the `ghosts` output table to audit all removals before
%   accepting the result.
%
%   NaN-count peaks are never flagged — unmeasured peaks cannot be ranked
%   as strong or weak, so the ghost test is skipped for them.  Only the
%   all-NaN case triggers the xrdc:peaks:noIntensity warning.
%
%   Inputs
%     twoTheta : peak positions, degrees (vector)
%     counts   : peak intensities (same size; all-NaN -> warning
%                xrdc:peaks:noIntensity and nothing is flagged, because
%                strong/weak cannot be ranked)
%     lambda   : main wavelength in A
%
%   Name-Value
%     GhostLambdas (default [1.3922 1.4763])  — Cu Kbeta1, W Lalpha1 (A)
%     PositionTol  (default 0.15)             — degrees 2-theta
%     MaxRatio     (default 0.3)              — ghost/parent intensity cap
%
%   Outputs
%     keep   : logical mask, true = real peak
%     ghosts : table (twoTheta, counts, parentTwoTheta, ghostLambda)
%
%   See also xrdc.lattice.identifyMaterial.

    arguments
        twoTheta (:,1) double
        counts   (:,1) double
        lambda   (1,1) double {mustBePositive}
        opts.GhostLambdas (1,:) double = [1.3922, 1.4763]
        opts.PositionTol  (1,1) double {mustBePositive} = 0.15
        opts.MaxRatio     (1,1) double {mustBePositive} = 0.3
    end

    n = numel(twoTheta);
    assert(numel(counts) == n, 'xrdc:peaks:sizeMismatch', ...
        'twoTheta and counts must be the same length.');
    keep = true(n, 1);
    ghosts = table('Size', [0 4], ...
        'VariableTypes', {'double','double','double','double'}, ...
        'VariableNames', {'twoTheta','counts','parentTwoTheta','ghostLambda'});

    if all(isnan(counts))
        warning('xrdc:peaks:noIntensity', ...
            ['No peak intensities provided - ghost filtering skipped ', ...
             '(cannot rank strong vs weak peaks).']);
        return
    end

    % Strongest first: a ghost can never claim its own parent.
    [~, order] = sort(counts, 'descend', 'MissingPlacement', 'last');
    for p = order(:).'
        if ~keep(p) || isnan(counts(p)), continue, end
        d = xrdc.lattice.twoThetaToD(twoTheta(p), lambda);
        for gl = opts.GhostLambdas
            g2t = xrdc.lattice.dToTwoTheta(d, gl);
            if isnan(g2t), continue, end
            hits = find(keep & (1:n).' ~= p ...
                & abs(twoTheta - g2t) <= opts.PositionTol ...
                & counts <= opts.MaxRatio * counts(p));
            for j = hits(:).'
                keep(j) = false;
                ghosts(end+1, :) = {twoTheta(j), counts(j), twoTheta(p), gl}; %#ok<AGROW>
            end
        end
    end
end
