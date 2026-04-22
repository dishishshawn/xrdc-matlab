function T = simulatePattern(lattice, hklRange, lambda, opts)
%SIMULATEPATTERN  Predict Bragg peak positions from lattice parameters.
%   T = xrdc.lattice.simulatePattern(lattice, hklRange, lambda)
%   T = xrdc.lattice.simulatePattern(lattice, hklRange, lambda, Name=Value)
%
%   Enumerates (h,k,l) over the given range, computes d-spacing from the
%   crystal system (xrdc.lattice.dSpacingFromHKL), and applies Bragg's
%   law to get 2θ. Duplicates (symmetry equivalents within 1e-6 Å of d)
%   are collapsed using the same tie-break the Delphi XRDC uses
%   (xrdc1.pas:4443+): first by min(h+k+l), then by min(l*10000+k*100+h).
%
%   Inputs
%     lattice  : struct with .system and lattice parameters
%                (see xrdc.lattice.dSpacingFromHKL).
%     hklRange : either
%                  - a scalar N (expands to h,k,l ∈ [-N, N])
%                  - [hmin hmax kmin kmax lmin lmax] (6-vector)
%     lambda   : wavelength in Å (scalar)
%
%   Name-Value options
%     KeepEquivalents (logical, default false)  — if true, keep all
%              symmetry-equivalent reflections instead of merging.
%     TwoThetaRange (1x2, default [0 180])  — restrict output by 2θ.
%
%   Output
%     T : table with columns h, k, l, d, twoTheta, label

    arguments
        lattice   (1,1) struct
        hklRange        double
        lambda    (1,1) double {mustBePositive}
        opts.KeepEquivalents (1,1) logical = false
        opts.TwoThetaRange   (1,2) double  = [0 180]
    end

    % expand hklRange
    if isscalar(hklRange)
        N = abs(round(hklRange));
        [H, K, L] = ndgrid(-N:N, -N:N, -N:N);
    elseif numel(hklRange) == 6
        r = round(hklRange);
        [H, K, L] = ndgrid(r(1):r(2), r(3):r(4), r(5):r(6));
    else
        error('xrdc:lattice:badRange', ...
            'hklRange must be a scalar or a 6-element vector.');
    end

    h = H(:); k = K(:); l = L(:);

    % drop (0,0,0)
    nonzero = ~(h == 0 & k == 0 & l == 0);
    h = h(nonzero); k = k(nonzero); l = l(nonzero);

    d = xrdc.lattice.dSpacingFromHKL(h, k, l, lattice);

    valid = isfinite(d) & d > 0;
    h = h(valid); k = k(valid); l = l(valid); d = d(valid);

    twoTheta = xrdc.lattice.dToTwoTheta(d, lambda);
    valid = isfinite(twoTheta) & twoTheta >= opts.TwoThetaRange(1) ...
                              & twoTheta <= opts.TwoThetaRange(2);
    h = h(valid); k = k(valid); l = l(valid);
    d = d(valid); twoTheta = twoTheta(valid);

    % sort by d descending
    [d, order] = sort(d, 'descend');
    h = h(order); k = k(order); l = l(order); twoTheta = twoTheta(order);

    % collapse duplicates within 1e-6 Å of d (matches Delphi XRDC)
    if ~opts.KeepEquivalents && ~isempty(d)
        keep = true(size(d));
        i = 1;
        while i <= numel(d)
            j = i + 1;
            while j <= numel(d) && abs(d(j) - d(i)) < 1e-6
                j = j + 1;
            end
            if j > i + 1
                % group is [i, j-1]; pick survivor by Delphi tie-break:
                %   1) min(h+k+l)   2) min(l*10000 + k*100 + h)
                groupIdx  = i:(j-1);
                sumHKL    = h(groupIdx) + k(groupIdx) + l(groupIdx);
                [~, primary] = min(sumHKL);
                primaries = groupIdx(sumHKL == sumHKL(primary));
                if numel(primaries) > 1
                    tieKey = l(primaries) * 10000 + k(primaries) * 100 + h(primaries);
                    [~, winIdx] = min(tieKey);
                    winner = primaries(winIdx);
                else
                    winner = primaries;
                end
                keep(groupIdx) = false;
                keep(winner)   = true;
            end
            i = j;
        end
        h = h(keep); k = k(keep); l = l(keep);
        d = d(keep); twoTheta = twoTheta(keep);
    end

    % build labels
    if strcmpi(lattice.system, 'hexagonal')
        % 4-index (h, -(h+k), k, l) matches Delphi's CalculateStructureLines
        labels = arrayfun(@(hh,kk,ll) ...
            sprintf('(%d %d %d %d)', hh, -(hh+kk), kk, ll), h, k, l, ...
            'UniformOutput', false);
    else
        labels = arrayfun(@(hh,kk,ll) sprintf('(%d %d %d)', hh, kk, ll), ...
            h, k, l, 'UniformOutput', false);
    end

    T = table(h, k, l, d, twoTheta, string(labels), ...
        'VariableNames', {'h','k','l','d','twoTheta','label'});
end
