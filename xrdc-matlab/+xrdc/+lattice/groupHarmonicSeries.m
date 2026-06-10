function [S, singles] = groupHarmonicSeries(d, opts)
%GROUPHARMONICSERIES  Cluster d-spacings into harmonic (00l) series.
%   [S, singles] = xrdc.lattice.groupHarmonicSeries(d)
%
%   In a (00l)-oriented scan every phase gives d(00l) = c/l. Greedy
%   clustering: every (peak, assumed order 1..MaxOrder) pair seeds a
%   candidate c; peaks whose nearest-integer order reproduces their d
%   within Tolerance join; the largest (then tightest) series is taken
%   and the loop repeats on the remainder. One peak per order (closest
%   wins). The greedy pass prefers the smallest c that explains a set —
%   order-doubling disambiguation belongs to the database-matching stage
%   (see identifyMaterial).
%
%   Inputs
%     d : vector of d-spacings (A)
%
%   Name-Value
%     Tolerance (default 0.005) — max relative |d - c/l| / (c/l)
%     MaxOrder  (default 4)
%     CRange    (default [2 8]) — plausible c window (A) for seeds
%
%   Outputs
%     S       : struct array, fields .members (indices into d), .orders,
%               .c (least-squares), .cSigma, .residRms (relative),
%               .evenOnly (all matched orders even)
%     singles : indices of peaks in no multi-peak series
%
%   See also xrdc.lattice.identifyMaterial.

    arguments
        d (:,1) double {mustBePositive}
        opts.Tolerance (1,1) double {mustBePositive} = 0.005
        opts.MaxOrder  (1,1) double {mustBeInteger, mustBePositive} = 4
        opts.CRange    (1,2) double = [2 8]
    end

    n = numel(d);
    unassigned = true(n, 1);
    S = struct('members', {}, 'orders', {}, 'c', {}, 'cSigma', {}, ...
               'residRms', {}, 'evenOnly', {});

    while true
        best = struct('members', [], 'orders', [], 'rms', Inf);
        idx = find(unassigned).';
        for j = idx
            for l = 1:opts.MaxOrder
                c0 = d(j) * l;
                if c0 < opts.CRange(1) || c0 > opts.CRange(2), continue, end
                [mem, ord, rms] = collectMembers(d, idx, c0, opts);
                better = numel(mem) > numel(best.members) || ...
                        (numel(mem) == numel(best.members) && rms < best.rms);
                if numel(mem) >= 2 && better
                    best = struct('members', mem, 'orders', ord, 'rms', rms);
                end
            end
        end
        if numel(best.members) < 2, break, end

        cVals  = d(best.members) .* best.orders(:);
        c      = mean(cVals);
        sigma  = std(cVals) / sqrt(numel(cVals));
        resid  = (d(best.members) - c ./ best.orders(:)) ./ (c ./ best.orders(:));
        S(end+1) = struct('members', best.members, 'orders', best.orders, ...
            'c', c, 'cSigma', sigma, 'residRms', sqrt(mean(resid.^2)), ...
            'evenOnly', all(mod(best.orders, 2) == 0)); %#ok<AGROW>
        unassigned(best.members) = false;
    end

    singles = find(unassigned);
end

function [mem, ord, rms] = collectMembers(d, idx, c0, opts)
%COLLECTMEMBERS  Peaks consistent with series constant c0; one per order.
    mem = []; ord = []; resid = [];
    for k = idx
        lk = round(c0 / d(k));
        if lk < 1 || lk > opts.MaxOrder, continue, end
        dPred = c0 / lk;
        r = abs(d(k) - dPred) / dPred;
        if r > opts.Tolerance, continue, end
        prev = find(ord == lk, 1);
        if isempty(prev)
            mem(end+1) = k; ord(end+1) = lk; resid(end+1) = r; %#ok<AGROW>
        elseif r < resid(prev)          % closer claimant wins the order
            mem(prev) = k; resid(prev) = r;
        end
    end
    rms = sqrt(mean(resid.^2));
    if isempty(mem), rms = Inf; end
end
