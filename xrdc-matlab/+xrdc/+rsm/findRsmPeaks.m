function peaks = findRsmPeaks(kPar, kPerp, intensity, options)
%FINDRSMPEAKS  Locate substrate and film peaks in an RSM point cloud.
%   peaks = xrdc.rsm.findRsmPeaks(kPar, kPerp, intensity, Name=Value)
%
%   Material-independent. Bins the scattered (kPar,kPerp,intensity) cloud onto
%   a regular grid (max per cell), lightly smooths, and takes regional maxima.
%   Substrate = brightest regional max (single crystal -> brightest); film =
%   the next-brightest regional max at least 3 grid cells away, above
%   NoiseFactor*median(intensity). Each is refined by an intensity-weighted
%   centroid. Base MATLAB only (no Image Processing Toolbox).
%
%   Name/Value
%     GridBins      grid resolution per axis (default 200).
%     NoiseFactor   a maximum must exceed NoiseFactor*median(intensity)
%                   (default 5).
%     NearThreshold film flagged "filmNearSubstrate" if within this distance
%                   (k-units) of the substrate (default 0.03).
%     RefineWindow  half-width of the centroid window (k-units, default 0.006).
%
%   Output peaks (struct)
%     .substrate .kPar .kPerp .intensity .found
%     .film      .kPar .kPerp .intensity .found
%     .flags     string array: "filmNearSubstrate", "filmNotBrighter"
%
%   Flags (never silently fail): filmNearSubstrate for the near-degenerate
%   hard case (result returned, wants a human glance); filmNotBrighter when no
%   clear secondary maximum exists (film.found=false -- the substrate is never
%   returned as the film).

    arguments
        kPar      (:,1) double
        kPerp     (:,1) double
        intensity (:,1) double
        options.GridBins      (1,1) double {mustBeInteger, mustBePositive} = 200
        options.NoiseFactor   (1,1) double = 5
        options.NearThreshold (1,1) double = 0.03
        options.RefineWindow  (1,1) double = 0.006
    end

    nb = options.GridBins;
    bg = median(intensity);
    kpEdges = linspace(min(kPar),  max(kPar),  nb+1);
    kzEdges = linspace(min(kPerp), max(kPerp), nb+1);
    ip = discretize(kPar,  kpEdges);
    iz = discretize(kPerp, kzEdges);
    valid = ~isnan(ip) & ~isnan(iz);
    G = accumarray([iz(valid) ip(valid)], intensity(valid), [nb nb], @max, 0);

    % light Gaussian smoothing via conv2 (base MATLAB; NO imgaussfilt)
    [xx, yy] = meshgrid(-2:2, -2:2);
    ker = exp(-(xx.^2 + yy.^2)/2); ker = ker / sum(ker(:));
    Gs = conv2(G, ker, 'same');

    % interior regional maxima: cell >= all 8 neighbours (NO ordfilt2)
    cc = Gs(2:end-1, 2:end-1);
    ge = true(size(cc));
    for dr = -1:1
        for dc = -1:1
            if dr == 0 && dc == 0, continue, end
            ge = ge & cc >= Gs(2+dr:end-1+dr, 2+dc:end-1+dc);
        end
    end
    isMax = false(size(Gs));
    isMax(2:end-1, 2:end-1) = ge & cc > options.NoiseFactor*bg;

    kpC = (kpEdges(1:end-1) + kpEdges(2:end))/2;
    kzC = (kzEdges(1:end-1) + kzEdges(2:end))/2;
    dk  = hypot(kpEdges(2)-kpEdges(1), kzEdges(2)-kzEdges(1));
    flags = string.empty(1,0);

    [ri, ci] = find(isMax);
    if isempty(ri)                              % degenerate: fall back to global max
        [mx, im] = max(intensity);
        peaks.substrate = struct('kPar', kPar(im), 'kPerp', kPerp(im), ...
            'intensity', mx, 'found', true);
        peaks.film = struct('kPar', NaN, 'kPerp', NaN, 'intensity', NaN, 'found', false);
        peaks.flags = "filmNotBrighter";
        return
    end
    vals = Gs(sub2ind(size(Gs), ri, ci));
    [vals, ord] = sort(vals, 'descend'); ri = ri(ord); ci = ci(ord);

    subSeed = [kpC(ci(1)), kzC(ri(1))];
    sub = refineCentroid(kPar, kPerp, intensity, subSeed, options.RefineWindow);
    peaks.substrate = struct('kPar', sub(1), 'kPerp', sub(2), ...
        'intensity', vals(1), 'found', true);

    peaks.film = struct('kPar', NaN, 'kPerp', NaN, 'intensity', NaN, 'found', false);
    for m = 2:numel(ri)
        cand = [kpC(ci(m)), kzC(ri(m))];
        if hypot(cand(1)-subSeed(1), cand(2)-subSeed(2)) > 3*dk
            fc = refineCentroid(kPar, kPerp, intensity, cand, options.RefineWindow);
            peaks.film = struct('kPar', fc(1), 'kPerp', fc(2), ...
                'intensity', vals(m), 'found', true);
            if hypot(fc(1)-sub(1), fc(2)-sub(2)) < options.NearThreshold
                flags(end+1) = "filmNearSubstrate"; %#ok<AGROW>
            end
            break
        end
    end
    if ~peaks.film.found
        flags(end+1) = "filmNotBrighter";
    end
    peaks.flags = flags;
end

function c = refineCentroid(kPar, kPerp, intensity, seed, win)
%REFINECENTROID  Intensity-weighted centroid within +/- win of seed (k-units).
    sel = abs(kPar - seed(1)) <= win & abs(kPerp - seed(2)) <= win;
    w = intensity(sel);
    w = max(w - min(w), 0);               % local-background subtract
    if sum(w) <= 0, c = seed; return, end
    c = [sum(kPar(sel).*w), sum(kPerp(sel).*w)] / sum(w);
end
