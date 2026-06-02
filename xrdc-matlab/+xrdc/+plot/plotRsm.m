function h = plotRsm(scans, options)
%PLOTRSM  Publication-quality reciprocal-space map.
%   h = xrdc.plot.plotRsm(scans)
%   h = xrdc.plot.plotRsm(scans, Name, Value, ...)
%
%   Matches the style of Schwaigert et al. JVST A 41, 022703 (2023),
%   Fig 2(e): filled contour on log intensity, decade-tick colorbar
%   (1, 10, 10^2, ..., 10^N), Arial 18 pt, Painters renderer.
%
%   For each slice calls xrdc.rsm.toReciprocalSpace to get (k_par, k_perp);
%   if every slice has the same number of 2θ points, treats the data as a
%   2-D grid and uses contourf directly (best quality, no interpolation).
%   Otherwise interpolates onto a uniform grid via scatteredInterpolant.
%
%   Input
%     scans : (1×N) struct array from xrdc.rsm.loadAreaScan.
%
%   Name/Value
%     'Lambda'         (1,1) double  — wavelength in Å (overrides scan.lambda)
%     'DeltaTheta'     (1,1) double  — 2θ offset, degrees  (default 0)
%     'DeltaOmega'     (1,1) double  — ω offset, degrees   (default 0)
%     'Flip'           (1,1) logical — negate k_par         (default false)
%     'Mode'           (1,1) string  — "contourf" | "imagesc" (default "contourf")
%     'NContours'      (1,1) double  — number of log-spaced contour levels (default 30)
%     'Imin'           (1,1) double  — min intensity for colorbar (default auto)
%     'Imax'           (1,1) double  — max intensity for colorbar (default auto)
%     'AxesLim'        (1,4) double  — [xmin xmax ymin ymax] (default auto = [])
%     'GridN'          (1,2) double  — interpolation grid for "imagesc" mode
%                                      (default [512 512])
%     'Colormap'       (1,1) string  — 'turbo' | 'parula' | 'jet' (default 'turbo')
%     'Smooth'         (1,1) logical — 3×3 mean filter before plotting (default false)
%     'ExportPath'     (1,1) string  — if non-empty, export PNG at 600 dpi
%     'TargetAxes'                  = []
%     'FontName'       (1,1) string  = "Arial"
%     'TickFontSize'   (1,1) double  = 18
%     'LabelFontSize'  (1,1) double  = 20
%     'TitleFontSize'  (1,1) double  = 22
%
%   Output h — struct with .contour / .image, .ax, .figure, .colorbar.

    arguments
        scans                          (1,:) struct
        options.Lambda                 (1,1) double  = NaN
        options.DeltaTheta             (1,1) double  = 0
        options.DeltaOmega             (1,1) double  = 0
        options.Flip                   (1,1) logical = false
        options.Mode                   (1,1) string  = "contourf"
        options.NContours              (1,1) double  = 30
        options.Imin                   (1,1) double  = NaN
        options.Imax                   (1,1) double  = NaN
        options.AxesLim                              = []
        options.GridN                  (1,2) double  = [512 512]
        options.Colormap               (1,1) string  = "turbo"
        options.AutoZoom               (1,1) logical = true
        options.Smooth                 (1,1) logical = false
        options.ExportPath             (1,1) string  = ""
        options.TargetAxes                           = []
        options.FontName               (1,1) string  = "Arial"
        options.TickFontSize           (1,1) double  = 18
        options.LabelFontSize          (1,1) double  = 20
        options.TitleFontSize          (1,1) double  = 22
    end

    if isempty(scans)
        error('xrdc:plot:emptyScans', 'scans array is empty.');
    end

    % --- Transform each slice to (kPar, kPerp) -----------------------
    nS = numel(scans);
    kParCells  = cell(1, nS);
    kPerpCells = cell(1, nS);
    intCells   = cell(1, nS);
    for i = 1:nS
        [kP, kZ] = xrdc.rsm.toReciprocalSpace(scans(i), ...
            'Lambda',      options.Lambda, ...
            'DeltaTheta',  options.DeltaTheta, ...
            'DeltaOmega',  options.DeltaOmega, ...
            'Flip',        options.Flip);
        kParCells{i}  = kP(:);
        kPerpCells{i} = kZ(:);
        intCells{i}   = double(scans(i).counts(:));
    end

    % --- Figure / axes setup -----------------------------------------
    if isempty(options.TargetAxes)
        fig = figure('Renderer', 'Painters');
        ax  = axes(fig);
    else
        ax  = options.TargetAxes;
        fig = ancestor(ax, 'figure');
        set(fig, 'Renderer', 'Painters');
    end
    hold(ax, 'on');

    % --- Check for uniform grid -> contourf direct ------------------
    nPts = cellfun(@numel, intCells);
    canGrid = all(nPts == nPts(1));
    useContour = options.Mode == "contourf";

    if canGrid && useContour
        % Each slice becomes one column; rows index 2θ within the slice
        Qx = cell2mat(kParCells);   % [nTT × nS]
        Qz = cell2mat(kPerpCells);  % [nTT × nS]
        Ig = cell2mat(intCells);    % [nTT × nS]
        Ig(Ig <= 0) = 1;            % log(0) guard, mirrors Barone RSMPlot.m

        if options.Smooth
            Ig = smooth3x3(Ig);
        end

        Imin = pickDefault(options.Imin, min(Ig, [], 'all'));
        Imax = pickDefault(options.Imax, max(Ig, [], 'all'));
        logContours = log(logspace(log10(Imin), log10(Imax), options.NContours));

        % Mask values below Imin so contourf leaves them as the axes
        % background (white/light) rather than the lowest colormap color.
        % Mirrors the JVST A 2023 Fig 2(e) appearance where weak/noise
        % regions appear unfilled rather than dark blue.
        logIg = log(Ig);
        logIg(Ig < Imin) = NaN;

        [~, hC] = contourf(ax, Qx, Qz, logIg, logContours, 'LineColor', 'none');
        plotHandle = struct('contour', hC, 'image', []);

        % Store for AutoZoom downstream
        zoomQx = Qx; zoomQz = Qz; zoomIg = Ig;
    else
        % --- Fall back: scatter → interpolate → imagesc -------------
        kParAll  = vertcat(kParCells{:});
        kPerpAll = vertcat(kPerpCells{:});
        intAll   = vertcat(intCells{:});
        intAll(intAll <= 0) = 1;

        Imin = pickDefault(options.Imin, min(intAll));
        Imax = pickDefault(options.Imax, max(intAll));

        kParEdge  = linspace(min(kParAll),  max(kParAll),  options.GridN(1));
        kPerpEdge = linspace(min(kPerpAll), max(kPerpAll), options.GridN(2));
        [PGrid, ZGrid] = meshgrid(kParEdge, kPerpEdge);
        F = scatteredInterpolant(kParAll, kPerpAll, log(intAll), 'linear', 'none');
        zGrid = F(PGrid, ZGrid);

        hI = imagesc(ax, kParEdge, kPerpEdge, zGrid);
        set(ax, 'YDir', 'normal');
        plotHandle = struct('contour', [], 'image', hI);

        % Store for AutoZoom: use the scatter data (more faithful to
        % the actual peak locations than the interpolated grid).
        zoomQx = kParAll; zoomQz = kPerpAll; zoomIg = intAll;
    end

    % --- Colorbar with decade ticks ---------------------------------
    colormap(ax, options.Colormap);
    % Guard: caxis requires strictly increasing limits. Uniform data (Imin==Imax)
    % occurs in synthetic tests and degenerate plots — pad by one decade.
    if ~(Imax > Imin)
        Imax = max(Imin * 10, Imin + 1);
    end
    [ticks, tickLabels] = decadeTicks(Imin, Imax);
    cb = colorbar(ax, 'Ticks', ticks, 'TickLabels', tickLabels);
    cb.Label.String   = 'Intensity (counts)';
    cb.Label.FontSize = options.LabelFontSize;
    cb.Label.FontName = char(options.FontName);
    caxis(ax, log([Imin Imax]));

    % --- Publication style ------------------------------------------
    xrdc.plot.publicationStyle(ax, ...
        'FontName',      options.FontName, ...
        'TickFontSize',  options.TickFontSize, ...
        'LabelFontSize', options.LabelFontSize, ...
        'TitleFontSize', options.TitleFontSize, ...
        'LogY',          false);
    set(ax, 'LineWidth', 1.5);
    % Use Unicode Å since MATLAB's default TeX interpreter does not accept \AA.
    xlabel(ax, '{\itQ}_x (Å^{-1})');
    ylabel(ax, '{\itQ}_z (Å^{-1})');

    if ~isempty(options.AxesLim) && numel(options.AxesLim) == 4
        xlim(ax, options.AxesLim(1:2));
        ylim(ax, options.AxesLim(3:4));
    elseif options.AutoZoom
        [xl, yl] = autoZoomBounds(zoomQx, zoomQz, zoomIg, Imin);
        xlim(ax, xl);
        ylim(ax, yl);
    end

    % Equal data aspect ratio: 1 Å⁻¹ in Qx covers the same screen distance
    % as 1 Å⁻¹ in Qz. This is the standard for reciprocal-space maps —
    % preserves angles and distances in q-space and matches the
    % JVST A 2023 Fig 2(e) appearance (tall, narrow window).
    daspect(ax, [1 1 1]);

    hold(ax, 'off');

    % --- Export ------------------------------------------------------
    if strlength(options.ExportPath) > 0
        exportgraphics(fig, char(options.ExportPath), 'Resolution', 600);
    end

    h = struct( ...
        'contour',  plotHandle.contour, ...
        'image',    plotHandle.image, ...
        'ax',       ax, ...
        'figure',   fig, ...
        'colorbar', cb);
end

% =====================================================================
function v = pickDefault(opt, def)
    if isnan(opt), v = def; else, v = opt; end
end

function M = smooth3x3(M)
    kernel = ones(3) / 9;
    M = conv2(M, kernel, 'same');
end

function [xl, yl] = autoZoomBounds(Qx, Qz, Ig, Imin)
%AUTOZOOMBOUNDS  Bounding box around the visible features of an RSM.
%   Thresholds in log-intensity space (captures streaks and weak peaks,
%   not just the bright core), pads by 40 %, and enforces a minimum
%   window of 15 % of the full data span so the bbox never collapses
%   onto a few pixels.
    qx = Qx(:); qz = Qz(:); ig = Ig(:);
    finite = isfinite(qx) & isfinite(qz) & isfinite(ig);
    qx = qx(finite); qz = qz(finite); ig = ig(finite);

    fullSpanX = max(qx) - min(qx);
    fullSpanZ = max(qz) - min(qz);

    % Log-intensity threshold at 30 % of the dynamic range above Imin.
    % Cu Kα RSMs typically span ~5 decades from noise floor to substrate
    % peak; the streak/CTR features sit at ~10²–10³ counts and must be
    % retained, not just the 10⁵-count core.
    logImin = log10(max(Imin, 1));
    logMax  = log10(max(max(ig), Imin * 10));
    logThr  = logImin + 0.30 * (logMax - logImin);
    thr     = 10 ^ logThr;
    mask    = ig >= thr;

    if nnz(mask) < 5
        % Fallback: centre on the brightest point, span = 20 % of data range
        [~, idx] = max(ig);
        cx = qx(idx); cz = qz(idx);
        xl = [cx - fullSpanX * 0.10, cx + fullSpanX * 0.10];
        yl = [cz - fullSpanZ * 0.10, cz + fullSpanZ * 0.10];
        return
    end

    qxB = qx(mask); qzB = qz(mask);
    xl = [min(qxB), max(qxB)];
    yl = [min(qzB), max(qzB)];

    % 40 % padding on each side
    xspan = max(diff(xl), eps);
    zspan = max(diff(yl), eps);
    xl = xl + [-1 1] * xspan * 0.40;
    yl = yl + [-1 1] * zspan * 0.40;

    % Enforce a minimum window of 15 % of full data range so we never
    % crop down to a sliver when only a few pixels exceed threshold.
    minSpanX = fullSpanX * 0.15;
    minSpanZ = fullSpanZ * 0.15;
    if diff(xl) < minSpanX
        cx = mean(xl);
        xl = [cx - minSpanX/2, cx + minSpanX/2];
    end
    if diff(yl) < minSpanZ
        cz = mean(yl);
        yl = [cz - minSpanZ/2, cz + minSpanZ/2];
    end
end

function [ticks, labels] = decadeTicks(Imin, Imax)
    %DECADETICKS  Build colorbar ticks at log(1), log(10), log(10^2)…
    %   Matches Barone RSMPlot.m colorbar style (JVST A 2023 Fig 2(e)).
    eLo = floor(log10(max(Imin, 1)));
    eHi = ceil (log10(Imax));
    exps = eLo:eHi;
    ticks = log(10 .^ exps);
    labels = cell(size(exps));
    for k = 1:numel(exps)
        if exps(k) == 0
            labels{k} = '1';
        elseif exps(k) == 1
            labels{k} = '10';
        else
            labels{k} = sprintf('10^{%d}', exps(k));
        end
    end
end
