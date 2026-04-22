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

        [~, hC] = contourf(ax, Qx, Qz, log(Ig), logContours, 'LineColor', 'none');
        plotHandle = struct('contour', hC, 'image', []);
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
    end

    % --- Colorbar with decade ticks ---------------------------------
    colormap(ax, options.Colormap);
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
    xlabel(ax, '{\itQ_x} (\AA^{-1})');
    ylabel(ax, '{\itQ_z} (\AA^{-1})');

    if ~isempty(options.AxesLim) && numel(options.AxesLim) == 4
        xlim(ax, options.AxesLim(1:2));
        ylim(ax, options.AxesLim(3:4));
    end

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
