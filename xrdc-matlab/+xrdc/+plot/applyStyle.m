function applyStyle(ax, style)
%APPLYSTYLE  Apply user plot-style overrides on top of analysis defaults.
%   xrdc.plot.applyStyle(ax, style)
%
%   The GUI's per-scan-type analysis sets sensible default title, labels,
%   axis scale, colours, etc. This applies the user's optional overrides on
%   TOP of those defaults so people can tune a figure for publication
%   without losing the defaults. Every field is opt-in: a blank text field
%   or an "auto" choice leaves whatever the analysis already set.
%
%   Input
%     ax    : target axes (uiaxes or axes).
%     style : struct of overrides. Recognised fields (all optional):
%       .title .xlabel .ylabel        char/string; ''  → leave as-is
%       .xmin .xmax .ymin .ymax       char/string number; '' → leave
%       .yscale                       "auto" | "linear" | "log"
%       .fontSize                     char/string number; '' → leave
%       .lineWidth                    char/string number; '' → leave
%       .lineColor                    "auto" | colour name (see colourFromName)
%       .markers                      "auto" | "on" | "off"
%       .grid                         "auto" | "on" | "off"
%
%   Notes
%     - The data trace is taken to be the line with the most points; peak
%       markers are the lines drawn with a marker. RSM (contour) plots have
%       no line, so line/marker options are simply skipped there.

    arguments
        ax    (1,1)
        style (1,1) struct
    end
    if ~isvalid(ax)
        return
    end

    if hasText(style, 'title'),  title(ax,  string(style.title));  end
    if hasText(style, 'xlabel'), xlabel(ax, string(style.xlabel)); end
    if hasText(style, 'ylabel'), ylabel(ax, string(style.ylabel)); end

    setLim(ax, 'X', field(style, 'xmin'), field(style, 'xmax'));
    setLim(ax, 'Y', field(style, 'ymin'), field(style, 'ymax'));

    ys = lower(field(style, 'yscale', 'auto'));
    if any(strcmp(ys, {'linear', 'log'}))
        set(ax, 'YScale', ys);
    end

    fs = parseNum(field(style, 'fontSize'));
    if ~isnan(fs) && fs > 0
        ax.FontSize = fs;
        if ~isempty(ax.Title),  ax.Title.FontSize  = fs + 1; end
        if ~isempty(ax.XLabel), ax.XLabel.FontSize = fs;     end
        if ~isempty(ax.YLabel), ax.YLabel.FontSize = fs;     end
    end

    switch lower(field(style, 'grid', 'auto'))
        case 'on',  grid(ax, 'on');
        case 'off', grid(ax, 'off');
    end

    lns = findobj(ax, 'Type', 'line');
    if ~isempty(lns)
        npts   = arrayfun(@(h) numel(h.XData), lns);
        [~, it] = max(npts);
        trace  = lns(it);

        lw = parseNum(field(style, 'lineWidth'));
        if ~isnan(lw) && lw > 0
            trace.LineWidth = lw;
        end
        lc = lower(field(style, 'lineColor', 'auto'));
        if ~strcmp(lc, 'auto')
            trace.Color = colourFromName(lc);
        end

        mk = lns(arrayfun(@(h) ~strcmp(h.Marker, 'none'), lns));
        if ~isempty(mk)
            switch lower(field(style, 'markers', 'auto'))
                case 'off', set(mk, 'Visible', 'off');
                case 'on',  set(mk, 'Visible', 'on');
            end
        end
    end
end

% ---------------------------------------------------------------------

function tf = hasText(s, f)
    tf = isfield(s, f) && strlength(strtrim(string(s.(f)))) > 0;
end

function v = field(s, f, default)
    if nargin < 3, default = ''; end
    if isfield(s, f), v = s.(f); else, v = default; end
    v = char(string(v));
end

function setLim(ax, axChar, mn, mx)
    lo = parseNum(mn); hi = parseNum(mx);
    if isnan(lo) && isnan(hi)
        return
    end
    cur = get(ax, [axChar 'Lim']);
    if ~isnan(lo), cur(1) = lo; end
    if ~isnan(hi), cur(2) = hi; end
    if cur(2) > cur(1)
        set(ax, [axChar 'Lim'], cur);
    end
end

function v = parseNum(x)
    if isnumeric(x)
        v = double(x);
        return
    end
    s = strtrim(string(x));
    if strlength(s) == 0
        v = NaN;
    else
        v = str2double(s);
    end
end

function c = colourFromName(name)
    switch lower(string(name))
        case "blue",   c = [0.10 0.40 0.80];
        case "black",  c = [0 0 0];
        case "red",    c = [0.85 0.20 0.20];
        case "green",  c = [0.20 0.60 0.20];
        case "orange", c = [0.90 0.50 0.10];
        case "purple", c = [0.50 0.20 0.60];
        case "gray",   c = [0.40 0.40 0.40];
        otherwise,     c = [0.10 0.40 0.80];
    end
end
