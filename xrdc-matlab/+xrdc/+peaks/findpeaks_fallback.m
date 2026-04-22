function [pks, locs, widths, proms] = findpeaks_fallback(y, varargin)
%FINDPEAKS_FALLBACK  Peak detection without Signal Processing Toolbox.
%   [pks, locs, widths, proms] = findpeaks_fallback(y, Name, Value, ...)
%
%   Detects peaks in y by finding local maxima and computing prominence.
%   Supports key options from MATLAB's findpeaks:
%     MinPeakProminence, MinPeakHeight, MinPeakDistance, MinPeakWidth, MaxPeakWidth
%
%   Returns:
%     pks: peak values
%     locs: peak indices (into y)
%     widths: peak FWHM (in samples)
%     proms: peak prominence

    % Parse options
    p = inputParser;
    addParameter(p, 'MinPeakProminence', 0);
    addParameter(p, 'MinPeakHeight', -Inf);
    addParameter(p, 'MinPeakDistance', 1);
    addParameter(p, 'MinPeakWidth', 0);
    addParameter(p, 'MaxPeakWidth', Inf);
    addParameter(p, 'WidthReference', 'halfheight');
    parse(p, varargin{:});
    opts = p.Results;

    y = y(:);
    n = numel(y);

    % Find local maxima (candidate peaks)
    if n < 3
        pks = []; locs = []; widths = []; proms = [];
        return
    end

    % Local maximum: y(i) > y(i-1) and y(i) > y(i+1)
    is_peak = [false; y(2:end-1) > y(1:end-2) & y(2:end-1) > y(3:end); false];
    locs = find(is_peak);

    if isempty(locs)
        pks = []; widths = []; proms = [];
        return
    end

    pks = y(locs);

    % Compute prominence (height above the lowest contour)
    proms = zeros(size(locs));
    widths = zeros(size(locs));

    for i = 1:numel(locs)
        loc = locs(i);
        pk_val = pks(i);

        % Left reference: lowest point between this peak and the previous peak
        if loc == 1
            left_ref = y(1);
        else
            left_idx = loc - 1;
            while left_idx > 1 && y(left_idx) < pk_val
                left_idx = left_idx - 1;
            end
            left_ref = min(y(left_idx:loc));
        end

        % Right reference: lowest point between this peak and the next peak
        if loc == n
            right_ref = y(n);
        else
            right_idx = loc + 1;
            while right_idx < n && y(right_idx) < pk_val
                right_idx = right_idx + 1;
            end
            right_ref = min(y(loc:right_idx));
        end

        % Prominence is height above the lower of the two contours
        contour = max(left_ref, right_ref);
        proms(i) = max(0, pk_val - contour);

        % FWHM: find left and right where y crosses half-height
        half_height = (pk_val + contour) / 2;

        % Find left half-max crossing
        left_half_idx = loc;
        for j = loc-1:-1:1
            if y(j) < half_height
                left_half_idx = j;
                break
            end
        end

        % Find right half-max crossing
        right_half_idx = loc;
        for j = loc+1:n
            if y(j) < half_height
                right_half_idx = j;
                break
            end
        end

        widths(i) = right_half_idx - left_half_idx;
    end

    % Apply filters
    keep = true(size(locs));
    keep = keep & pks >= opts.MinPeakHeight;
    keep = keep & proms >= opts.MinPeakProminence;

    if opts.MinPeakWidth > 0
        keep = keep & widths >= opts.MinPeakWidth;
    end
    if isfinite(opts.MaxPeakWidth)
        keep = keep & widths <= opts.MaxPeakWidth;
    end

    % MinPeakDistance: remove peaks that are too close together
    if opts.MinPeakDistance > 1
        to_remove = false(size(locs));
        for i = 1:numel(locs)
            if to_remove(i)
                continue
            end
            % Remove nearby lower peaks
            for j = i+1:numel(locs)
                if to_remove(j)
                    continue
                end
                if locs(j) - locs(i) < opts.MinPeakDistance
                    if pks(j) > pks(i)
                        to_remove(i) = true;
                    else
                        to_remove(j) = true;
                    end
                end
            end
        end
        keep = keep & ~to_remove;
    end

    % Return filtered results
    locs = locs(keep);
    pks = pks(keep);
    widths = widths(keep);
    proms = proms(keep);
end
