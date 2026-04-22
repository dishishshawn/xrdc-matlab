function scans = loadAreaScan(source, options)
%LOADAREASCAN  Load a set of θ-2θ slices for reciprocal-space mapping.
%   scans = xrdc.rsm.loadAreaScan(folderPath)
%   scans = xrdc.rsm.loadAreaScan(folderPath, 'Pattern', '*.xrdml')
%   scans = xrdc.rsm.loadAreaScan(fileList)
%
%   Reads every file matched by folderPath+Pattern (or the explicit list),
%   dispatches each through xrdc.io.readScan, and returns a struct array
%   sorted by secondAxis ascending.
%
%   Input
%     source : string  — path to a folder, OR
%              string array / cell array of file paths
%
%   Name/Value
%     'Pattern'  (1,1) string  — glob pattern used when source is a folder
%                               (default "*" = everything readScan accepts)
%     'Lambda'   (1,1) double  — override wavelength for all scans (Å)
%                               (default NaN = use each file's own lambda)
%
%   Output
%     scans : (1×N) struct array, each element as per xrdc.io.emptyScan,
%             sorted by .secondAxis ascending.  .scanType is set to "area"
%             on every element so callers can distinguish area-scan members.

    arguments
        source                  (1,:)           % string, string array, or cell
        options.Pattern         (1,1) string  = "*"
        options.Lambda          (1,1) double  = NaN
    end

    % Collect file paths
    if ischar(source) || (isstring(source) && isscalar(source))
        % Folder or single string — treat as folder + pattern
        folder = char(source);
        if ~isfolder(folder)
            error('xrdc:rsm:notAFolder', ...
                'Source "%s" is not a folder. To load one file, pass a 1-element string array.', folder);
        end
        listing = dir(fullfile(folder, char(options.Pattern)));
        listing = listing(~[listing.isdir]);
        if isempty(listing)
            error('xrdc:rsm:noFiles', ...
                'No files matching pattern "%s" found in "%s".', options.Pattern, folder);
        end
        filePaths = fullfile({listing.folder}, {listing.name});
    elseif iscell(source) || (isstring(source) && ~isscalar(source))
        filePaths = cellstr(source);
    else
        error('xrdc:rsm:badSource', ...
            'source must be a folder path (string) or a list of file paths.');
    end

    nFiles = numel(filePaths);
    scans  = repmat(xrdc.io.emptyScan(), 1, nFiles);

    for i = 1:nFiles
        s = xrdc.io.readScan(filePaths{i});
        s.scanType = "area";      % mark as area-scan member
        if ~isnan(options.Lambda)
            s.lambda = options.Lambda;
        end
        scans(i) = s;
    end

    % Sort by secondAxis ascending (handles unsorted folder listings)
    [~, idx] = sort([scans.secondAxis]);
    scans = scans(idx);
end
