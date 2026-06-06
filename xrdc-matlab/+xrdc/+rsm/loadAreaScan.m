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
%     'UseParallel' (1,1) logical — read files in parallel via parfor
%                               (default true). Requires the Parallel
%                               Computing Toolbox; silently falls back to
%                               a serial for-loop otherwise.
%
%   Output
%     scans : (1×N) struct array, each element as per xrdc.io.emptyScan,
%             sorted by .secondAxis ascending.  .scanType is set to "area"
%             on every element so callers can distinguish area-scan members.

    arguments
        source                  (1,:)           % string, string array, or cell
        options.Pattern         (1,1) string  = "*"
        options.Lambda          (1,1) double  = NaN
        options.UseParallel     (1,1) logical = true
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

    nFiles    = numel(filePaths);
    collected = cell(1, nFiles);
    lambdaOv  = options.Lambda;
    useParfor = options.UseParallel && nFiles > 1 && hasParallelToolbox();

    if useParfor
        parfor i = 1:nFiles
            collected{i} = readOneFile(filePaths{i}, lambdaOv);
        end
    else
        for i = 1:nFiles
            collected{i} = readOneFile(filePaths{i}, lambdaOv);
        end
    end

    scans = [collected{:}];

    % Sort by secondAxis ascending (handles unsorted folder listings)
    [~, idx] = sort([scans.secondAxis]);
    scans = scans(idx);
end

function fileScans = readOneFile(fp, lambdaOv)
%READONEFILE  Dispatch a single area-scan file through the right reader.
    [~, ~, ext] = fileparts(fp);
    if strcmpi(ext, '.xrdml') && isMultiScanXrdml(fp)
        fileScans = xrdc.io.readXrdmlArea(fp);
    else
        fileScans = xrdc.io.readScan(fp);
    end
    for k = 1:numel(fileScans)
        fileScans(k).scanType = "area";
        if ~isnan(lambdaOv)
            fileScans(k).lambda = lambdaOv;
        end
    end
end

function tf = hasParallelToolbox()
%HASPARALLELTOOLBOX  True if the Parallel Computing Toolbox is licensed.
    persistent cached
    if isempty(cached)
        cached = ~isempty(ver('parallel')) && license('test', 'Distrib_Computing_Toolbox');
    end
    tf = cached;
end

function tf = isMultiScanXrdml(path)
%ISMULTISCANXRDML  Quickly decide whether an XRDML file has multiple <scan> blocks.
%   Streams the first ~2 MB looking for a second <scan> tag so large files
%   don't need full XML parsing just for dispatch.
    tf = false;
    fid = fopen(path, 'r');
    if fid < 0
        return
    end
    cleanup = onCleanup(@() fclose(fid));
    count = 0;
    while ~feof(fid)
        chunk = fread(fid, 1e6, '*char')';
        if isempty(chunk)
            break
        end
        count = count + numel(regexp(chunk, '<scan\s', 'start'));
        if count >= 2
            tf = true;
            return
        end
    end
end
