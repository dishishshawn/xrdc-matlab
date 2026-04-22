function scan = readScan(path)
%READSCAN  Format-sniffing dispatcher for XRD scan files.
%   scan = xrdc.io.readScan(path)
%
%   Dispatches to the right reader based on file contents (first non-
%   empty line) and extension. Mirrors ReadFile in xrdc1.pas:438, plus
%   Rigaku formats (not yet implemented — see README for the blocker).
%
%   Input
%     path : path to file (char or string)
%
%   Output
%     scan : struct in xrdc.io.emptyScan() shape

    arguments
        path (1,1) string
    end

    if ~isfile(path)
        error('xrdc:io:notFound', 'File not found: %s', path);
    end

    % Peek at the first ~200 bytes to sniff format
    fid = fopen(path, 'r');
    cleanup = onCleanup(@() fclose(fid));
    head = fread(fid, 512, '*char')';
    clear cleanup

    headUpper = upper(strtrim(head));

    if startsWith(headUpper, "<?XML") || contains(headUpper, "<XRDMEASUREMENT")
        scan = xrdc.io.readXrdml(path);
    elseif contains(headUpper, "HR-XRDSCAN")
        scan = xrdc.io.readPhilipsX00(path);
    elseif startsWith(headUpper, "*RAS_DATA_START") || contains(headUpper, "*MEAS_COND")
        error('xrdc:io:notImplemented', ...
            'Rigaku .ras parser not yet implemented. Send sample files to unblock Phase 1.');
    elseif startsWith(headUpper, "RAW1.0")
        error('xrdc:io:notImplemented', ...
            'Rigaku binary .raw parser not yet implemented. Send sample files to unblock Phase 1.');
    elseif startsWith(headUpper, "SUPERDATA")
        error('xrdc:io:notSupported', ...
            'Picker .596/.1035 format is not ported (obsolete).');
    elseif isRigakuTxt(path, headUpper)
        scan = xrdc.io.readRigakuTxt(path);
    else
        scan = xrdc.io.readTextScan(path);
    end
end

% =====================================================================
function tf = isRigakuTxt(path, headUpper)
    %ISRIGAKUTXT  Heuristic: does this look like a Rigaku SmartLab .txt export?
    %
    %   Markers:
    %     - file extension .txt
    %     - header contains "2Θ" (U+03B8) or the byte-sequence "2θ,"
    %       appearing in Rigaku's ASCII column labels
    %     - OR filename contains "TR_" prefix (Rigaku SmartLab sample naming)
    [~, stem, ext] = fileparts(path);
    ext = lower(string(ext));
    if ext ~= ".txt"
        tf = false;
        return;
    end
    % Header label cue — "INTENSITY, CPS" is Rigaku's unique column label
    if contains(headUpper, "INTENSITY, CPS")
        tf = true;
        return;
    end
    % Filename cue (TR_ prefix commonly used by the Paik lab's Rigaku)
    if startsWith(upper(string(stem)), "TR_")
        tf = true;
        return;
    end
    tf = false;
end
