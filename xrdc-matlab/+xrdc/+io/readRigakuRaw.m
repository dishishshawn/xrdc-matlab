function scan = readRigakuRaw(path)  %#ok<STOUT,INUSD>
%READRIGAKURAW  Read a Rigaku binary .raw scan — NOT IMPLEMENTED.
%
%   Blocked on sample files from Dr. Paik. Binary format variants:
%     - RAW1.01 / RAW1.02 (older SmartLab)
%     - RAW4.x (newer)
%
%   See ALGORITHM_SPEC.md §2.5 for the target contract, and
%   docs/RIGAKU_NOTES.md (once it exists) for reverse-engineering notes.

    error('xrdc:io:notImplemented', ...
        ['Rigaku binary .raw parser not yet implemented. ', ...
         'Send 3–5 sample .raw files to unblock.']);
end
