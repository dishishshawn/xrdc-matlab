function scan = readRigakuRas(path)  %#ok<STOUT,INUSD>
%READRIGAKURAS  Read a Rigaku RAS (ASCII) scan — NOT IMPLEMENTED.
%
%   Blocked on sample files from Dr. Paik. Target spec is documented in
%   ../../docs/RIGAKU_NOTES.md once that exists. The parser needs to
%   handle:
%     - Header keys prefixed *MEAS_COND_*
%     - Data block between *RAS_INT_START and *RAS_INT_END
%     - Mapping MEAS_COND_AXIS_NAME to xrdc scanType
%
%   See ALGORITHM_SPEC.md §2.5 for the full target contract.

    error('xrdc:io:notImplemented', ...
        ['Rigaku .ras parser not yet implemented. ', ...
         'Send 3–5 sample .ras files to unblock.']);
end
