function cfg = defaults()
%DEFAULTS  Canonical lab-preset settings for xrdcApp.
%   cfg = xrdc.config.defaults()
%
%   Returns the factory-default settings struct that the app's "Settings…"
%   menu edits and persists (see xrdc.config.load / xrdc.config.save). The
%   values here reproduce the historical hardcoded behaviour of the Paik
%   lab's Rigaku, so an unconfigured install behaves exactly as before.
%
%   Fields
%     wavelength       (double, Å) — λ applied to scans whose file format
%                       drops it (Rigaku .txt/.hgx, plain text). Cu Kα1 by
%                       default. PANalytical .xrdml / Philips .x00 carry a
%                       real λ and are never overridden (see applyWavelength).
%     applyWavelength  (logical)   — master switch for the λ override.
%     substrate        (string)    — default substrate for material ID and
%                       strain-RSM analysis (was hardcoded "SrTiO3").
%     rcShape          (string)    — rocking-curve default fit shape:
%                       "gauss" | "lorentz" | "pseudoVoigt".
%     useFilenameRules (logical)   — honour filename cues (TR_ prefix, RC /
%                       XRR / 2theta tokens) when routing scan type. Labs
%                       with a different naming convention turn this off.
%     defaultScanType  (string)    — scan-type route used when nothing else
%                       disambiguates: "twoThetaOmega" | "omega" | "phi" |
%                       "xrr" | "rsm".
%
%   The struct is intentionally flat and JSON-friendly (no nested structs,
%   no function handles) so jsonencode/jsondecode round-trips it cleanly.

    cfg = struct( ...
        'wavelength',       1.5406, ...
        'applyWavelength',  true, ...
        'substrate',        "SrTiO3", ...
        'rcShape',          "gauss", ...
        'useFilenameRules', true, ...
        'defaultScanType',  "twoThetaOmega");
end
