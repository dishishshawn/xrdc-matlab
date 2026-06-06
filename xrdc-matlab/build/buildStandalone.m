function out = buildStandalone(opts)
%BUILDSTANDALONE  Compile the XRDC GUI into a standalone Windows executable.
%
%   buildStandalone() compiles xrdcApp.m and the +xrdc package into a
%   double-clickable XRDC.exe that runs WITHOUT a MATLAB license — end
%   users only need the (free, royalty-free) MATLAB Runtime.
%
%   Usage (from the repo root):
%     >> addpath build
%     >> buildStandalone                 % build into build/standalone/
%     >> buildStandalone(Embed=true)     % also bundle the Runtime installer
%
%   Requirements (build machine only):
%     - MATLAB Compiler license (check with `ver` — "MATLAB Compiler").
%     - The toolboxes whose *better* code paths you want baked in
%       (Signal Processing, Optimization, Curve Fitting). The toolkit has
%       pure-MATLAB fallbacks, but the compiler bundles whatever path THIS
%       machine resolves at build time — so build on a fully-licensed box.
%
%   Name-Value options
%     OutputDir (string)  Destination folder. Default: <repo>/build/standalone.
%     AppName   (string)  Executable base name. Default: "XRDC".
%     Embed     (logical) If true, package the MATLAB Runtime installer into
%                         the output so peers don't download it separately
%                         (much larger artifact). Default: false (web installer).
%     Verbose   (logical) Print progress. Default: true.
%
%   Output
%     out : the compiler.build.Results object (paths, included files, log).
%
%   Notes
%     - The +xrdc/+data/*.json files are added explicitly: the compiler's
%       dependency analyzer cannot see files read via dynamically-built
%       `fileread` paths, so they would otherwise be omitted.
%     - Produces a Windows app (no console window). For macOS/Linux, swap
%       standaloneWindowsApplication -> standaloneApplication below and run
%       the build on that platform (the Runtime is platform-specific).

    arguments
        opts.OutputDir (1,1) string = ""
        opts.AppName   (1,1) string = "XRDC"
        opts.Embed     (1,1) logical = false
        opts.Verbose   (1,1) logical = true
    end

    % --- Resolve repo root relative to this file (build/ is under root) ---
    thisDir  = fileparts(mfilename("fullpath"));
    repoRoot = fileparts(thisDir);

    appFile = fullfile(repoRoot, "xrdcApp.m");
    pkgDir  = fullfile(repoRoot, "+xrdc");
    dataDir = fullfile(pkgDir, "+data");

    if ~isfile(appFile)
        error("xrdc:build:appNotFound", ...
            "Could not find xrdcApp.m at %s — run this from the repo root.", appFile);
    end

    % --- Preflight: confirm MATLAB Compiler is available ---
    if isempty(ver("compiler"))
        error("xrdc:build:noCompiler", ...
            "MATLAB Compiler is not installed/licensed on this machine. " + ...
            "Check `ver` — you need the ""MATLAB Compiler"" product to build.");
    end

    if opts.OutputDir == ""
        opts.OutputDir = fullfile(repoRoot, "build", "standalone");
    end
    if ~isfolder(opts.OutputDir)
        mkdir(opts.OutputDir);
    end

    % --- Data files the dependency analyzer can't discover on its own ---
    extraFiles = strings(0, 1);
    for f = ["substrates.json", "xrayLines.json"]
        p = fullfile(dataDir, f);
        if isfile(p)
            extraFiles(end+1, 1) = p; %#ok<AGROW>
        end
    end

    if opts.Verbose
        fprintf("Compiling %s -> %s\n", appFile, opts.OutputDir);
        fprintf("  bundling %d data file(s); +xrdc resolved by dependency analysis\n", ...
            numel(extraFiles));
    end

    if opts.Embed
        runtimeOpt = "installer";   % web-based MCR installer bundled
    else
        runtimeOpt = "none";        % peers fetch the Runtime themselves
    end

    % --- Build (Windows GUI app: no console window) ---
    out = compiler.build.standaloneWindowsApplication(appFile, ...
        "ExecutableName",   opts.AppName, ...
        "AdditionalFiles",  cellstr([pkgDir; extraFiles]), ...
        "OutputDir",        char(opts.OutputDir), ...
        "EmbedArchive",     true, ...
        "RuntimeIncludeOption", runtimeOpt, ...
        "Verbose",          matlab.lang.OnOffSwitchState(opts.Verbose));

    if opts.Verbose
        fprintf("\nDone. Executable written to:\n  %s\n", ...
            fullfile(char(opts.OutputDir), opts.AppName + ".exe"));
        fprintf("Peers run it after installing the matching MATLAB Runtime " + ...
            "(see deploy notes in README).\n");
    end
end
