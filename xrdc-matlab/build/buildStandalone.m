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
%     SingleFile (logical) If true (default), embed the CTF archive INSIDE the
%                         exe for a single-file deliverable. If false, the code
%                         ships as a companion XRDC.ctf that must travel next to
%                         the exe. Set false if antivirus blocks the build (see
%                         Notes). Default: true.
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
%     - ANTIVIRUS RACE: with SingleFile=true the builder re-opens the freshly
%       written exe to embed the archive; Windows Defender can grab a handle on
%       the new unsigned exe in that window, causing a "process cannot access
%       the file ... being used by another process" failure (often leaving a
%       small partial exe). Fixes, in order of preference: (1) add a Defender
%       exclusion for this build folder, or (2) build with SingleFile=false
%       (single write, no re-open) and distribute exe + XRDC.ctf together.

    arguments
        opts.OutputDir  (1,1) string  = ""
        opts.AppName    (1,1) string  = "XRDC"
        opts.Embed      (1,1) logical = false
        opts.SingleFile (1,1) logical = true
        opts.Verbose    (1,1) logical = true
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

    % --- Build (Windows GUI app: no console window) ---
    try
        out = compiler.build.standaloneWindowsApplication(appFile, ...
            "ExecutableName",   opts.AppName, ...
            "AdditionalFiles",  cellstr([pkgDir; extraFiles]), ...
            "OutputDir",        char(opts.OutputDir), ...
            "EmbedArchive",     matlab.lang.OnOffSwitchState(opts.SingleFile), ...
            "Verbose",          matlab.lang.OnOffSwitchState(opts.Verbose));
    catch err
        if contains(err.message, "being used by another process") && opts.SingleFile
            error("xrdc:build:fileLocked", ...
                "Build failed: the exe was locked mid-embed (usually antivirus " + ...
                "scanning the new file). Either add a Defender exclusion for %s, " + ...
                "or rebuild with buildStandalone(SingleFile=false) to ship exe + " + ...
                "XRDC.ctf together.\nOriginal error: %s", opts.OutputDir, err.message);
        end
        rethrow(err);
    end

    % --- Optionally wrap in an installer that bundles the MATLAB Runtime ---
    % compiler.build.* only produces the raw exe; Runtime packaging is a
    % SEPARATE step via compiler.package.installer (the build call has no
    % runtime option — passing one errors). RuntimeDelivery="installer"
    % embeds the Runtime (large, offline); "web" fetches it at install time.
    if opts.Embed
        if opts.Verbose
            fprintf("Packaging installer with bundled MATLAB Runtime...\n");
        end
        compiler.package.installer(out, ...
            "ApplicationName", char(opts.AppName), ...
            "RuntimeDelivery", "installer", ...
            "OutputDir",       char(fullfile(opts.OutputDir, "installer")));
    end

    if opts.Verbose
        fprintf("\nDone. Executable written to:\n  %s\n", ...
            fullfile(char(opts.OutputDir), opts.AppName + ".exe"));
        fprintf("Peers run it after installing the matching MATLAB Runtime " + ...
            "(see deploy notes in README).\n");
    end
end
