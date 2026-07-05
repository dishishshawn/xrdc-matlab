# Standalone executable — release checklist

Compiling `XRDC.exe` is the **last** step before a binary release. The exe freezes
whatever the source does at build time, so everything else should be settled first.
Work top to bottom; don't start the compile until every box above it is checked.

## Why last
- The compiler bakes in the resolved code paths and a fixed version — any source fix
  after the build means a full rebuild + redistribute.
- A binary is only worth cutting once the source is stable and validated.
- It's the one step that needs extra licensed products (MATLAB Compiler) and a large
  Runtime download, so there's no reason to pay that cost early.

## Phase A — source is frozen (do these first, in MATLAB)
- [ ] `runtests` is **green** on the lab MATLAB (the README's stated handoff gate).
- [ ] GUI smoke test: `xrdcApp`, load one file of each type (Rigaku `.txt`,
      PANalytical `.xrdml`, Philips `.x00`), confirm analysis + 600 dpi export.
- [ ] Version/citation strings in `README.md` and the app are current.
- [ ] No pending source edits you'd be unhappy to ship frozen.

## Phase B — build machine prerequisites (one-time)
- [ ] **MATLAB Compiler installed** — `ver("compiler")` returns non-empty.
      (Installed on the lab build machine as of R2026a.)
- [ ] The optional toolboxes installed so the *better* paths get baked in, not the
      pure-MATLAB fallbacks: **Signal Processing, Optimization, Curve Fitting** at
      minimum (Global Optimization + Parallel Computing if you want those features).
- [ ] Note the build MATLAB release (`version`) — end users need the matching Runtime.

## Phase C — build
- [ ] From the repo root:
      ```matlab
      addpath build
      buildStandalone               % → build/standalone/XRDC.exe
      ```
- [ ] For a self-contained artifact (bundles the Runtime installer, much larger):
      ```matlab
      buildStandalone(Embed=true)
      ```
- [ ] Build completes with no unresolved-dependency warnings.

## Phase D — verify the binary (ideally on a machine WITHOUT a MATLAB license)
- [ ] Install the matching **MATLAB Runtime** on a clean machine.
- [ ] Launch `XRDC.exe` — GUI opens, no console/error window.
- [ ] Repeat the Phase A smoke test (load each format, analyze, export) against the
      exe, not MATLAB. Confirm exported figures match the MATLAB output.

## Phase E — distribute
- [ ] Package `XRDC.exe` + a short run note (Runtime version + download link — see
      README §7) for collaborators.
- [ ] Tag the release in git so the binary is traceable to a commit.

## Phase F — source / GitHub release (installer + source)

The repo delivers to two audiences: no-MATLAB users get the web installer, MATLAB
users get the source. Cut both together as one tagged release.

- [ ] Bump the version in **one** place: `+xrdc/version.m`. Confirm it flows to
      the app title (`xrdcApp.m`) and any README references.
- [ ] `runtests` green on lab MATLAB (Phase A).
- [ ] Build the web installer from the frozen source:
      ```matlab
      addpath build
      buildStandalone(Embed=true)   % → build/standalone/installer/XRDCInstaller.exe
      ```
- [ ] **Runtime-verify** `XRDCInstaller.exe` on a clean machine WITHOUT MATLAB
      (Phase D): it installs, fetches the Runtime, launches, and passes the Phase A
      smoke test. The installer bakes in source at build time — rebuild it whenever
      the version changes.
- [ ] Push `main`; create an annotated tag `vX.Y.Z`.
- [ ] Create the GitHub Release (body = the `CHANGELOG.md` section). GitHub
      auto-attaches the source zip.
- [ ] Upload the installer asset:
      `gh release upload vX.Y.Z build/standalone/installer/XRDCInstaller.exe`.
