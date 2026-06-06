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
      (As of this writing it is licensed but **not yet installed** — install via the
      MATLAB installer / Add-Ons before building.)
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

---

**Housekeeping note:** `CLAUDE.md` still points to `../PROJECT_PLAN.md` and
`../ALGORITHM_SPEC.md`, which were removed in commit `453dfd5`. Those references are
stale — worth cleaning up or restoring separately from this work.
