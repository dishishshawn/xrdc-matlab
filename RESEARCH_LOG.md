# Research log — xrdc-matlab

Summer research project (Paik group). Goal: **log 160 research hours**; time spent
developing this MATLAB port of XRDC counts toward that total. This file is the running
record of what was done each session so the hours are defensible and traceable to commits.

**How to use:** add one entry per work session. Put your real hours in the Hours column,
keep the summary short, and list concrete changes. Where possible, name the commit(s) so an
entry maps to a diff. Update the running total at the top.

| Metric | Value |
| ------ | ----- |
| Target hours | 160 |
| Logged so far | _fill in_ |
| Remaining | _fill in_ |

---

## 2026-06-05 — Toolboxes, standalone-exe plan, validation, feature audit
**Hours:** _fill in_

**Summary:** Established what's needed to install/run the toolkit, scaffolded a standalone
executable build path (planned, not built), set up a validation workflow and cross-checked
the app against another researcher's (Tushar) data, and produced an accurate feature checklist.

**Changes:**
- Documented required/recommended MATLAB toolboxes and how to install them in bulk
  (README install section already covered this; confirmed the 3 must-haves: Signal
  Processing, Optimization, Curve Fitting).
- Added standalone-executable build path: `xrdc-matlab/build/buildStandalone.m`
  (compiles `xrdcApp` + `+xrdc` to `XRDC.exe`), `build/RELEASE_CHECKLIST.md` (compile is
  the final gate), and a new README §7 documenting build + MATLAB Runtime requirements.
  Verified the build script parses in MATLAB; confirmed MATLAB Compiler is licensed but
  **not yet installed** (`ver("compiler")` empty) — install before building.
- Set up `validation/tushar/` (input / reference / output, git-ignored payloads) with a
  README + RESULTS template for cross-checking app output against researcher data.
- **Validated sample S25 (SRO on STO(100)) against Tushar's graphs:**
  - θ-2θ: app peak positions (STO 100/200/300 = 22.766 / 46.493 / 72.595°) match theory
    and his figure to <0.03°.
  - Rocking-curve FWHM matches his numbers to ~1% **with a Gaussian fit** (sub 0.0304° vs
    0.0308°, film 0.0911° vs 0.0920°); the app's default Lorentzian reads ~16% low. Logged
    full fit-shape sweep in `validation/tushar/RESULTS.md`.
  - XRR: could not validate — he sent the XRR figure but not the XRR `.txt` (follow-up).
- Added `xrdc-matlab/docs/FEATURES.md` — accurate feature checklist reconciled against the
  `+xrdc/` source (done / gated / wanted / stubbed). Confirmed only `readRigakuRaw` and
  `readRigakuRas` are stubs.
- Housekeeping: freed local disk space for the toolbox install; noted that `CLAUDE.md`
  still references the deleted `PROJECT_PLAN.md` / `ALGORITHM_SPEC.md`.

**Open follow-ups:**
- Confirm with Tushar that his FWHM used a Gaussian fit; decide whether RC analysis should
  default to/offer Gaussian.
- Get the XRR `.txt` for S25 to validate the Kiessig-thickness path.
- Run the full `runtests` suite on lab MATLAB (stated handoff gate).
- Install MATLAB Compiler when ready to cut the standalone exe.
