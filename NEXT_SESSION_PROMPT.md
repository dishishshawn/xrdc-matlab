# Prompt for the next cowork session

Paste the block below (between the `---` lines) into a new cowork session when you're ready to continue the port. Adjust the "Specifically, this session" paragraph based on what you want to tackle — suggestions follow the prompt.

---

I'm continuing work on the XRDC → MATLAB port for Dr. Paik's group. The project lives at `XDRC MATLAB Project/`. Before doing anything, please read, in this order:

1. `XDRC MATLAB Project/xrdc-matlab/CLAUDE.md` — conventions you must follow.
2. `XDRC MATLAB Project/ALGORITHM_SPEC.md` — the algorithm reference for every port decision.
3. `XDRC MATLAB Project/PROJECT_PLAN.md` — the phased roadmap and what's blocked on Dr. Paik.
4. `XDRC MATLAB Project/xrdc-matlab/README.md` — current implementation status.

Quick state summary:
- `+xrdc/+lattice/` (Bragg, d-spacing for all 7 crystal systems, Nelson-Riley, Kiessig thickness, structure simulation): **done with tests**.
- `+xrdc/+signal/` (smoothing, background, Savitzky-Golay derivatives): **done with tests**.
- `+xrdc/+io/` (XRDML, Philips `.x00`, plain-text): **done with tests**; Rigaku `.ras`/`.raw` **stubbed** — blocked on sample files from Dr. Paik.
- `+xrdc/+peaks/` (`findPeaks` prominence-based, `findPeaksLegacy` for ScanPeak1/ScanPeak2 ports, `adjustPeaks` FWHM area-bisector refinement, `fitPeak` Lorentz/Gauss/pseudoVoigt with `lsqcurvefit` + Jacobian SEs and a `'Method','bruteforce'` legacy grid option): **done with tests**.
- `+xrdc/+plot/` (`publicationStyle`, `plotScan`, `plotStack`): **scaffolded with tests**. Style defaults come from ALGORITHM_SPEC §10 (Arial 18 pt, log-Y, 3^(j) / 10^(20·j/N) stacking multiplier). Visual parity against a paper figure is pending Dr. Paik's reference.
- Data files `+xrdc/+data/xrayLines.json` + `substrates.json`: **ported from the Delphi `.def` files**.
- Not yet started: RSM (`+rsm`), any GUI.

Verification notes from previous passes (don't redo):
- `nelsonRiley.m` uses the textbook OLS standard-error formulas, **not** the Delphi expression at `xrdc3.pas:281`. That Delphi expression is off by a factor of 1/√n (latent bug). A `fitlm`-gated cross-check test lives in `tests/testLattice.m::testNelsonRileySEvsFitlm`. See ALGORITHM_SPEC §6.3 and the comment block in `nelsonRiley.m` for the reasoning — do not "restore Delphi parity" on this formula.
- Peak detection has two paths on purpose: `findPeaks` (modern prominence-based via MATLAB's `findpeaks`) is the default; `findPeaksLegacy` is a direct port of `ScanPeak1` / `ScanPeak2` from `xrdc1.pas:1482-1564` for bit-for-bit reproduction of old XRDC analyses. Don't merge them.
- `fitPeak` defaults to `lsqcurvefit` with a Jacobian-derived covariance (`sigma2 * pinv(J'*J)`) for parameter SEs. The Delphi 20³ grid search is preserved under `'Method','bruteforce'` — leave it there.
- `+xrdc/+plot/` funnels every style decision through `publicationStyle.m`. When Dr. Paik picks a reference paper, change defaults there only; `plotScan` / `plotStack` shouldn't need edits.
- Delphi's `plotStack` uses a random RGB per trace. We intentionally replaced this with a deterministic palette (ALGORITHM_SPEC §10 rationale) — don't restore the random behaviour.
- `plotScan` clamps counts ≤ 0 → 1 when `LogY=true` to match the Delphi `PaintGraph` clamp at `xrdc1.pas:1901-1904`. `LogY=false` leaves them alone.

Specifically, this session: **<FILL IN — see suggestions below>**

Ground rules:
- Follow the conventions in `CLAUDE.md` (arguments blocks, camelCase, `xrdc:<subpkg>:<reason>` error IDs, units documented in docstrings).
- Every numerical function needs a unit test. Run `runtests` at the end and confirm green.
- Match algorithms, not implementations. Prefer MATLAB's modern primitives (`findpeaks`, `lsqcurvefit`, `fitlm`) over ports of the Delphi brute-force code, but keep a `'legacy'` option where the spec calls for reproducibility.
- If something's blocked on Dr. Paik's answers (Rigaku files, GUI decision, plot style), stub it and note it — don't guess and build the wrong thing.

Use the task tools to track progress. When done, update `NEXT_SESSION_PROMPT.md` with what was accomplished and what's next, and leave a concise summary at the end of the chat.

---

## Suggestions for what to tackle next (pick one, or ask Shawn)

### A. Rigaku `.raw` / `.ras` parsers (Phase 1 finish — blocked unless files exist) — top priority if files arrived

Check `XDRC MATLAB Project/test-data/` first. If Shawn has dropped 3–5 Rigaku sample files there, implement:

- `+xrdc/+io/readRigakuRas.m` — ASCII format with `*MEAS_COND_*` keys and `*RAS_INT_START` block.
- `+xrdc/+io/readRigakuRaw.m` — binary format (RAW1.01/1.02 first, newer variants only if present in the samples).
- Tests that read the sample files and check `twoTheta`/`counts` shape and reasonable values.
- Update `readScan` dispatcher to call the real parsers instead of erroring.

Reference: ALGORITHM_SPEC §2.5 for the target contract.

### B. Structure-simulation + Kiessig overlay demo (half a day — low effort, high polish)

Now that `+xrdc/+plot/` is in place, build an `examples/demoOverlay.m` script that:
- Loads an XRDML scan via `xrdc.io.readScan`.
- Runs `xrdc.peaks.findPeaks` + `xrdc.peaks.adjustPeaks` to mark the observed peaks.
- Overlays `xrdc.lattice.simulatePattern` predictions for the substrate (pulled from `substrates.json`).
- Overlays Kiessig-fringe predictions from `xrdc.lattice.kiessigThickness`.
- Exports a figure with `exportgraphics` at 600 dpi.

This doubles as an end-to-end regression check that the packages compose cleanly.

### C. Finalize plot style once Dr. Paik picks a reference paper

The defaults in `publicationStyle.m` come from ALGORITHM_SPEC §10 and are deliberately conservative. When Dr. Paik names a reference figure, the changes are:
- Font family, tick/label/title sizes → Name/Value defaults in `publicationStyle.m`.
- Palette → change the default `Palette` arg in `plotStack.m` (or ship a new colormap under `+xrdc/+plot/palettes/`).
- Legend placement, tick density, minor-tick behaviour → all live in `publicationStyle.m`.
- Don't touch `plotScan.m` / `plotStack.m` callers — the whole point of `publicationStyle` is that it's the one knob.

### D. RSM package (`+xrdc/+rsm/`) — Phase 4, larger effort

Reciprocal-space map reader + plotter. Start from `xrdcrsm.pas` in the Delphi source. Likely ~1 week of work and the first real case that needs 2-D data handling. Skip until peaks/plotting have been exercised on real data.

### E. GUI decision

Still open. PROJECT_PLAN §5 calls out App Designer vs. scripts-only. Don't start building a GUI without explicit confirmation from Shawn — the package is script-friendly already.
