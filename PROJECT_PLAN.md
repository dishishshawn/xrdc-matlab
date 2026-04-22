---
title: XRDC → MATLAB Port — Project Plan
companion: ALGORITHM_SPEC.md
status: draft — awaiting Dr. Paik's review of §1 decision points
---

# XRDC → MATLAB Port — Project Plan

This plan turns the algorithm spec (`ALGORITHM_SPEC.md`) into an ordered delivery schedule. Each phase is a standalone deliverable Dr. Paik can use in the lab; later phases add capability without breaking the API of earlier ones.

Design principles (carried over from the spec):

- **Reproduce where it matters, modernise where it helps.** Match the original's numerical results on canonical test sets, but prefer `findpeaks` / `lsqcurvefit` / `fitlm` / `Savitzky–Golay` over the hand-rolled brute-force code where that gives better uncertainties or robustness for free.
- **Scripting API first, GUI later.** A clean `+xrdc` package that works from the command line is the primary deliverable. A GUI wrapper is a nice-to-have and a late-phase decision.
- **Every algorithm has a test.** Phase "done" means the unit tests pass *and* a representative real dataset from the lab has been processed end-to-end.

---

## 1. Decision points for Dr. Paik

These are the questions that shape the scope of the port. Getting answers before Phase 0 ends avoids rework later.

1. **MATLAB version and toolboxes.** The plan assumes R2022b+ with *Signal Processing*, *Optimization*, *Curve Fitting*, and *Statistics and Machine Learning* toolboxes. If any of these are missing from the lab's license, some algorithms need alternative implementations (e.g. no Curve Fitting Toolbox → keep the brute-force peak fit as the default).

2. **Parity vs. improvement.** When the original algorithm is merely adequate, do we preserve it exactly or upgrade? Recommendation: upgrade by default, but keep a `'legacy'` option on peak detection and peak fitting so old analyses can be reproduced bit-for-bit if ever needed.

3. **Rigaku `.raw` variants.** Which files do you have? The port must handle at least:
   - the lab's current instrument output (almost certainly `RAS` ASCII or `RAW1.01`/`RAW1.02` binary)
   - and any historical `.raw` files that should still open.
   
   Sharing 3–5 sample files covering all variants you care about is the fastest way to scope this.

4. **Legacy formats.** PANalytical `.xrdml` and `.x00` are required (Dr. Heeg's original datasets). Picker `.596`/`.1035` is obsolete — drop?

5. **GUI scope.** Options, in increasing effort:
   - **(a) None** — scripting API only. Users call `xrdc.io.readScan` and friends from scripts / the command window.
   - **(b) Light App Designer GUI** — one window for loading a scan, plotting, peak-finding, and RSM. ~1 week.
   - **(c) Full port of the 17-form Delphi UI.** ~3–4 weeks. Not recommended; the Delphi UI is dated and the scripting API is more flexible for publication workflows.

6. **Publication-plot style guide.** The spec's §10 proposes defaults (Arial 18 pt ticks, Helvetica axis labels, etc.). Does Dr. Paik have a preferred style? If there is an existing paper whose figures we should match, the fastest path is to match it exactly.

7. **Deliverable form.** Git repo the lab maintains? A zipped package handed off once? Integrated with Dr. Paik's group's existing MATLAB code?

---

## 2. Phase roadmap

Rough sizing: 8–10 weeks of focused part-time work for a single developer, front-loaded on the unfamiliar bits (Rigaku binary parsing, RSM). Each phase below is a checkpoint — stop, share with Dr. Paik, iterate.

### Phase 0 — Project setup (3–5 days)

**Goal:** a working skeleton you can run `runtests` against and a shared understanding of scope.

- [ ] Create git repo layout per §12 of `ALGORITHM_SPEC.md`:
  ```
  xrdc-matlab/
    +xrdc/ …                      (code as per the spec's §12)
    tests/                        (matlab.unittest)
    examples/                     (small scripts showing each feature)
    data/test/                    (~10 MB of reference files)
    docs/                         (README.md + ALGORITHM_SPEC.md + this plan)
    CLAUDE.md                     (conventions: naming, error strategy, etc.)
  ```
- [ ] `runtests` CI via GitHub Actions on a MATLAB Online runner (or documented manual run if the lab doesn't use GitHub).
- [ ] Collect test data: 3× PANalytical `.xrdml`, 1× `.x00`, 3–5× Rigaku `.raw` (covers all variants), and the Delphi-program-produced output for each (peak list, lattice params, RSM image) for parity checks.
- [ ] Pick and commit a code-style (MathWorks defaults + 4-space indent + camelCase public + `x_private` or `pPrivate`).
- [ ] Get answers to §1 decision points from Dr. Paik.

**Exit criteria:** empty but passing test suite; at least one sample file per format confirmed to load in the Delphi program for comparison.

### Phase 1 — I/O and data model (5–8 days)

**Goal:** every file the lab owns can be read into a uniform struct.

- [ ] `+xrdc/+io/readScan.m` — dispatcher that sniffs format.
- [ ] `readXrdml.m` — PANalytical XML. Unit test compares element-for-element against a hand-decoded reference.
- [ ] `readTextScan.m` — two-column text.
- [ ] `readPhilipsX00.m` — `HR-XRDSCAN` header + SCANDATA block.
- [ ] **`readRigakuRas.m` — ASCII `.ras` parser.** New territory. Parse `*MEAS_COND_*` keys, `*RAS_INT_START` data block.
- [ ] **`readRigakuRaw.m` — binary `.raw` parser** (RAW1.01 / RAW1.02). New territory. Reference: Rigaku SmartLab application notes + any existing community implementations (Python's `xylib`, `rigaku-parser` crate, etc. — for cross-reference only; the port should not take dependencies).
- [ ] Both Rigaku parsers must correctly detect scan axis ("2θ/θ", "ω", etc.) and populate `secondAxis` for θ-2θ scans.
- [ ] Unit tests for every parser: `twoTheta`, `counts`, `scanType`, `secondAxis` all match expected values.

**Exit criteria:** at least one file of each format loads and plots (raw, no processing) on screen with axes labelled.

**Risk:** the Rigaku binary `.raw` layout is proprietary and under-documented. If a sample file breaks the parser, we may need to ask Rigaku or reverse-engineer from a hex dump. Budget 2 extra days for this.

### Phase 2 — Signal processing + single-scan plot (3–5 days)

**Goal:** load a scan, smooth, subtract background, show it the way a paper figure would look.

- [ ] `+xrdc/+signal/smooth.m` — moving-average + Savitzky-Golay options.
- [ ] `+xrdc/+signal/subtractBackground.m` — rolling mean (legacy) + `msbackadj`-style modern option.
- [ ] `+xrdc/+signal/slopes.m` — Savitzky-Golay first/second derivative (replaces `CalcSlopes`).
- [ ] `+xrdc/+plot/plotScan.m` with the invariants from spec §10 (Arial, log-Y, deterministic colours, `exportgraphics` for PDF/PNG).
- [ ] `+xrdc/+plot/publicationStyle.m` — applies the style to any existing `axes` handle, so users can adopt the style on their own plots.

**Exit criteria:** a script `examples/plotSingleScan.m` loads a real Rigaku file and produces a PDF figure ready for a lab notebook.

### Phase 3 — Peak detection + fitting (5–8 days)

**Goal:** find peaks and fit their positions/widths with reportable uncertainties.

- [ ] `+xrdc/+peaks/findPeaks.m` — built-in `findpeaks` with prominence/width options, sensible defaults for XRD.
- [ ] `+xrdc/+peaks/findPeaksLegacy.m` — ports `ScanPeak1` and `ScanPeak2` behind a `'algorithm','slope'|'slope2'` option. Only needed if decision point §1.2 lands on "parity".
- [ ] `+xrdc/+peaks/adjustPeaks.m` — FWHM-based refinement with `interp1(...,'pchip')`.
- [ ] `+xrdc/+peaks/fitPeak.m` — `lsqcurvefit` with Lorentz/Gauss/Voigt (bonus) shapes; returns parameters + Jacobian-derived uncertainties.
- [ ] Parity test: on a synthetic Lorentzian, both the legacy brute-force and the modern `lsqcurvefit` paths should recover the known FWHM to within 0.001°.
- [ ] Parity test: on a real θ-2θ scan, peak positions from `findPeaks` should match the Delphi program's output within 0.01° for well-isolated peaks.

**Exit criteria:** an example script fits a single Bragg peak of a real substrate (e.g. SrTiO₃(002)) and reports position, FWHM, intensity with uncertainties.

### Phase 4 — Lattice and Nelson–Riley (3–5 days)

**Goal:** extract lattice parameters from a peak list with proper uncertainties.

- [ ] `+xrdc/+lattice/bragg.m` — `d ↔ 2θ ↔ λ` with energy-to-wavelength conversion (ships with `xrayLines.json`).
- [ ] `+xrdc/+lattice/dSpacingFromHKL.m` — all seven crystal systems from spec §8.
- [ ] `+xrdc/+lattice/nelsonRiley.m` — `fitlm` on `a_i` vs. `cos²θ/sinθ + cos²θ/θ_deg`. Returns intercept, slope, and 95 % CI.
- [ ] Port `Substrates.def` → `+xrdc/+data/substrates.json` (JSON because MATLAB parses it natively and it's diff-friendly).
- [ ] Port `XRAY.def` → `+xrdc/+data/xrayLines.json`.

**Exit criteria:** run on a real Si(400)/Si(331)/Si(422)/Si(511)/Si(440) multi-peak set; intercept lattice parameter agrees with Delphi output to 1e-4 Å.

### Phase 5 — Structure simulation and film thickness (2–4 days)

**Goal:** overlay predicted Bragg peaks on a real scan; estimate film thickness from fringes.

- [ ] `+xrdc/+lattice/simulatePattern.m` — vectorised hkl expansion + crystal-system dispatch. Returns `(h,k,l,d,2θ)` table.
- [ ] Duplicate-elimination rule from spec §8 (merge within 1e-6 Å, tie-break by `min(h+k+l)` then `min(l·10000+k·100+h)`).
- [ ] `+xrdc/+lattice/thicknessFromFringes.m` — `polyfit(sin(θ), 1:N)` form; supports half-integer fringe-index mode.
- [ ] Overlay predicted peaks on `plotScan` with a flag.

**Exit criteria:** script that overlays SrTiO₃ predicted peaks on a real SrTiO₃ substrate scan; and a script that calculates thickness from a Kiessig-fringe set and reports thickness ± uncertainty.

### Phase 6 — Reciprocal-space maps (5–8 days)

**Goal:** convert an area scan into a publication-quality RSM, with interactive offset-setting.

- [ ] `+xrdc/+rsm/toReciprocalSpace.m` — the transform in spec §7.1, preserving the θ asymmetry called out there.
- [ ] `+xrdc/+rsm/loadAreaScan.m` — ingest a folder/pattern of θ-2θ slices with varying ω.
- [ ] `+xrdc/+plot/plotRsm.m` — `imagesc` or `contourf` on `(k_par, k_perp)`; colour map matches spec §10 (turbo or a custom-defined `rcol/gcol/bcol`).
- [ ] `+xrdc/+rsm/setOffsetsInteractive.m` — figure with a `ginput`-style callback: user clicks the expected substrate peak, enters the known (2θ, ω), and the function returns `(ΔΘ, ΔΩ)`.
- [ ] Parity test: a published RSM dataset (a previously-analysed sample from Dr. Paik) produces a map visually indistinguishable from the reference.

**Exit criteria:** one-command RSM from a folder of slices; interactive offset-setting works; PNG export at 600 dpi.

### Phase 7 — App wrapper (optional, 7–14 days)

**Goal:** a GUI for non-scripting users, if Dr. Paik decides it is worth the effort (§1.5).

- [ ] App Designer single-window app: file load, scan preview, peak finder, fit panel, RSM viewer, export.
- [ ] Wraps — does not reimplement — the `+xrdc` scripting API. The scripting API remains the system of record.
- [ ] File history (mimics the `xrdc.INI` `[FILEHISTORY]` behaviour).
- [ ] No auto-updater (§11 of the spec; dropped).

**Exit criteria:** someone who has never used MATLAB can load a scan, find peaks, and export a figure without writing code.

### Phase 8 — Validation and docs (4–7 days)

**Goal:** ship-ready.

- [ ] Run the §13 verification plan end-to-end on ≥ 3 datasets.
- [ ] `docs/USER_GUIDE.md` — how to install, how to load a scan, how to do each workflow.
- [ ] `docs/RIGAKU_NOTES.md` — what we learned about the binary format, for future maintenance.
- [ ] `examples/` — 5–7 scripts that together cover every public API.
- [ ] Handoff meeting with Dr. Paik walking through the above.

**Exit criteria:** Dr. Paik (or a student) can go from a fresh clone to a publication figure on one of the lab's real datasets, using only the docs.

---

## 3. Dependencies and risks

| Risk | Phase | Mitigation |
|---|---|---|
| Rigaku binary `.raw` layout is undocumented | 1 | Have sample files + reference Python `xylib` layout; budget 2 extra days |
| No Curve Fitting Toolbox in lab license | 3 | Fall back to `lsqcurvefit` (Optimization Toolbox, more common) |
| No Optimization Toolbox either | 3 | Keep the legacy brute-force as the default; document the accuracy trade-off |
| PANalytical XML namespace quirks | 1 | Use local-name matching instead of fully-qualified tags |
| RSM parity with Delphi binary on pre-existing samples | 6 | The θ asymmetry in the transform is a known trap (spec §7.1); dedicated parity test on a known-good dataset |
| Publication plot style disagreement | 2 | Expose the style as a name–value arg; ship three presets (`'paper'`, `'talk'`, `'notebook'`) |

---

## 4. What I'd build first if I only had a week

In case Dr. Paik wants a quick demo before committing to the whole plan:

1. Day 1: `readRigakuRas` or `readRigakuRaw` (whichever matches the lab's files) + `readXrdml`.
2. Day 2: scan struct + `plotScan` with the publication style.
3. Day 3: `findPeaks` + `fitPeak` (Lorentz via `lsqcurvefit`).
4. Day 4: `bragg` + `dSpacingFromHKL` + `simulatePattern` overlay on the plot.
5. Day 5: a demo script that takes a `.raw` file and produces a PDF figure with peak positions and simulated substrate reflections overlaid.

That's ~70 % of the day-to-day workflow a grad student needs, and it establishes the package shape for the rest of the work.

---

## 5. Open questions for the next conversation with Dr. Paik

- Can we see a few sample `.raw` files? (Blocker for Phase 1.)
- Which MATLAB toolboxes does the lab have?
- Is there a reference paper whose figure style we should match?
- Is a GUI a hard requirement, or is a scripting API acceptable?
- Who will maintain this code after the initial port — Shawn? A student? Dr. Paik himself?
- What is the realistic timeline? Is "working demo in 1 week, production-quality in 8" about right?

---

*The plan lives in the repo and is expected to change. Update it as answers to §1/§5 come in and as phases complete.*
