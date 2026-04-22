# Conventions for AI contributors to xrdc-matlab

This file gives future Claude / AI sessions enough context to continue the port without re-deriving conventions from scratch. Read this first.

## Context

- Project goal: MATLAB port of Dr. Tassilo Heeg's Delphi `XRDC` tool, for the Paik group.
- Owner on the human side: Shawn (shawnagarwal0@gmail.com). Stakeholder: Dr. Paik.
- Three goals (from Dr. Paik):
    1. Read Rigaku `.raw` natively — the original can't.
    2. Cleaner data handling than the Delphi version.
    3. Consistent publication-ready figure formatting.
- Everything algorithmic is documented in `../ALGORITHM_SPEC.md` with line references to the Delphi source. **Read the spec before changing any numerical routine.**
- Phase-by-phase roadmap lives in `../PROJECT_PLAN.md`. Respect the phasing — don't start RSM work before the I/O and peak-detection phases are solid.

## Code conventions

**Language level.** Assume MATLAB R2022b+. Use `arguments` blocks for every public function (type + size constraints). Use `string` for identifiers/names and `char` only when MATLAB's API forces it (XML DOM, `sprintf` format strings).

**Package layout.** Everything lives under `+xrdc/` in sub-packages (`+io`, `+lattice`, `+signal`, `+peaks`, `+rsm`, `+plot`, `+data`). One function per file. The package path name is always `xrdc.<subpkg>.<fn>`.

**Naming.**
- Functions: `camelCase` verbs (`readXrdml`, `dSpacingFromHKL`).
- Struct fields: `camelCase` (`twoTheta`, not `two_theta`).
- Crystal system strings: lowercase (`'cubic'`, `'tetragonal'`, …).
- Scan type strings (the `scanType` field): `"twoThetaOmega"`, `"omega"`, `"phi"`, `"psi"`, `"theta"`, `"area"`, `"unknown"`.

**Error IDs.** Use `xrdc:<subpkg>:<reason>` for every `error()` call — makes errors greppable and test-able via `verifyError`.

**Units.**
- Angles in public APIs: degrees for 2θ and ω; radians only internally.
- Wavelengths: Å.
- Lattice parameters: Å.
- Film thicknesses: nm (to match how the original reports, even though other units come from Å).

**No unit-printing in code.** Don't put `°` or `Å` inside function names; keep them in docstrings and plot labels.

**Scan struct shape.** Every I/O function returns the shape defined by `xrdc.io.emptyScan()`. Extend that function (not ad-hoc field additions) when a new field is needed everywhere.

## Testing

- One `tests/testXxx.m` file per package. Use `functiontests(localfunctions)`.
- Every new numerical function must have at least:
    - a round-trip or known-answer test,
    - a vectorised/broadcast test if the function takes arrays,
    - an error-path test for expected failure modes.
- Run via `runtests` from the xrdc-matlab root.
- Known-answer tests for lattice math should cite the reference (Substrates.def row, tabulated literature value, etc.).

## Porting rules from the Delphi source

1. **Match the algorithm, not the implementation.** The Delphi code is often Pascal-constrained (no numerical libraries); MATLAB has better primitives. Don't port brute-force grid searches when `lsqcurvefit` or `findpeaks` is the obvious modern answer — but leave a `'legacy'` flag option on peak detection / fitting for reproducing old analyses bit-for-bit (§1.2 of PROJECT_PLAN).
2. **Preserve known asymmetries.** The RSM transform has a deliberate θ-asymmetry documented in ALGORITHM_SPEC §7.1. Don't "fix" that.
3. **German-locale tolerance on input only.** We read files that use `,` as decimal separator; we never write them.
4. **Drop the auto-updater, the Picker format, and gnuplot.**

## What each cowork session should do first

1. Re-read this file and the algorithm spec.
2. Run `runtests` to confirm the suite is green on the current state.
3. Check `../PROJECT_PLAN.md` to see which phase is in flight.
4. Check outstanding tasks with `TaskList` if using the task tools.
5. If starting Rigaku I/O: confirm sample files are available in `../test-data/` (or equivalent) — if not, the phase is still blocked.

## Don't

- Don't commit tests that require real Rigaku files to pass — gate them on file existence via `assumeTrue(isfile(...))`.
- Don't introduce dependencies on unlicensed toolboxes. Each toolbox dependency should be listed in README "Requirements" and have a fallback path.
- Don't reformat the Delphi source — it's reference material, not a target.
