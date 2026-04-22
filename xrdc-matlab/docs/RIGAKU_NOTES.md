# Rigaku SmartLab format notes

What we know about the files the Paik lab's Rigaku SmartLab produces, and how `xrdc.io.readRigakuTxt` handles them. Written 2026-04-22 from the first drop of real files in `rexdrctomatlabport_rigakudatasets/`.

---

## File types observed

| Extension | What it is                                     | In scope? |
| --------- | ---------------------------------------------- | --------- |
| `.txt`    | SmartLab ASCII export — 2-column intensity vs angle | **Yes** — primary ingest path |
| `.hgx`    | GenX reflectivity-fit project file (HDF5)      | No — third-party fit software save state; contains simulated curves, fit parameters, optimiser logs. Not raw data. |
| `.raw` / `.ras` | Rigaku binary / RAS ASCII as documented in ALGORITHM_SPEC §2.5 | Not observed in this lab. Parsers remain stubbed with a clear `xrdc:io:notImplemented` error. |

Only `.txt` is actually read by the pipeline. `.hgx` files are ignored.

---

## `.txt` format anatomy

Two variants observed in a single data drop.

### Variant A — headered (most files)

```
<UTF-8 BOM: EF BB BF><SampleName>\t\n
2θ, °\tIntensity, cps\n
<angle>\t<intensity>\n
<angle>\t<intensity>\n
...
```

Characteristics:
- UTF-8 BOM at byte 0 (`EF BB BF`). `readRigakuTxt` strips it on both `char(65279)` and byte-sequence forms.
- Line 1: sample name, tab-terminated, CRLF or LF line ending.
- Line 2: column labels. The `, °` in `"2θ, °"` is a **label decoration**, not a decimal separator — do not interpret as German locale.
- Column labels say `2θ` even when the data are ω (rocking curves) or φ (phi scans). This is a Rigaku quirk; the actual scan axis lives in the filename, not the header.
- Data: tab-delimited, dot decimal, integer or float counts.

### Variant B — headerless (some XRR files)

```
<angle>\t<intensity>\n
<angle>\t<intensity>\n
...
```

No BOM, no sample name, no column labels — data starts at line 1. Observed on several `_XRR_` files in the S04, S05, S06, S07, S10, S11 series (not all XRRs; the headered form appears too).

Detection: `readRigakuTxt` tries to parse line 1 as two floats. Success → headerless (line 1 is data). Failure → treat first 2 lines as header, start data at line 3.

---

## Scan-type inference

Rigaku's ASCII export always labels the x-axis `"2θ, °"` regardless of the physical scan axis. The filename carries the real information.

`readRigakuTxt` maps the filename stem (case-insensitive) to `scanType`:

| Filename substring                    | `scanType`        | Physical axis |
| ------------------------------------- | ----------------- | ------------- |
| `RC`, `rocking`                       | `"omega"`         | ω detuning   |
| `phi`                                 | `"phi"`           | φ            |
| `psi`                                 | `"psi"`           | ψ            |
| everything else (`2theta omega`, `XRR`, `th2th`, ...) | `"twoThetaOmega"` | coupled 2θ / ω |

If your filename convention differs, override by assigning `scan.scanType` after `readScan`.

---

## Dispatcher routing

`xrdc.io.readScan` identifies Rigaku `.txt` via **two cues**:

1. **Header marker**: the string `"INTENSITY, CPS"` (case-insensitive) anywhere in the first 512 bytes. Unique to Rigaku exports.
2. **Filename prefix**: `TR_` (the Paik lab's Rigaku naming convention). Catches the headerless variant.

Either cue routes to `readRigakuTxt`. Otherwise the file falls through to the plain-text parser — which will parse the data but not set the Rigaku scan-type mappings.

---

## Wavelength

The SmartLab export has no in-file wavelength record. `readRigakuTxt` sets `scan.lambda = 1.5406` (Cu Kα₁), the instrument's standard source. If you're on a non-Cu anode or want higher precision, override:

```matlab
scan.lambda = 1.5405929;   % Cu Kα₁ NIST
```

---

## Step uniformity

The Delphi XRDC (xrdc1.pas:1235) warns when the 2θ step is non-uniform. `readRigakuTxt` does the same (tolerance 1e-3°). Non-uniform data aren't rejected — just noted — because Rigaku sometimes dwell-corrects detector counts in ways that shift the nominal angle by µ-degrees without breaking the analysis.

---

## What's NOT parsed

The SmartLab export is lossy relative to the instrument's internal state. Not available in `.txt`:

- Goniometer zero offset
- Tube voltage / current / monochromator
- Scan mode (step vs continuous)
- Per-point counting time
- Omega offset for RSM slices (each slice is a separate file; we infer omega from the secondAxis set by the I/O layer, which is currently `NaN` for `.txt` — user must set it manually for RSM work with Rigaku data)

For RSM with Rigaku `.txt` slices, the workflow is:

```matlab
scans = xrdc.rsm.loadAreaScan(folder, 'Pattern', 'TR_*_omega-*.txt');
for i = 1:numel(scans)
    scans(i).secondAxis = parseOmegaFromFilename(scans(i).identifier);
end
xrdc.plot.plotRsm(scans);
```

If this becomes common, add an `OmegaFromFilename` option to `loadAreaScan`.

---

## `.hgx` (GenX) — why it's ignored

Magic bytes `89 48 44 46 0D 0A 1A 0A` — HDF5. Probe structure:

```
current/                          NXentry
  config                          object (pickled config)
  data/datasets/0/
    x                             float64[N]  — 2θ
    y                             float64[N]  — observed reflectivity
    y_raw                         float64[N]  — pre-reduction counts
    y_sim                         float64[N]  — simulated curve
    y_fom                         float64[N]  — figure-of-merit
    error                         float64[N]  — reduced error bars
  optimizer/
    fom_log                       float64[K,2]
    solver, solver_module         object
  parameters/
    data col 0..5                 fit parameter table
  script                          object (Python reduction script)
```

This is a full GenX project save — raw data, fit model, optimiser state, the lot. The raw data is available in `y_raw` but the matching `.txt` file sits next to it, so there's no reason to parse the HDF5.

If we ever need `.hgx` support: MATLAB's `h5read` handles HDF5 natively. Read `current/data/datasets/0/x` and `y_raw`, treat as Variant B (headerless text-equivalent).

---

## Adding support for a new lab's Rigaku files

Steps if another group drops a different Rigaku variant:

1. `head -5` on 3–5 sample files. Note BOM, header lines, column labels, decimal separator.
2. Add a new cue to `isRigakuTxt` in `readScan.m` if the existing "INTENSITY, CPS" / `TR_` heuristics miss.
3. If the layout diverges substantially from both variants A and B, branch inside `readRigakuTxt` on a detected signature.
4. Add a `testReadRigakuTxt<LabName>` integration test gated on file existence.

Do not add a new reader function unless the format is fundamentally different. `.ras` (key=value + `*RAS_INT_START`) and binary `.raw` (magic bytes `RAW1.0x`) are separate parsers because the framing differs; `.txt` variants from different instruments should stay in one reader.
