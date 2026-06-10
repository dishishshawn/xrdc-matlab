# SP79 superlattice — `xrdc.lattice.superlatticePeriod` validation

Sample: **SP79 WO₃/LAO superlattice on LAO(100)**, nominally 3 uc × 26 stacks,
~60 nm total (700 °C, 50 mT). Cu Kα θ–2θ scan, 12–80°.
Data: `data/tushar superlattice data/SP79_..._2theta_omega_04232026.txt`.

Run: `matlab -batch "run('validation/tushar_superlattice/runSuperlattice.m')"`

## Results

| Test | Λ (fit) | R² | n peaks | Notes |
|------|---------|----|---------|-------|
| Synthetic, Λ=15 nm (known answer) | **15.000 ± 0.000 nm** | 1.000000 | 7 | exact recovery, error 0.00% |
| Real, LAO(001) film tail (21.0–23.45°) | **59.6 ± 3.8 nm** | 0.992 | 6 | fringes at 22.65–23.41° |
| Real, LAO(002) film tail (45.6–47.9°) | **58.2 ± 34.8 nm** | 0.992 | 3 | fewer fringes resolved |

## Interpretation

1. **The routine is correct.** The synthetic known-answer case recovers the
   15 nm period exactly (R² = 1).

2. **Both orders independently give ≈ 59 nm**, matching the **nominal 60 nm
   total film thickness** — not the ~2.3 nm chemical bilayer period
   (60 nm / 26 stacks). The series `findPeaks` locked onto is the
   **total-thickness Laue/Kiessig fringe train** (spacing ≈ 0.13–0.16° in 2θ),
   which is what is physically resolved in this scan.

3. **No resolved superlattice (bilayer) satellites exist in this scan.** A
   2.3 nm period would place 1st-order satellites ~±4° from each Bragg peak
   (≈ 19.5° / 27.3° around the 001); the wide scan shows only background there.
   The WO₃/LAO scattering contrast is too low to produce visible satellites.

This is exactly the separation the routine's own header documents: it returns a
*period* from evenly-(sin θ)-spaced peaks and cannot, by itself, distinguish a
chemical bilayer repeat from a total-thickness fringe train — that the recovered
value equals the total thickness is a property of *this* sample's data, not a
routine error. Cross-check total thickness against XRR
(`xrdc.xrr.analyzeFringes`) as the header advises.
