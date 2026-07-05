# Changelog

## v1.1.0 — 2026-07-05

First tagged release. Ships a runnable Windows installer for users without
MATLAB, plus the MATLAB source. Notable capabilities since the original hand-off:

- **XRR slab-model fitting** — Parratt reflectivity with Névot–Croce roughness;
  seeded multi-start fit of thickness / density / roughness with parameter
  standard errors.
- **RSM strain & composition** — `analyzeStrainRSM`: auto-finds RSM peaks,
  derives in-plane / out-of-plane lattice parameters, biaxial strain, the relaxed
  lattice parameter, and PZT composition. Geometry validated on real TiO₂ data;
  strain/composition validated on synthetic data this iteration.
- **Material identification** — `identifyMaterial` from lattice parameters.
  A known PTO-vs-PZT (00l) ranking caveat remains and is documented.
- **Automatic peak prominence** — data-driven minimum-prominence selection for
  peak detection (conservative default).
- **Expanded I/O and materials** — reads Rigaku SmartLab Studio II `.hgx`;
  added LaAlO₃, TiO₂, and WO₃ with formula + bulk density.

Requires (installer path) Windows + internet once for the MATLAB Runtime R2026a;
(source path) MATLAB R2022b+.
