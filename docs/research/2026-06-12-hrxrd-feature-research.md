# Prioritized Feature Report: High-Value Analysis Additions to an In-House MATLAB HRXRD Toolkit for Epitaxial Perovskite/Oxide Thin Films

> Deep-research output, 2026-06-12. Generated from a research brief about beneficial features to add to xrdc-matlab. Roadmap reference — not yet specced or scheduled. Pairs with `docs/session-logs/2026-06-12.md`.

## TL;DR
- The single highest-value addition is **RSM-based strain/composition decomposition** (Area 1): peak-picking an asymmetric RSM (103/113/204 on STO) into q∥/q⊥, decomposing into in-plane/out-of-plane strain via elastic constants, recovering the relaxed lattice parameter, and only then applying Vegard's law — this closes the fundamental gap that (00l)-only ID cannot resolve (PTO-vs-PZT degeneracy; films compressed below bulk c). It reuses data they already collect.
- The next tier is **full Parratt/Abelès XRR fitting** (Area 2, replacing fragile fringe-spacing) and **quantitative RSM analytics** (Area 5: relaxation line, degree of relaxation, tilt, mosaic-vs-size broadening) — both high-value, medium effort, and built on existing measurements.
- Lower-priority but worthwhile: **uncertainty propagation** (Area 7, low effort, high credibility payoff), **line-profile microstructure analysis** (Area 3, caveated for epitaxial films), **dynamical superlattice simulation** (Area 4, high effort), and **geometry corrections** (Area 6, low effort, accuracy-critical).

## Key Findings
1. **(00l)-only identification is fundamentally under-determined for epitaxial films.** A symmetric scan yields only the out-of-plane spacing d⊥. For a clamped film, the measured c reflects *both* composition and epitaxial strain, which are confounded. The minimal fix is one asymmetric reflection or RSM to obtain the in-plane parameter, enabling strain decomposition and recovery of the strain-free lattice constant before composition is inferred.
2. **The PTO-vs-PZT near-degeneracy is real and quantitatively severe.** Bulk tetragonal PbTiO₃ is a≈3.904 Å, c≈4.152 Å (c/a≈1.064). Pertsev et al. (arXiv:0704.2401) report that PZT 52/48 "grow[s] closely lattice-matched on SrRuO₃-electroded SrTiO₃ up to about 40 nm in thickness (a≈0.391 nm, c/a≈1.09)" and "above 40 nm PZT thickness, the films relax to a bulk-like state (a≈0.395 nm, c/a≈1.05)" — i.e., strain alone moves c/a by ~0.04 at fixed composition. A c-axis-only measurement therefore cannot distinguish a Ti-rich PZT from a strained PTO, nor a strained film from a relaxed one of different composition.
3. **Films can be compressed below bulk c.** Lichtensteiger, Triscone, Junquera & Ghosez (*Phys. Rev. B* 71, 014103 (2005); arXiv:cond-mat/0404228) found by high-resolution XRD "a systematic decrease of the c-axis lattice parameter with decreasing film thickness below 200 Å," which they relate to "a reduction of the polarization attributed to the presence of a residual unscreened depolarizing field" (films below 50 Å show significantly reduced polarization but remain ferroelectric). So a naive "c larger than bulk = compressive in-plane strain" heuristic fails; only a measured in-plane parameter plus a relaxed reference disentangles this.
4. **Full XRR modeling (Parratt/Abelès) is the standard the fringe-spacing method only approximates.** GenX (differential evolution + Parratt) and REFL1D (Abelès optical matrix + Levenberg-Marquardt/MCMC) are the reference implementations; both use Névot–Croce roughness. Fringe spacing fails for thin/multilayer/graded/rough films and gives no density or roughness.
5. **Quantitative RSM analytics are standard in commercial tools** (PANalytical Epitaxy/AMASS, Rigaku SmartLab Studio II) and xrayutilities: relaxation line/triangle, degree of relaxation, layer tilt, and separation of mosaic (angular) from size (radial) broadening.

## Details

### AREA 1 — Strain & Composition from a Single Oriented Film (TOP PRIORITY)

**Recommended feature:** An RSM/asymmetric-reflection strain-and-composition module that (a) picks film and substrate peaks in an asymmetric RSM, (b) converts to reciprocal coordinates q∥ (q_x) and q⊥ (q_z), (c) computes in-plane a∥ and out-of-plane a⊥, (d) tests pseudomorphic vs relaxed by comparing film q_x to substrate q_x, (e) decomposes strain using a Poisson-ratio/elastic-constant biaxial model to recover the relaxed (free) lattice parameter a₀, and (f) maps a₀ → composition via Vegard's law with explicit strain/composition-confound handling.

**Method and references.** For a (001) biaxially strained pseudocubic/tetragonal film, the relaxed cubic parameter is recovered from a₀ = [(1−ν)a⊥ + 2ν a∥]/(1+ν) with ε⊥ = [(1−ν)/(1+ν)]·(a∥−a⊥)/a₀ and ε∥ = −[2ν/(1+ν)]·(a∥−a⊥)/a₀ (the standard biaxial-strain relations; see the Co₂FeAl/MgO derivation, arXiv:1305.0714, and Birkholz Ch. 7 on HRXRD strain). Equivalently ε₃ = −(2C₁₃/C₃₃)·ε₁ when single-crystal elastic constants are used. Reciprocal coordinates follow q_x = (2π/λ)(cos ω − cos(2θ−ω)), q_z = (2π/λ)(sin ω + sin(2θ−ω)). Standard references: Fewster, *X-Ray Scattering from Semiconductors*; Pietsch, Holý & Baumbach, *High-Resolution X-Ray Scattering*; Birkholz, *Thin Film Analysis by X-Ray Scattering* (Ch. 7, with Fewster). The full-unit-cell-from-reciprocal-space-vectors approach (combining 002 + 103 + 013) is given by Yang, Liu, Chen, Chen & Wang, *J. Appl. Cryst.* (2014) (arXiv:1311.3382), who recover a relaxed PZT 52/48 cell a=4.079±0.008 Å, c=4.109±0.002 Å.

**What (00l)-only misses.** Only d⊥; no in-plane parameter, no strain state, no relaxation, and a composition estimate that is entangled with strain. The PTO/PZT confound: tetragonal PbTiO₃ (c/a≈1.064) and strained Ti-rich PZT span overlapping c values, and the same PZT 52/48 composition shows c/a 1.09→1.05 across strain relaxation (Pertsev et al., arXiv:0704.2401). Films compressed below bulk c (ultrathin PTO; Lichtensteiger et al., *Phys. Rev. B* 71, 014103) further break c-based heuristics.

**Standard asymmetric reflections for perovskites on STO:** {103}/{013} pseudocubic is the workhorse (used routinely for SrRuO₃, SrIrO₃, PZT on STO); 113 and 204 are also standard. Acquiring the {103} family at φ = 0/90/180/270° additionally diagnoses tetragonal vs orthorhombic/monoclinic distortion (used for SrRuO₃: identical Q_z across the family = tetragonal; split Q_z = orthorhombic — e.g., the Sr₁₋ₓBaₓRuO₃ study, arXiv:2212.00267).

**Composition via Vegard with confound treatment.** Bulk end members: PbTiO₃ a=3.904/c=4.152 Å (JCPDS 06-0452; Shirane & Suzuki, *J. Phys. Soc. Jpn.* 7, 333 (1952); Mabud & Glazer, *J. Appl. Cryst.* 12, 49 (1979)); PbZrO₃ orthorhombic Pbam a=5.882, b=11.783, c=8.228 Å → pseudocubic ≈4.16/4.11 Å (Jona, Shirane, Mazzi & Pepinsky, *Phys. Rev.* 105, 849 (1957)); MPB at Zr/Ti≈52/48 (Jaffe, Cook & Jaffe, *Piezoelectric Ceramics*, 1971; refined to x≈0.47 by Noheda et al., *Phys. Rev. B* 61, 8687 (2000)). Critical caveat the toolkit must encode: Vegard interpolation is only piecewise-valid within a single phase field and breaks down at the MPB discontinuity (first-order tetragonal↔rhombohedral transition); and composition must be derived from the *relaxed* a₀, never from the clamped c. The group's existing pseudomorphic-strain PZT composition estimate is a special case (assumes a∥ = a_substrate); the RSM module generalizes it to partially/fully relaxed films and provides a direct relaxation check.

**Measurements required:** asymmetric RSM (or at least one asymmetric reflection) — **the group already collects RSMs and φ-scans, so no new measurement type is needed.**

**MATLAB complexity: LOW–MEDIUM.** Reciprocal-space conversion, peak picking (reuse existing prominence/profile fitting), and the closed-form biaxial algebra are straightforward. Medium effort goes into a clean elastic-constants/Poisson database (PbTiO₃ C₁₁≈237, C₁₂≈90, C₁₃≈70–100, C₃₃≈60–90 GPa, with the explicit caveat that experimental C₁₃/C₃₃ are not uniquely determined; ν≈0.3 bulk) and confound-aware composition logic.

**How leading tools do it.** PANalytical Epitaxy/AMASS extracts "strain and relaxation, mismatch and composition" from labelled peak positions; SmartLab Studio II HRXRD derives "strain, composition, and thickness" from RSM/rocking data; xrayutilities provides the reciprocal-space conversion and Materials classes underpinning this.

### AREA 2 — Thickness/Roughness/Density via Full XRR Modeling (HIGH PRIORITY)

**Where fringe-spacing breaks down.** The Δθ→thickness method assumes a single, smooth, uniform-density layer with well-resolved fringes. It fails for: very thin layers (too few fringes), multilayers (beating/overlapping periods), density gradients, and significant roughness (fringe damping), and it returns *no* density or roughness and no per-layer thickness.

**Recommended feature.** A slab-model XRR simulator/fitter using the **Parratt recursion** (or **Abelès transfer matrix**) with **Névot–Croce** interface roughness, fitting thickness, electron/mass density (SLD), and roughness per layer. Parameterization: box/slab model with d, ρ (or SLD), σ per interface; optional SLD-profile (graded) layers.

**Optimizers.** Levenberg-Marquardt for local refinement from good initial guesses (PANalytical X'Pert Reflectivity uses Parratt + LM); differential evolution for global search (GenX); MCMC/Bayesian for posterior uncertainties and to expose the well-known XRR phase-problem non-uniqueness (REFL1D). Best practice is DE/MCMC for initial global search, LM for final polish.

**References:** Parratt, *Phys. Rev.* 95, 359 (1954); Abelès (1950); Névot & Croce, *Rev. Phys. Appl.* 15, 761 (1980); Als-Nielsen & McMorrow, *Elements of Modern X-Ray Physics*; Birkholz (reflectivity chapter); Tolan, *X-Ray Scattering from Soft-Matter Thin Films*.

**Measurements required:** specular XRR — **already collected.** No new hardware.

**MATLAB complexity: MEDIUM.** Parratt recursion is a compact, vectorizable loop over layers; the engineering effort is the model-builder UI, robust optimizer wrapping (MATLAB has lsqnonlin for LM; DE and MCMC are a few hundred lines or via File Exchange), and constraint handling.

**How leading tools do it.** GenX (Parratt + differential evolution, arbitrary models, simultaneous multi-dataset); REFL1D/refnx (Abelès + LM/DE/MCMC); PANalytical X'Pert Reflectivity (Parratt recursion); SmartLab Studio II reflectivity plugin.

### AREA 3 — Microstructure from Peak Shape (Line Profile Analysis) (MEDIUM, CAVEATED)

**Recommended features (in priority order for epitaxial work):** (1) instrumental-broadening deconvolution + Scherrer for out-of-plane coherent thickness; (2) Laue/finite-thickness fringe fitting for high-quality films (often more reliable than Scherrer for the coherent thickness — the Laue period gives the number N of coherently diffracting unit cells, which times c gives the coherent thickness, while Kiessig fringes give the *total* thickness; comparing the two diagnoses dead/amorphous layers, per Miller et al., *Z. Naturforsch. B* 2022); (3) Williamson–Hall (and modified WH with dislocation contrast) and Warren–Averbach Fourier methods for size/microstrain — but only where applicable.

**Critical applicability caveat.** WH/WA are *powder* methods assuming many randomly oriented crystallites. A high-quality epitaxial film is a single (mosaic) crystal: "size" broadening is the coherent domain thickness (better captured by Laue fringes or rocking-curve/RSM analysis), and conventional WH largely fails. For epitaxial layers the physically correct route is RSM-based separation of radial (size/strain) vs angular (mosaic) broadening (see Area 5) and dislocation-density analysis via mosaic broadening (the GaN threading-dislocation literature, e.g., Kaganer et al., arXiv:cond-mat/0410510). WH/WA remain useful for polycrystalline/textured oxide-electrode films (e.g., polycrystalline PtO₂/RuO₂) the group may also grow.

**Instrumental resolution function calibration.** Use NIST SRM 660c (LaB₆) — engineered with large crystallite size and minimal defects specifically to define the instrument profile function — for the divergent-beam configuration; for HRXRD a near-perfect substrate reflection (the STO/Si substrate peak) gives the instrumental resolution directly (substrate broadening reported as ≥100× smaller than film broadening, per IOP *J. Phys. D* 56, acc597). References: Birkholz (line profile analysis chapter, with Genzel); Warren, *X-Ray Diffraction*; Williamson & Hall (1953); Warren & Averbach (1950); Ungár & Borbély (1996) for dislocation contrast; NIST SRM 660c certificate (Cline et al.).

**Measurements required:** existing θ-2θ and rocking curves; SRM 660c is a one-time calibration acquisition (new sample, not new hardware).

**MATLAB complexity: LOW (Scherrer/WH, Laue-fringe fit) to MEDIUM-HIGH (Warren–Averbach Fourier deconvolution, modified WH/WA with contrast factors).**

**How leading tools do it.** SmartLab Studio II has a crystallite-size plugin and uses SRM 660c for instrument correction; powder suites (HighScore, GSAS-II, TOPAS) implement fundamental-parameters profile fitting and WH/WA.

### AREA 4 — Superlattice/Multilayer Modeling (MEDIUM-LOW PRIORITY, HIGH EFFORT)

**Recommended feature.** A (00l) superlattice simulator beyond satellite-spacing: start with a **kinematical** sum-over-layers model (fast, valid for thin/low-reflectivity SLs) to fit individual layer thicknesses, bilayer period, and interface grading/intermixing; optionally add **dynamical (Takagi–Taupin / Darwin recursion)** for thick or strongly diffracting stacks.

**Method and references.** Kinematical structure-factor summation reproduces satellite intensities and lets a graded interdiffusion layer be fit — exactly as PANalytical AMASS does for GaN/AlN SLs, extracting per-layer d and an interdiffusion thickness d_AlGaN (arXiv:2102.06443). The dynamical standard is the Takagi–Taupin recursion for multilayers/superlattices (Bartels, Hornstra & Lobeek, *Acta Cryst. A* 42, 539 (1986)); Fewster and Pietsch/Holý/Baumbach both treat this. Per the IUCr literature, kinematical theory saves time only for an *ideal* superlattice with reflectivity below ~10% (so multiple reflections can be neglected); dynamical is needed otherwise.

**Measurements required:** existing symmetric θ-2θ around (00l) with satellites — already collected.

**MATLAB complexity: MEDIUM (kinematical) to HIGH (dynamical Takagi–Taupin with strain/grading and a fitting loop).**

**How leading tools do it.** PANalytical Epitaxy/AMASS and Rigaku SmartLab Studio II HRXRD both simulate and fit SL rocking curves/coupled scans (Bartels-type dynamical theory) to extract period, individual thicknesses, and interface quality; academic codes include the UConn MCDDM (Takagi–Taupin with dislocation/mosaic broadening), the Wisconsin kinematical nanobeam simulator, and GenX-style approaches.

### AREA 5 — Reciprocal-Space-Mapping Analytics (HIGH PRIORITY)

**Recommended features.** On top of existing load/transform/contour: (1) **relaxation line/triangle** construction (locus of fully strained point at substrate q_x vs fully relaxed point on the origin–bulk line); (2) **degree of relaxation** R from the film peak's position along that line; (3) **film tilt** from asymmetric-peak azimuthal/q_x offset; (4) **in-plane strain** extraction; (5) **separation of mosaic (angular) from size (radial) broadening** by decomposing the RLP shape along ω (tilt/lateral coherence) vs 2θ–ω (radial: vertical coherence + heterogeneous strain).

**Method and references.** Fewster's RSM reviews; Pietsch, Holý & Baumbach (RLP shape analysis — tilt and lateral coherence length from the ω-scan; vertical coherence and heterogeneous strain from the 2θ–ω scan); the relaxation-triangle construction is standard (Bruker/PANalytical app notes). Important caveat from the literature: the cubic/hexagonal "standard" relaxation formulas can give significantly wrong R for low-symmetry or miscut substrates — Kryśko et al. (arXiv:2403.02213) found degrees of relaxation "varied significantly even after 180° specimen rotation, e.g., from 8.5% to 36% for [equivalent] reflections... This makes such measurements unreliable... the hexagonal symmetry approximation is invalid." The toolkit should warn when symmetry assumptions are violated and ideally use ≥2 reflections.

**Measurements required:** RSMs — already collected.

**MATLAB complexity: LOW–MEDIUM.** Geometry/line construction and peak-shape decomposition reuse existing RSM-transform and fitting code.

**How leading tools do it.** PANalytical Epitaxy/AMASS (strain, relaxation, mismatch, mosaic spread, lateral correlation length from peak positions); SmartLab Studio II (layer tilt, relaxation degree, strain/stress, mosaicity, with non-cubic-symmetry models explicitly beyond the single-relaxation-parameter cubic assumption); xrayutilities (reciprocal-space conversion, gridding, and line-cut extraction: get_omega_scan, get_qx_scan, get_arbitrary_line for radial/angular cuts).

### AREA 6 — Reflection Geometry & Corrections (LOW EFFORT, ACCURACY-CRITICAL)

**Corrections that matter, ranked by impact on perovskite lab work:**
- **Sample displacement / zero-offset (highest impact on lattice parameters).** A height error produces a systematic 2θ shift; the group's existing Nelson–Riley refinement addresses part of this, but using the *substrate peak as an internal reference* (offset all film peaks by the known substrate position) is the single most effective lab-source correction for epitaxial films and should be built in.
- **Refraction/Kiessig correction at grazing incidence/exit** — shifts apparent peak positions for asymmetric reflections and the XRR critical edge; matters for precise in-plane parameters (Resel et al., *J. Synchrotron Rad.* 2016, on refraction/multiple-scattering effects on lattice-constant determination).
- **Kα₁/Kα₂ stripping (Rachinger).** Important for accurate FWHM and for resolving closely spaced film/substrate peaks; Rachinger (*J. Sci. Instrum.* 25, 254 (1948)) is the classic algorithm, though for HRXRD with a Ge monochromator delivering near-pure Kα₁ it is often unnecessary — the toolkit should detect/flag the configuration.
- **Lorentz–polarization and absorption.** Affect *intensities* (and thus any structure-factor or integrated-intensity work) more than peak positions; commonly omitted but needed for quantitative intensity analysis (Pike & Ladell thin-film angular-factor corrections, arXiv:1612.01805).

**Most important for precise lattice parameters and thickness:** substrate-referenced zero-offset and refraction corrections. **Commonly omitted to the detriment of accuracy:** refraction correction on asymmetric reflections and zero-offset when no internal reference is used.

**References:** Fewster (lattice-parameter accuracy, Bond method); Birkholz (corrections chapter); Rachinger (1948); Pike & Ladell.

**MATLAB complexity: LOW** for all of these (algebraic corrections); they are high-leverage relative to effort.

### AREA 7 — Uncertainty & Robustness (LOW-MEDIUM EFFORT, HIGH CREDIBILITY)

**Recommended feature.** Standardized error propagation from Poisson counting statistics into all derived quantities (peak positions, FWHM, lattice parameters, thickness, composition, relaxation). Offer two engines:
- **Analytic Jacobian/covariance propagation** — the fit already returns parameter covariances (the toolkit's profile fits report parameter errors); propagate via σ_f² = JᵀCJ for downstream quantities. Fast, exact for near-linear transforms.
- **Monte Carlo / parametric bootstrap** — resample counts with Poisson noise (or residual bootstrap), re-run the full pipeline N times, take percentile intervals. Robust for nonlinear steps (e.g., strain decomposition, XRR fits with correlated parameters, MCMC posteriors) and for non-Gaussian/asymmetric uncertainties.

**Best practice.** Use analytic covariance for cheap near-linear steps and bootstrap/MC for nonlinear or strongly correlated derived quantities; report both when they disagree (a diagnostic of nonlinearity). This is the route to defensible publication-grade error bars on c, a₀, composition, and thickness.

**References:** standard error-propagation treatments (Taylor-series Jacobian, e.g., as applied to diffraction in arXiv:2401.00412); bootstrap (Efron); MCMC for reflectivity (REFL1D/Bumps).

**MATLAB complexity: LOW** (analytic, given existing covariances) to **MEDIUM** (MC/bootstrap wrapper around the pipeline).

## Recommendations

**Stage 1 (build first — highest value/effort ratio, all reuse existing data):**
1. **RSM strain/composition module (Area 1)** — the flagship capability; closes the (00l)-only gap, generalizes the existing pseudomorphic PZT estimator, and directly serves PTO/PZT/SRO-on-STO work. *Benchmark to advance:* validate recovered a₀ and c/a against a known SRO/STO RSM and against the Yang et al. PZT 52/48 relaxed cell (a=4.079, c=4.109 Å).
2. **Quantitative RSM analytics (Area 5)** — relaxation line/degree, tilt, mosaic-vs-size broadening; small increment on existing RSM code.
3. **Uncertainty propagation (Area 7)** — analytic first, then MC/bootstrap; immediately raises every reported number to publication grade.
4. **Geometry corrections (Area 6)** — substrate-referenced zero-offset + refraction; low effort, high accuracy payoff.

**Stage 2 (build next — medium effort, high value):**
5. **Full XRR modeling (Area 2)** — Parratt + Névot–Croce, LM + differential evolution; replaces fragile fringe-spacing and adds density/roughness. *Threshold:* adopt MCMC only if reviewers/projects need posterior uncertainties or the phase-problem ambiguity bites.
6. **Laue-fringe thickness + instrumental-deconvolution Scherrer (Area 3, subset)** — reliable coherent-thickness for high-quality films; the Laue-vs-Kiessig comparison diagnoses dead layers.

**Stage 3 (build if demand justifies — higher effort or narrower applicability):**
7. **Dynamical superlattice simulation (Area 4)** — start kinematical; escalate to Takagi–Taupin only for thick/strong SLs (reflectivity >~10%). Highest implementation cost.
8. **Full Williamson–Hall/Warren–Averbach (Area 3)** — only for polycrystalline/textured oxide-electrode films; flag clearly that it is inappropriate for single-crystal epitaxial layers.

**Decision thresholds that change priorities:** if the group begins growing thick (>~40–100 nm) or partially relaxed PZT routinely, Area 1's elastic decomposition and Area 5's relaxation analytics become essential rather than nice-to-have (the Pertsev et al. PZT 52/48 data show relaxation onset near 40 nm). If multilayer/superlattice device stacks dominate, promote Area 4. If quantitative intensity/structure-factor work begins, promote the LP/absorption corrections in Area 6.

## Caveats
- **Elastic-constant uncertainty:** experimental PbTiO₃ C₁₃ and C₃₃ are not uniquely determined, and effective Poisson ratios in highly strained films deviate strongly from bulk — the *Adv. Sci.* study of highly strained PbZr₀.₂Ti₀.₈O₃/STO (ε=−1.2%; PMC8061395) reports "abnormally large negative thermal expansion (α=−1.9×10⁻⁴/°C) and Poisson's ratios (ν=0.57)." So elastic-based back-calculation of a₀ carries real uncertainty — prefer *direct* in-plane measurement (RSM) over elastic back-calculation wherever possible.
- **Vegard's law for PZT is only piecewise valid;** it breaks at the MPB (≈52/48) where a first-order tetragonal↔rhombohedral transition makes a, c, and c/a change abruptly. Composition estimates near the MPB should be flagged as unreliable.
- **WH/WA line-profile methods are powder methods** and are largely inappropriate for single-crystal epitaxial films; RSM-based mosaic/size separation is the correct epitaxial analogue.
- **XRR and RSM fits are non-unique** (phase problem; multiple SLD profiles fit the same curve) — global optimizers and uncertainty quantification are essential, not optional.
- Several quantitative film values cited (e.g., strained-vs-relaxed PZT c/a, ν=0.57) come from specific samples in the literature and are illustrative of magnitude, not universal constants.
- Commercial-tool feature descriptions are from vendor documentation and may not reflect every version; exact algorithm details (e.g., which dynamical formalism) are sometimes not publicly specified.

*Note: no new measurement types or hardware are required for the four Stage-1 recommendations or for Stage-2 XRR modeling — all reuse the symmetric scans, rocking curves, φ-scans, XRR, and RSMs the group already acquires. The only genuinely new acquisition flagged is a one-time NIST SRM 660c (LaB₆) standard for instrument-profile calibration if powder-style line-profile analysis (Area 3) is pursued for polycrystalline electrode films.*
