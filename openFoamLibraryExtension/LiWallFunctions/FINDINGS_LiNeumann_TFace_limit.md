# LiNeumann diagnostics findings and approximate `TFace` analysis

Context: `M6Tw076`, `Tw = 300 K`, HiSA + SA, LiNeumann wall-function cases under `WF/LiNeumann_*`.

This note collects three things:

- the run-by-run review of `postProcessing/LiNeumannWallFunctionDiagnostics`
- the interpretation of why `TFace` clamping happens
- a back-of-the-envelope reconstruction of the approximate `TFace` required to deliver the requested wall heat flux in `LiNeumann_low_yplus/Li_1_extralong`

## Executive summary

- The low-y+ LiNeumann cases up to `Li_6e-1_extralong` are only locally clamp-limited near the leading edge.
- `Li_1_extralong` converges, but the upstream third of the wall is still clamp-limited and should not be treated as validated.
- All currently available high-y+ LiNeumann cases are effectively clamp-limited over nearly the full wall and are not yet credible for wall-heat-transfer validation.
- The clamp is not on the physical wall temperature `Tw`. It is on the extrapolated patch-face temperature `TFace`, which is a bookkeeping state used to make the discrete diffusion operator deliver the Li target flux.
- In the strongly clamped region, the `TFace` needed to hit the target can fall below `100 K`, and in part of the upstream region it appears to fall below `0 K`. Lowering `TFaceMin` can help a little in the mild part of the clamp region, but cannot rescue the strongly upstream region.

## Important interpretation points

- `Tw` in the boundary-condition dictionary is the physical wall temperature used by the LiNeumann closure.
- `TFace` in the diagnostics is the extrapolated patch-face temperature seen by generic OpenFOAM/HiSA field machinery.
- In `LiNeumann`, patch `T` is therefore **not** the physical wall temperature.
- The diagnostics use OpenFOAM’s flux sign convention. In these runs, negative `qTargetOF` corresponds to heat being extracted from the boundary layer and going into the wall.

So there are two different “wall temperatures” in play:

- physical wall temperature for the model: `Tw`
- face/patch temperature used by the solver machinery: `TFace`

That distinction matters because the clamp acts on `TFace`, not on `Tw`.

## Why the clamp exists

The LiNeumann coupling computes a target wall heat flux and then solves backward for a face state that will make HiSA’s diffusive energy term deliver that flux. In simplified form:

```text
q_OF = alphaEff * deltaCoeffs * (he_face - he_1)
```

where `he_1` is the first-cell enthalpy and `he_face` is the patch-face enthalpy implied by `TFace`.

If the first cell is too far from the wall, the required face state can become extremely low even when the physical wall temperature is perfectly ordinary (`Tw = 300 K`). This is a representability problem of the coarse wall-normal spacing, not a prediction that the real wall is cold.

The floor is there because `TFace` is not a harmless algebraic ghost variable. It is still carried through real thermo/transport machinery:

- the BC inverts thermo to get `T` from the needed face enthalpy
- patch `T` is then seen by generic OpenFOAM/HiSA boundary-field logic
- negative `T` would therefore propagate into EOS/thermo/transport evaluations

So a negative `TFace` is mathematically understandable as a bookkeeping state, but numerically unsafe in the present implementation.

## Diagnostics review

### Coverage

- `LiNeumann_5e-2_extralong`: no diagnostics folder was available at the time of review
- reviewed cases:
  - `LiNeumann_1e-1_extralong`
  - `LiNeumann_2e-1_extralong`
  - `LiNeumann_6e-1_extralong`
  - `LiNeumann_1_extralong`
  - `LiNeumann_1p8_extralong`
  - `LiNeumann_3p5_extralong`
  - `LiNeumann_5_extralong`
  - `LiNeumann_10_extralong`

### Latest-time summary

| case | latest time | faces | `TFaceMin=100 K` active | furthest clamped `x` [m] | faces with `|qError|/|qTargetOF| > 10%` | wrong-sign flux faces | non-healthy faces | max relative error on unclamped faces |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `LiNeumann_1e-1_extralong` | 5000 | 617 | 9 | 0.0018538 | 10 | 4 | 0 | 9.2913 |
| `LiNeumann_2e-1_extralong` | 5000 | 617 | 15 | 0.0064840 | 15 | 8 | 1 | 0.02267 |
| `LiNeumann_6e-1_extralong` | 5000 | 617 | 26 | 0.0440517 | 26 | 17 | 0 | 0.0009468 |
| `LiNeumann_1_extralong` | 5000 | 617 | 229 | 1.0984982 | 220 | 21 | 0 | 0.0002850 |
| `LiNeumann_1p8_extralong` | 3000 | 617 | 614 | 2.9975013 | 617 | 27 | 0 | 0.6818 |
| `LiNeumann_3p5_extralong` | 3000 | 617 | 598 | 2.9975013 | 617 | 42 | 0 | 0.8680 |
| `LiNeumann_5_extralong` | 3000 | 312 | 247 | 2.9950044 | 312 | 32 | 0 | 0.9828 |
| `LiNeumann_10_extralong` | 3000 | 312 | 297 | 2.9950044 | 312 | 62 | 0 | 0.9823 |

### Reading the table

- For `LiNeumann_2e-1_extralong`, the single non-healthy face is a trailing-edge `status = 4` face with `qTargetOF = 0`, so it is not a material heat-flux problem.
- For `LiNeumann_1e-1_extralong`, the large “max relative error on unclamped faces” is driven by one tiny-target trailing-edge oddity. The main leading-edge behaviour is otherwise localized.
- For `LiNeumann_6e-1_extralong` and `LiNeumann_1_extralong`, the unclamped part of the wall is extremely clean.
- For the high-y+ cases, even many of the unclamped faces still miss the target badly, so these are not simply “a few bad leading-edge cells”.

### Practical interpretation by case

- `LiNeumann_1e-1_extralong`: basically healthy; clamp problem confined to the first ~1.9 mm of wall plus one small trailing-edge oddity.
- `LiNeumann_2e-1_extralong`: still healthy overall; clamp problem confined to the first ~6.5 mm.
- `LiNeumann_6e-1_extralong`: degraded but still mostly usable downstream; clamp problem reaches to ~44 mm.
- `LiNeumann_1_extralong`: converged, but the upstream third of the wall is still clamp-limited and not validated.
- `LiNeumann_1p8_extralong`, `LiNeumann_3p5_extralong`, `LiNeumann_5_extralong`, `LiNeumann_10_extralong`: not credible yet for wall-heat-transfer validation; the clamp-limited / large-miss region occupies essentially the full wall.

## Approximate `TFace` reconstruction for `LiNeumann_1_extralong`

### Purpose

The question was why `TFaceMin = 100 K` is not always sufficient if the physical wall temperature is `300 K` and the first-cell temperature is “around 260 K”.

The answer is that the needed `TFace` is the ghost/patch-face temperature required to produce the requested **discrete** wall flux through the coarse wall-normal spacing. It is not the physical wall temperature.

### Approximation used

Using the diagnostics at time `5000`, and approximating `he ~ cp*T` locally with `cp = 1004.703 J/(kg K)`, one can estimate

```text
q ~ alphaEff * deltaCoeffs * cp * (TFace - T1)
```

and therefore, from the actually delivered and target fluxes,

```text
TFace_needed ≈ T1 + (qTargetOF / qOpenFOAM) * (TFace_diag - T1)
```

This is only a back-of-the-envelope estimate:

- it assumes locally linear enthalpy-temperature behaviour
- it treats the diagnostic `TFace` as the active face state of the converged run
- it does not attempt to re-solve the coupled thermo inversion exactly

But it is enough to show whether the required face state lies above `100 K`, below `100 K`, or below `0 K`.

### Results at representative `x` locations

| `x` [m] | approx `T1` [K] | diagnostic `TFace` [K] | `qTargetOF` [W/m²] | `qOpenFOAM` [W/m²] | approx `TFace_needed` [K] | interpretation |
|---|---:|---:|---:|---:|---:|---|
| 0.0490 | 199.18 | 100.00 | -16863.15 | -2000.49 | -636.85 | far upstream: floor is nowhere near enough |
| 0.1990 | 271.74 | 99.99 | -5000.94 | -3463.89 | 23.78 | positive, but below `100 K` |
| 0.4988 | 267.06 | 100.00 | -5516.53 | -3369.66 | -6.44 | slightly below `0 K` |
| 0.7487 | 267.86 | 100.00 | -4643.39 | -3385.77 | 37.65 | positive, but below `100 K` |
| 0.9986 | 268.50 | 100.00 | -3922.53 | -3398.82 | 74.04 | positive, but still below `100 K` |
| 1.1984 | 268.29 | 163.30 | -3397.16 | -3397.16 | 163.30 | clean: no clamp needed |
| 1.4983 | 270.44 | 178.09 | -3228.77 | -3228.77 | 178.09 | clean |
| 1.9980 | 273.60 | 191.70 | -3053.91 | -3053.91 | 191.70 | clean |

### What this shows

- The intuition “`T1 ≈ 260 K`, so `TFace = 100 K` should be enough” only works in the milder part of the plate and only if the target flux is modest.
- Very far upstream, `T1` is not around `260 K`; at `x ≈ 0.049 m` the inferred first-cell temperature is only about `199 K`, while the target flux magnitude is very large.
- In the middle of the clamp region, the required `TFace` can still be below `100 K`, and in part of that region it appears to drop below `0 K`.

This is why lowering `TFaceMin` can only help partially:

- around `x ≈ 1.0 m`, reducing the floor below `100 K` could help
- around `x ≈ 0.75 m`, it would need to fall toward `40 K`
- around `x ≈ 0.50 m`, even `0 K` would not be enough
- around `x ≈ 0.05 m`, the required ghost state is far outside any physical temperature range

So the floor is not the core problem. The core problem is that the current coarse wall-normal spacing cannot represent the requested wall flux through a physically admissible face temperature state over the entire upstream region.

## Implications

- `TFaceMin` and `gradientRelaxation` are startup/stabilization aids, not a cure for a fundamentally under-resolved coupling.
- Lowering `TFaceMin` may recover part of the mild clamp region, but cannot rescue the strongly upstream region of `LiNeumann_1_extralong`.
- The strong trend with wall-normal spacing remains the dominant message:
  - fine low-y+ LiNeumann cases are usable
  - `LiNeumann_1_extralong` is only partly usable
  - currently available high-y+ LiNeumann cases are not validated

## Recommended practical reading of the current sweep

- acceptable / mostly acceptable:
  - `LiNeumann_1e-1_extralong`
  - `LiNeumann_2e-1_extralong`
  - `LiNeumann_6e-1_extralong`
- not acceptable as a full-wall validated result:
  - `LiNeumann_1_extralong`
- not acceptable for wall-heat-transfer validation in their current form:
  - `LiNeumann_1p8_extralong`
  - `LiNeumann_3p5_extralong`
  - `LiNeumann_5_extralong`
  - `LiNeumann_10_extralong`

## Source locations used for this note

- LiNeumann diagnostics folders under:
  - `WF/LiNeumann_low_yplus/*/postProcessing/LiNeumannWallFunctionDiagnostics`
  - `WF/LiNeumann_high_yplus/*/postProcessing/LiNeumannWallFunctionDiagnostics`
- LiNeumann wall-temperature BC setup example:
  - `WF/LiNeumann_low_yplus/Li_1_extralong/0.org/T`
- design / findings notes:
  - `/home/michael_gnadlinger/myOpenfoamLibs/LiWallFunctions/FINDINGS_M6Tw076_energy_coupling.md`
  - `/home/michael_gnadlinger/myOpenfoamLibs/LiWallFunctions/DESIGN_LiNeumann_energy_coupling.md`

---

# Addendum (2026-07-07): what the solver actually does with `TFace`, and the wall-`alphat` conditioning fix

This addendum records a follow-up discussion (during the `WangWallFunction`
implementation work) that sharpens three points of this note: exactly where
`TFace` enters the solver, why the clamp problem is largely an artefact of
`alphat = 0` at the wall rather than an intrinsic mesh-representability limit,
and why a drifted-but-physical `TFace` does not corrupt the interior solution.

## 1. `TFace` is a live thermodynamic state, not bookkeeping

Confirmed at the source: in `hePsiThermo::calculate()`
(`/opt/of-v2412/src/thermophysicalModels/basic/psiThermo/hePsiThermo.C:129-134`),
for a patch whose T BC does not fix a value (fixedGradient falls in this
branch), the boundary loop does

```cpp
pT[facei]     = mixture_.THE(phe[facei], pp[facei], pT[facei]);
ppsi[facei]   = mixture_.psi(pp[facei], pT[facei]);
pmu[facei]    = mixture_.mu(pp[facei], pT[facei]);
palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
```

So `TFace` directly sets the face compressibility `psi` (hence face density
`rho = psi*p`), face molecular `mu`, and face laminar `alpha` on every
`thermo.correct()`, and from there reaches HiSA's boundary flux assembly, the
flux Jacobian/preconditioner, and boundary contributions to `fvc::grad`/`fvc::div`
in the first cell.

The consumers split into two classes:

- **Self-cancelling by construction (while the clamp is inactive).** The energy
  coupling back-solves `heFace` against the *current* `alphaEff`, so whatever
  `alphaEff(TFace)` becomes, the delivered flux returns to the target on the
  next outer iteration (one-iteration lag only). Likewise the nut inversion
  deliberately uses the drifted patch `mu`/`rho`, so the delivered shear
  re-cancels each update.
- **Not cancelling.** Face `psi`/`rho` in the boundary linearization,
  `mu(TFace)` in the viscous wall flux between BC updates, and boundary
  contributions to near-wall gradients see the state directly.

Once the clamp binds, both classes fail at once: the delivered flux no longer
equals the target (the cancellation argument is void) *and* the face
thermodynamic state is fiction. This is the qualitatively dangerous regime.

## 2. The clamp limit is a laminar-`alphaEff` artefact: wall-`alphat` conditioning

The required face-state excursion is

```text
TFace_needed - T1  ~  qTargetOF / (alphaEff * deltaCoeffs * cp)
```

and `alphaEff` cancels identically from the delivered flux (Section "Why the
clamp exists"), so it is a **free conditioning parameter**. The Gao-pattern
`alphat = fixedValue 0` pins `alphaEff` to the laminar value
(~2.5e-5 kg/m/s); setting e.g. `fixedValue alphat 0.05` at the wall raises it
by a factor ~2000 and shrinks the required excursion by the same factor.

Quantified against the reconstruction table above (`Li_1_extralong`, time 5000):
the worst face at `x ~ 0.049 m` needs `TFace ~ -637 K` under laminar
`alphaEff`, but only `T1 - 0.4 K ~ 198.6 K` with `alphat 0.05`. This was
confirmed empirically in the Wang library smoke tests (Mach-3 flat plate,
HiSA): with `alphat 0`, delivery collapsed to 0.5% of a 600 kW/m2 target with
the clamp bound; with `alphat 0.05`, `TFace` stayed in 274-376 K and
`qError` mean was ~3 W/m2. See
`~/OpenFOAM/mgnadlinger-v2412/src/WangWallFunction/DESIGN_AND_DEVIATIONS.md`,
section "Wall alphat conditioning (important)".

**Consequence for this note's conclusions:** the statement "the coarse
wall-normal spacing cannot represent the requested flux through a physically
admissible face temperature" holds only at laminar `alphaEff`. The
clamp-limited high-y+ cases (`Li_1p8` through `Li_10`) and the upstream third
of `Li_1_extralong` are very likely rescueable by the one-line wall-`alphat`
BC change, without touching the library code.

Two clarifications on what the fix does and does not do:

- It **removes** the unphysical-state pathology rather than hiding it:
  `TFace` lands within ~1 K of `T1`, a fully admissible state, so the
  non-cancelling consumers (EoS, Jacobian, gradients) see nothing pathological.
  It also largely severs the destabilizing feedback loop
  `TFace -> alpha(TFace) -> TFace`, because `alphaEff` is then dominated by the
  prescribed constant.
- It does **not** make patch T the physical wall temperature — nothing in the
  fixedGradient coupling ever does. The residual drifted-face-state regime is
  the same one the already-clean faces operate in (e.g. `x = 1.2-2.0 m` above:
  `TFace = 163-192 K` against `Tw = 300 K`, with exact delivery).

Note also: rescuing *delivery* on the high-y+ cases says nothing about whether
the Li closure's q_w target is itself accurate at those spacings — the fix
unmasks that validation question instead of burying it under clamp artefacts.

## 3. Why a drifted-but-physical `TFace` does not corrupt the interior solution

In the finite-volume conservation formulation, the wall communicates with the
interior **exclusively through fluxes**: mass (zero at no-slip), momentum
(tau_w plus pressure), energy (q_w). "Wall temperature" never appears in the
interior equations; a hot or cold wall affects the flow only *via* the flux it
produces. If q_w and tau_w are delivered correctly, two different
(`TFace`, `alphaEff`) pairs producing the same flux are indistinguishable to
the interior solve.

Honest list of where this leaks:

1. **Second-order gradient pollution in cell 1** — the drifted face value
   enters gradient reconstructions in the first cell, perturbing higher-order
   interpolation at that cell's interior faces. Bounded, confined to the
   modeled near-wall zone, smaller than the wall-function modeling error
   itself.
2. **Face `rho`/`mu` in the linearization** — affects convergence path and
   robustness, not the converged flux balance.
3. **Anything that reads patch T directly as "the wall temperature"** — this
   is the real trap: Stanton numbers / heat-transfer coefficients computed
   from patch T, radiation models sampling wall T^4, conjugate heat transfer
   (coupled patches exchange T), or any other BC/function object consuming
   patch T would be silently wrong. The physical `Tw` lives only in the BC
   dictionary; wall quantities must be extracted from the model diagnostics,
   never from patch T. (The wall-function closures themselves are immune by
   design: they evaluate wall properties from dictionary `Tw` +
   `wallTransport`, precisely because patch fields drift.)

### Contaminated wall-patch fields (state-type vs flux-type outputs)

The `TFace` contamination propagates through `hePsiThermo::calculate()` into
the **entire thermo patch-field family on the wall patch** — everything below
is evaluated at `TFace` instead of `Tw` (interior cells, including cell 1, are
unaffected; and patch `p` is trustworthy, it comes from the pressure BC):

| patch field | how contaminated | typical error at `TFace ~ 190 K`, `Tw = 300 K` |
|---|---|---|
| `T` | is `TFace` itself | ~110 K low |
| `psi` | `psi(p, TFace) = 1/(R*TFace)` | ~60% high |
| `rho` | `psi*p = p/(R*TFace)` | ~60% high |
| `mu` | Sutherland at `TFace` | ~30% low |
| `alpha` (laminar) | `alphah(p, TFace)` | ~30% low |
| `nut` | fiction by design: absorbs the drifted `mu`/`rho` so the *total* stress is right | split mu/rho*nut meaningless |
| `alphat` | prescribed conditioning constant, not physics | meaningless as a state |

Derived / post-processed quantities that inherit the contamination if built
from wall-patch fields:

- Stanton number or heat-transfer coefficient using patch `T` or patch `rho`
- skin-friction coefficient `Cf` normalized with patch `rho` (use freestream
  or dictionary-`Tw` wall density instead)
- a hand-computed `y+`/`y*` from patch `nu = mu/rho` (doubly wrong: `mu` low,
  `rho` high; use the nut BC's `yPlus()` output, which is computed from the
  physical wall state)
- recovery factor / adiabatic-wall-temperature estimates read off patch `T`
- sampled wall surface fields (`surfaces`/`boundaryCloud` of `T`, `rho`, `mu`)
  fed to external tools
- any mass-flux-like normalization using wall `rho*U_tau`

Flux-type function objects remain **correct** despite the drift:

- `wallHeatFlux` = `alphaEff*snGrad(he)` — this *is* the delivered target, by
  construction of the back-solve
- `wallShearStress` — the nut inversion calibrates the total
  `(mu + rho*nut)*snGrad(U)` to the model tau_w, so the total is right even
  though the mu / rho*nut split is fiction

Recovery rule: reconstruct wall-state quantities from dictionary `Tw` and
patch `p` (e.g. `rho_w = p/(R*Tw)`, `mu_w = mu_Sutherland(Tw)`), or read them
from the library diagnostics — never from the wall-patch thermo fields.

### Why patch `p` is the safe anchor

Patch `p` escapes the contamination for two independent reasons:

1. **The dependency at the patch flows one way, and `p` is upstream of it.**
   The wall `p` BC (`zeroGradient` in these cases) just copies the interior
   cell value, `p_face = p_1`; it is never evaluated by the thermo. In the
   `hePsiThermo::calculate()` patch loop, `psi`/`mu`/`alpha` are computed
   *from* patch `(p, T)` — pressure is an input to the contaminated
   evaluations, never an output. There is no code path in which face `rho`,
   `mu`, or `psi` is used to compute face `p`, so the contamination cannot
   flow backwards into it.
2. **`p_face = p_1`, and the interior (including `p_1`) is protected by the
   flux argument** (Section 3 above): the wall's flux contributions are mass
   zero, pressure force `p_face*n` (itself interior data), calibrated tau_w,
   calibrated q_w. Face density has almost nowhere to act at a wall even in
   principle — it multiplies the convective flux, which is identically zero
   through a no-slip face except for the pressure term. So `p_1`, hence
   `p_face`, is clean to the same discretization order as the rest of the
   interior, with no `TFace`-specific error on top.

**Where this firewall breaks** — the "`p` is safe" statement is specific to
the present setup (density-based HiSA, `zeroGradient` `p` at a static no-slip
wall). Configurations in which the contaminated face state *would* reach the
pressure boundary or the flux balance:

- **Pressure-based buoyant solvers**: `p_rgh` BCs such as `fixedFluxPressure`
  explicitly consume the face density in their `rho*g*h` terms — a ~60%-high
  face `rho` would contaminate the pressure boundary directly.
- **Moving walls**: the boundary viscous-work term `tau . U_wall` is nonzero
  and uses face `mu`; at a static no-slip wall it vanishes identically.
- **Radiation / CHT** (already listed above): wall `T^4` sampling and coupled
  patch T exchange consume the drifted patch T directly.
- **Any custom `p` (or other) BC that references patch `rho` or patch `T`.**

If the Neumann coupling is ever reused outside the HiSA static-wall context,
these paths must be re-audited before trusting patch `p` — or any interior
field — as uncontaminated.

For a pure aero-thermal RANS solve (no radiation, no CHT, wall quantities from
diagnostics), the interior is protected to discretization order.

**Falsifiable check** (recommended before rolling the `alphat` fix across the
sweep): since `alphat` is pure conditioning, run one case twice with e.g.
`alphat 0.05` and `alphat 0.5`. `TFace` will differ substantially between the
runs (excursion scales as `1/alphaEff`); if the arbitrary face temperature
leaked into the physics, the converged fields would differ. Matching q_w,
tau_w and interior profiles to convergence tolerance demonstrates directly
that the interior never sees `TFace`.

## Cross-references

- `~/OpenFOAM/mgnadlinger-v2412/src/WangWallFunction/DESIGN_AND_DEVIATIONS.md`
  ("Wall alphat conditioning (important)") — quantified smoke-case numbers.
- Sign-convention warning while touching the sibling libraries: Gao/Nichols
  q_w is positive fluid->wall (OpenFOAM target `-q_w`); Wang/Craft q_w is
  positive wall->fluid (target `+q_w`). Never copy the sign mapping between
  the families.

# Addendum (2026-07-10): HiSA source audit — which solver paths actually consume `TFace`

The 2026-07-07 addendum established that `TFace` sets the face `psi`, `mu`
and `alpha` on every `thermo.correct()`. This addendum traces where HiSA
(`~/OpenFOAM/mgnadlinger-v2412/modules/hysa`) then *consumes* those face
values, path by path. Headline: HiSA already special-cases wall patches in
most of the acoustically dangerous places, so the exposure is narrower than
the previous addendum's "reaches HiSA's boundary flux assembly, the flux
Jacobian/preconditioner" suggests — but the surviving paths are exactly the
ones that matter for delivered shear and implicit conditioning.

`TFace` enters through `updateFields.H:186` (`T.correctBoundaryConditions()`).
HiSA then rebuilds the face states from patch `T`: face energy
`e = he(p, T)` (line 194), `psi` via `thermo.correct()` (line 196), face
`rho` and `rhoE` (lines 198-225). Note HiSA has its own solver-level patch-T
clamp at lines 188-192 (`TMin`/`TMax` in the `pseudoTime` dict), but the
defaults are `ROOTSMALL`/`GREAT`, so the BC-level `TFaceMin`/`TFaceMax` is
the operative fence in practice.

## Paths verified `TFace`-INDEPENDENT (no impact)

| Path | Why it is safe | Source |
|---|---|---|
| Local pseudo-time step at walls | wall faces use the *cell* sound speed and velocity, not the face state | `hisaModuleTimestepping.C:96-120` |
| Wall face density (with `characteristicWallPressure` on `p`) | HiSA substitutes the *cell* density on those patches | `updateFields.H:199-217` |
| `characteristicWallPressure` itself | builds `rho*c` from `patchInternalField()` (cell) values | `characteristicWallPressureFvPatchScalarField.C:196-213` |
| Wall inviscid flux (AUSM+up etc.) | reconstruction does set L = R = patch value at walls (`reconstructionScheme.C:254-258`), but at a no-slip face `un_L = un_R = 0` and `p_L = p_R`, so every `TFace`-dependent term (mass flux, low-Mach p-diffusion, velocity diffusion, energy flux) carries an exact zero factor; the survivor is `p_patch*Sf` | `ausmPlusUpFluxScheme.C:82-241` |
| Inviscid boundary Jacobian | same cancellation: all `d(flux)/d(p,T)` sensitivities carry `U_B & Sf = 0`, and the `d/dU` terms chain through the no-slip BC's zero `valueInternalCoeffs` | `laxFriedrichs.C:177-273` |
| Viscous implicit treatment | `boundaryJacobian` returns zeros and the viscous stab excludes physical boundaries (`includePhysicalBoundaries` defaults `false`) — viscous wall terms are explicit-only | `spectralRadius.C:89-112`, `fvjOperators.H:127` |
| **Delivered wall shear** | the nut inversion deliberately divides by the solver's patch `mu`/`rho` (`LiNeumannWallFunction.C:948-966`, comment at 944-947), so `(muPatch + rhoPatch*nut)*|gradU| = tauw` identically — the same cancellation as the energy inversion | see also "one-sided failure mode" below |
| Wall heat flux target | the closure computes `q_w` from `wallTransport` at the physical `Tw`; `TFace` never enters the target | `LiNeumannWallFunction.C:620-676` |

Consequence for **accuracy/general flow**: with local time stepping,
`characteristicWallPressure` wall pressure, no-slip walls, and both Li
inversions cancelling, a drifted-but-clamped `TFace` does **not** leak into
the interior solution through the inviscid discretization, the delivered
shear, or the delivered heat flux. This independently confirms (and
mechanistically explains) the falsifiable `alphat 0.05` vs `0.5` check above.

## Paths verified `TFace`-DEPENDENT

1. **Implicit Lax-Friedrichs stabilization diagonal — performance, not
   accuracy.** `laxFriedrichs.C:129` passes `includePhysicalBoundaries=true`,
   so wall faces add `0.5*(c_face + |U_face|)*magSf` to the wall-cell diagonal
   (`fvjOperatorsTemplates.C:711-731`) with
   `c_face = sqrt(gamma/psi_patch) = sqrt(gamma*R*TFace)`
   (`laxFriedrichs.C:326-333`). The Jacobian only steers the pseudo-transient
   path, so the converged solution is untouched; a high `TFace` over-damps
   wall-cell updates (slower convergence), a low clamp (`TFaceMin ~ 1 K`)
   under-damps them by up to `sqrt(T1/TFaceMin)`. Impact: convergence
   rate/robustness only.

2. **One-sided nut clip — the only accuracy leak, upper side only.** The
   shear cancellation requires `nut >= 0`. When `mu(TFace)*|gradU| > tauw`
   (Sutherland: `mu(10000 K)/mu(300 K) ~ 8`), the inversion clips at zero and
   the delivered shear **exceeds** the target by exactly
   `mu(TFace)*|gradU| - tauw`, uncompensated. This only engages when `TFace`
   drifts *high* — the cold-wall clamping regime of this note (low `TFace`,
   deflated `mu`) is safe because `nut` simply compensates upward. It can
   also engage physically near leading edges where the laminar term alone
   carries most of `tauw`. This is the momentum-side argument for a real
   `TFaceMax` on heated-wall cases. Diagnostics: see below.

3. **SA/turbulence wall-face transport — small accuracy exposure.** The
   `nuTilda` diffusion coefficient and anything reading wall
   `nu = mu(TFace)/rho_face` from the patch sees the drifted state; no
   inversion compensates it. Second order (wall-face diffusion coefficient of
   one equation), but nonzero.

4. **Global (non-local) timestep branch** (`hisaModuleTimestepping.C:131`)
   uses interpolated face values including `c(TFace)`. Only relevant if
   `localTimestepping` is off — performance impact only (timestep tightened
   by the max face wave speed).

5. **Post-processing.** The genuinely drifted outputs are the raw thermo
   patch fields (`T`, `rho`, `psi`, `mu`, `alpha`) and anything hand-derived
   from them (Stanton number from patch `T`, `Cf` normalized by patch `rho`,
   a hand-computed `y+` from patch `nu`, sampled wall surfaces of
   `T`/`rho`/`mu`), consistent with the 2026-07-07 addendum. The standard
   function objects, however, are protected — an earlier revision of this
   item wrongly listed `wallShearStress` and `yPlus` among the drifted
   outputs; verified against the source (2026-07-10), both are correct:

   - **`wallShearStress` FO — correct (total), as long as the `nut >= 0`
     clip is inactive.** The FO recomputes stress from patch
     `devRhoReff()`, i.e. from `muEff = mu(TFace) + rho(TFace)*nut`
     (`/opt/of-v2412/src/functionObjects/field/wallShearStress/wallShearStress.C:198`).
     That split *is* the wrong drifted state — `mu(TFace)` is ~30% low and
     `rho(TFace)` ~60% high in the cold-clamp regime, and the mu / rho*nut
     decomposition is fiction. But this has no effect on the reported total,
     because the nut inversion back-solved `nut` against the *same* drifted
     patch `mu`/`rho` so that `(muPatch + rhoPatch*nut)*|snGrad(U)| = tauw`
     identically — the drift cancels in exactly the combination the FO
     evaluates, the same algebraic cancellation that protects the delivered
     shear itself. The cancellation fails only when the inversion clips at
     `nut = 0` (`mu(TFace)*|gradU| > tauw`, i.e. `TFace` drifted *high* or
     near-laminar leading-edge faces): then the FO faithfully reports the
     over-delivered stress, not the model target. In the cold-wall sweep the
     clip cannot engage (low `TFace` deflates `mu`); the `tauError` /
     clipped-faces diagnostics below flag it when it does. Residual caveats:
     one-outer-iteration lag during transients (exact at convergence), and
     the FO's full `devRhoReff` tensor vs. the calibration's plain
     `muEff*snGrad(U)` form — a small generic wall-function discrepancy that
     scales with the calibrated total, not a `TFace` contamination.
   - **`yPlus` FO — correct.** The FO delegates to the nut BC's own
     `yPlus()` for any patch whose nut BC
     `isA<nutWallFunctionFvPatchScalarField>`
     (`/opt/of-v2412/src/functionObjects/field/yPlus/yPlus.C:152-157`);
     `nutLiNeumannWallFunction` derives from that base and returns the
     core's `data.yPlus = rhoWall*uTau*yFace/muWall`
     (`LiNeumannWallFunction.C:903`), built from `rhoEoS(pWall, Tw)` and
     Sutherland-at-`Tw` — the physical wall state, no `TFace` anywhere. Only
     the FO's fallback branch for wall patches *without* a wall-function nut
     BC computes from patch `nu` and would report the drifted state.

   Flux-type outputs mirrored by the BC (`wallHeatFlux` via the `he`
   gradient) remain correct. The physical wall temperature lives only in
   `Tw` (dictionary/diagnostics), never in the `T` patch values.

## New tau diagnostics (implemented 2026-07-10)

The diagnostics writer now emits, per face, alongside the q columns:

- `magGradU` — tangential `|snGrad(U)|` as used by the nut inversion.
- `tauDeficit` — `min(0, tauw - muPatch*magGradU)` recorded *inside* the
  inversion (`LiNeumannWallFunction.C`): negative exactly where the
  `nut >= 0` clip engaged, with the over-delivery margin in stress units.
  Lag-free and unambiguous.
- `tauOpenFOAM` — reconstructed delivered stress
  `(muPatch + rhoPatch*nut)*magGradU` with the current patch state.
- `tauError` — `tauOpenFOAM - tauw`. Because the cancellation is algebraic,
  this is identically zero wherever the clip is inactive (modulo
  one-outer-iteration lag), so a persistent nonzero `tauError` *is* the clip
  signature — unlike `qError`, it cannot be polluted by relaxation.

The `diagnosticLog` summary prints `tauTarget`/`tauOpenFOAM`/`tauError`
min/max/mean and a `nut clipped faces = N / total` counter. Expected
signature in the current cold-wall M6 sweep: zero clipped faces everywhere
(the clip cannot engage from the low-`TFace` side); the counter earns its
keep on heated walls, high-`TFace` excursions, and near-laminar leading-edge
faces.
