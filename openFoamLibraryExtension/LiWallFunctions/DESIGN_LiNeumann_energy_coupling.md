# LiNeumann Wall Functions - Neumann Energy Coupling Design Notes

This documents the `LiNeumann*` boundary-condition family, the "option 1"
rework specified by `IMPLEMENTATION_GUIDE_OPTION1_neumann_energy_BC.md` in
response to `FINDINGS_M6Tw076_energy_coupling.md`. It is a copy of the Li
frozen-composition wall function (see `DESIGN_AND_DEVIATIONS.md`, which
remains authoritative for the closure itself) with one change: the wall
heat flux is imposed directly through a prescribed-gradient temperature
boundary condition instead of being packaged into a non-negative `alphat`.
The original `*LiWallFunction` conditions are unchanged.

## Classes

| field  | boundary condition | notes |
|--------|--------------------|-------|
| `T`       | `LiNeumannWallHeatFluxTemperature` | derives from `fixedGradient`; computes q_w each outer iteration and prescribes the wall-normal T gradient |
| `alphat`  | plain `fixedValue 0`               | no Li alphat class in this variant |
| `nut`     | `nutLiNeumannWallFunction`         | identical momentum logic to `nutLiWallFunction` |
| `nuTilda` | `nuTildaLiNeumannWallFunction`     | defines the identical Eq. (39) callback; stock v2412 SA does not invoke it without an explicit `boundaryManipulate()` call |

Shared evaluator: namespace `LiNeumannWallFunction`
(`LiNeumannWallFunction/`), a copy of `LiWallFunction` extended with the
wall-state plumbing below. All conditions on a patch must carry identical
coefficients, including `Tw` and `wallTransport`.

## Why the flux can now be delivered for T_1 < T_w

The discrete wall energy flux remains `alphaEff * snGrad(he)` with plain
non-negative diffusivities, but `snGrad(he)` is now the *prescribed*
gradient (the tangent), not the resolved cell-to-wall secant. A cooled-wall
flux therefore coexists with a first cell colder than the wall; the
pinned-`T_1 = T_w` attractor of the alphat inversion cannot form.

## The inversion (deviation from the implementation guide)

The guide prescribed inverting `gradientEnergyFvPatchScalarField`'s
expression (`g_he = Cpv*g_T + mixture-corr`). Inspection of the actual
v2412 + HiSA code path shows that expression is never exercised by HiSA:

- `createFields.H`: HiSA's `e` *is* `thermo.he()`.
- `updateFields.H`: each outer iteration HiSA calls
  `T.correctBoundaryConditions()` (which runs this BC and extrapolates the
  face T along the prescribed gradient), then **clobbers the he boundary
  values** with `e.boundaryFieldRef() == thermo.he(p, T)->boundaryField()`;
  `e.correctBoundaryConditions()` is never called, so
  `gradientEnergy::updateCoeffs()` never runs.
- `residualsUpdate.H`: the energy diffusion is
  `fvc::laplacian(alphaEff, eCalc)` with `eCalc = rhoE/rho - 0.5|U|^2`, a
  calculated-boundary field whose face values reduce to `he(p_f, T_f)`
  minus a uniform reference (`eZero`, uniform, cancels in `snGrad`).

The flux HiSA assembles is therefore the secant through the *evaluated*
face state:

```text
q_OF = alphaEff * deltaCoeffs * (he(p_f, T_f) - he_1),  T_f = T_1 + g_T/deltaCoeffs
```

The exact-by-construction inversion places the face enthalpy directly:

```text
he_f,needed = he_1 + qTargetOF/(alphaEff*deltaCoeffs)
T_f,needed  = THE(he_f,needed, p_f)          // thermo's own inverse
g_T         = (T_f,needed - T_1)*deltaCoeffs
```

This is exact for arbitrary Cp(T) under HiSA's assembly and identical to
the guide's Cpv formula for a calorically perfect gas. Under a solver that
does run `gradientEnergy` (e.g. `rhoSimpleFoam`), the delivered gradient is
its Cpv translation of the same `g_T`; the two differ only by the second-
order Cp variation across the face-cell enthalpy gap, and coincide for
`hConst`. `alphaEff` appears once in the inversion and once in the
assembly, so its evaluation at the drifted patch state cancels
(bookkeeping, as with alphat, but sign-free). All quantities lag one outer
iteration; at convergence the delivered flux equals `-data.qw` exactly
(unless the face-T clamp or relaxation is active, below).

## wallHeatFlux reporting

`wallHeatFlux` evaluates `alphaEff * heBf.snGrad()`; for a
`gradientEnergy` patch `snGrad()` returns the stored gradient, which HiSA
never refreshes. `updateCoeffs()` therefore mirrors the delivered energy
gradient onto the he patch (`fixedGradient::gradient()`), so the function
object reports the imposed q_w. Solvers that do run `gradientEnergy`
overwrite the mirror with the equivalent value; both paths are consistent.

## Wall state: Tw and wallTransport

`heBoundaryTypes()` pairs any `fixedGradient`-derived T condition with
`gradientEnergy` on he (verified in v2412 `basicThermo.C`); the condition
must NOT fix a value. Consequently the patch T value is the tangent-
extrapolated face temperature — on M6Tw076 at y+_1 ~ 10 it reads ~224 K
against a 300 K wall. This is a representation artifact, not a bug, but it
leaks into everything evaluated at patch T (patch mu, alpha, psi, rho,
yPlus-style post-processing, ParaView wall temperatures). Hence:

- `Tw` (dictionary, required for isothermal walls) supplies the physical
  wall temperature; `hw`, `Cp_w`, `rho_w`, `he_w` are evaluated at
  `(p_wall, Tw)` through the existing thermo overloads.
- `wallTransport` supplies `mu_w`/`lambda_w` since `basicThermo` has no
  `mu(p,T,patchi)` overload: `sutherland` (`As`, `Ts`; `kappa` from `Pr`
  if given, else modified Eucken `mu*Cv*(1.32 + 1.77*R/Cv)`) or `constant`
  (`mu`, `kappa`). The default `fromPatch` is refused for isothermal walls
  under this coupling (fails loudly) and remains valid for adiabatic mode,
  where the floating patch T is physical.
- The momentum inversion to `nut` deliberately keeps the *patch-field*
  `mu` and `rho` (the values OpenFOAM's stress operator uses), while
  `uTau`, `yPlus`, `Pr_w`, and `q_w` use the physical wall state.
- `thermalMode auto` recognises the new T condition (its
  `adiabaticIntent()`), so nut/nuTilda do not silently misclassify the
  wall as adiabatic.

The physical wall temperature is recorded in the diagnostics output
(`Tw` column) and in the boundary dictionary — not in the T patch values.

## Flux target

`thermalFluxTarget auto` resolves to `qw` (the physically and discretely
consistent target for a boundary face at y=0, FINDINGS Sec. 3). `qHalf`
remains selectable as a legacy/diagnostic target and is still rejected for
adiabatic walls. Adiabatic mode prescribes `g_T = 0` and is identical to
`zeroGradient` (verified in the integration tests).

## Numerical guards and startup

- First start (no `gradient` entry in `0/T`): gradient initialises to
  zero — an adiabatic start — instead of a tangent from garbage fields.
- `TFaceMin`/`TFaceMax` (defaults 1/GREAT) clamp the extrapolated face
  temperature; without a floor a violent startup transient can demand a
  negative face temperature and crash the thermo. While a clamp is active
  the delivered flux undershoots the target (visible as `qError` in the
  diagnostics).
- `gradientRelaxation` (default 1 = off) under-relaxes the gradient
  update.
- Invalid faces (laminar fallback, `status != 0`) carry `q_w = 0`, so the
  face relaxes toward the adiabatic state rather than extrapolating from
  garbage.

A prescribed gradient enters the energy matrix as an explicit flux with no
diagonal feedback (`gradientInternalCoeffs() = 0`; HiSA's viscous Jacobian
takes the same already-working branch as zeroGradient walls). Two
consequences observed on the deliberately crude integration toy case:
startup overshoot (mitigate with `gradientRelaxation` + `TFaceMin`, or the
solver's own pseudo-transient damping), and, when the delivered cooling
drives `T_1` below `T_w`, target flip-flopping through the pre-existing
`rg <= 0` guard (heated-wall states fall back to `q_w = 0`). The latter is
the known heated-wall limitation of the underlying implementation, not a
property of the coupling.

## Diagnostics

Same writer as the alphat variant, moved onto the T condition, written to
`postProcessing/LiNeumannWallFunctionDiagnostics/<patch>/<time>.dat` with
the added columns `gradT TFace alphaEff heGradient` replacing `alphat
alphaEff heSnGrad`. `qOpenFOAM = alphaEff*heGradient` is the flux the
assembly will deliver; `qError = qOpenFOAM - qTargetOF` is nonzero only
where the clamp or relaxation intervened.

## Validation status

- Equation tests (`tests/equations`): LiNeumann law evaluations match the
  Li originals; coefficient parsing; Neumann inversion round trip for both
  signs of `(T_1 - T_w)`. All pass.
- Integration smoke tests (`tests/integration/Allrun`): isothermal,
  adiabatic (reduces to zeroGradient exactly), and 2-rank decomposed runs
  of the new conditions, alongside the unchanged legacy cases. All pass.
- The acceptance runs of the implementation guide (Gao Ma=3 healthy-regime
  regression against the alphat coupling; M6Tw076 at y+_1 ~ 10 against
  wall-resolved SA and DNS) are pending — they require the HiSA case
  family under `~/Xue_et_al` and are the intended next step.
