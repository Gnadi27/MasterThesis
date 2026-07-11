# Implementation Guide: Neumann (Prescribed-Gradient) Energy Coupling for the Li Wall Function

This is a prompt-style work order for an AI (or human) implementer. It
specifies "option 1" of the option survey in
`FINDINGS_M6Tw076_energy_coupling.md`: replace the alphat-inversion energy
coupling with a prescribed-gradient (Neumann) temperature boundary
condition, so the Li wall-function heat flux is imposed directly and the
T_1 < T_w regime (near-adiabatic walls at coarse y+) becomes representable.

Read before writing any code:

- `FINDINGS_M6Tw076_energy_coupling.md` (Sec. 2: why the alphat idiom
  fails; "Mental model" subsection: tangent vs. secant; option survey)
- `LiWallFunction/LiWallFunction.{H,C}` (shared evaluator `calcPatch`)
- `alphatLiWallFunction/alphatLiWallFunctionFvPatchScalarField.C`
  (current coupling + diagnostics — the diagnostics must survive)
- OpenFOAM v2412 sources: `heThermo::heBoundaryTypes()`,
  `gradientEnergyFvPatchScalarField`, `hePsiThermo::calculate()`,
  `externalWallHeatFluxTemperatureFvPatchScalarField` (pattern reference)
- HiSA module: `~/OpenFOAM/mgnadlinger-v2412/modules/hysa/src/modules/hisa/residualsUpdate.H`
  (energy diffusion is `fvc::laplacian(alphaEff, e)`) and
  `src/finiteVolume/jacobians/jacobianViscous/laplacian/laplacian.C`
  (boundary Jacobian treatment)

## Goal

A new boundary condition `LiWallHeatFluxTemperature` on `T` that:

1. each outer iteration, obtains q_w from `LiWallFunction::calcPatch`
   (target `qw`, per FINDINGS Sec. 3), and
2. prescribes the wall-normal temperature gradient such that the discrete
   energy flux through the wall face equals exactly -q_w (OpenFOAM sign),
   with plain non-negative diffusivities and **independent of the sign of
   (T_1 - T_w)**.

The wall's alphat patch becomes trivial (`fixedValue 0`). nut/nuTilda Li
BCs are unchanged.

## Sign conventions (get these right first)

- Paper: q_w > 0 means heat flows *into* the wall (cold wall).
- OpenFOAM: patch normal points out of the fluid; the face diffusive flux
  is `alphaEff * snGrad(he)`; the target is `qTargetOF = -q_w(paper)`.
- The prescribed gradient is the **tangent** (true wall gradient), e.g.
  M6Tw076: g_T = -q_w/kappa_w ~ -83,000 K/m. Face-value extrapolation
  along the tangent from the cell centre gives T_face ~ 224 K, *not* T_w =
  300 K. This drift is expected and is the central integration issue (see
  Pitfalls 2-4).

## Design

### New class: `LiWallHeatFluxTemperatureFvPatchScalarField`

- Derive from `fixedGradientFvPatchScalarField` (NOT mixed, NOT
  fixedValue). Reason: `heThermo::heBoundaryTypes()` inspects the T
  boundary condition and pairs a `fixedGradient`-derived T BC with
  `gradientEnergyFvPatchScalarField` on `he`, which forwards the T
  gradient to the energy field each iteration. Verify this pairing in the
  v2412 sources before relying on it (Pitfall 1).
- Dictionary entries: the shared Li `Coeffs` (kappa, B, Prt, omegaH, cv1,
  enthalpyOffset, maxIter, tolerance, strict), `Tw` (required, scalar or
  scalarField — the *physical* wall temperature, since the patch value can
  no longer carry it), `wallTransport` sub-dict (see below), plus the
  diagnostics switches currently on the alphat BC (move/duplicate them
  here; the diagnostics writer can be lifted almost verbatim).
- `updateCoeffs()`:
  1. `calcPatch(...)` -> `data.qw` (and all diagnostics fields).
  2. Compute the gradient so that the *solver's own* face assembly
     reproduces the target exactly. Do not hand-derive kappa_w; invert the
     identical expression `gradientEnergy` will apply:

     ```
     g_he,needed = qTargetOF / alphaEff_he(patch)        // same alphaEff field HiSA uses
     g_T = (g_he,needed - corr) / Cpv(patch)
     // corr = the he(p,T,patch) vs he(p,T,cells) delta term used by
     // gradientEnergyFvPatchScalarField::updateCoeffs(); copy that
     // expression verbatim so the round trip is exact by construction.
     ```

     Because alphaEff appears once in this inversion and once in the flux
     assembly, any error in the patch transport properties (they are
     evaluated at the drifted patch T) cancels — the delivered flux is
     exact regardless (bookkeeping, as with alphat, but sign-free).
  3. `gradient() = g_T`.
- Adiabatic mode (`thermalMode adiabatic`): `g_T = 0` — must reduce to
  zeroGradient behaviour identically.

### Changes to `LiWallFunction::Coeffs` / `calcPatch`

The patch T field no longer equals T_w (Pitfall 3), so every wall-state
evaluation must switch from patch values to the dictionary `Tw`:

- `Tw`, `hw`: from dict Tw. `cpw = thermo.Cp(pWall, TwField, patchi)`,
  `heWall = thermo.he(pWall, TwField, patchi)`,
  `rhoW = thermo.rhoEoS(pWall, TwField, cells)` — the (p,T,...) overloads
  already exist and are already used; only the T argument changes.
- `muw`, `lambdaW`: `basicThermo` has **no** `mu(p,T,patchi)` overload, so
  add a small `wallTransport` evaluator with modes:
  - `fromPatch` (default; current behaviour — only valid when the patch T
    is pinned at Tw, i.e. with the legacy alphat coupling),
  - `sutherland` (`As`, `Ts`, kappa via Eucken `mu*Cv*(1.32+1.77*R/Cv)` or
    via a `Pr` entry: `kappa = mu*Cp/Pr`),
  - `constant` (`mu`, `kappa`).
  Fail loudly if the new BC is used with `fromPatch`.
- Add the Tw override and wallTransport to `Coeffs` so nut/nuTilda BCs use
  the *same* wall state (all three BCs must keep identical coefficients —
  extend the existing README rule).
- `isAdiabaticMode` "auto" currently tests `Tpatch.fixesValue()`; the new
  BC does not fix a value, so "auto" would silently classify the wall as
  adiabatic. Either (a) require explicit `thermalMode isothermal;` when
  the new BC is detected, or (b) make the mode detection recognise
  `LiWallHeatFluxTemperature` and read its intent. Do (b), keep (a) as a
  fatal-error fallback for inconsistent setups.

## Pitfalls — each one must be explicitly verified, not assumed

1. **heBoundaryTypes mapping.** Confirm in v2412 that the check order is
   `fixesValue()` -> fixedEnergy, `isA<fixedGradient>` -> gradientEnergy.
   Consequently `fixesValue()` must remain `false` on the new BC. Do NOT
   "helpfully" override it to true: that pairs he with `fixedEnergy`,
   whose wall flux is `alphaEff*(he_face - he_1)/Delta` — the secant — and
   silently reinstates the entire pathology this guide exists to remove.
2. **Patch T value drift.** `gradientEnergy` extrapolates the face he
   along the prescribed gradient; `hePsiThermo::calculate()` then
   back-computes the patch T from that he for non-fixesValue patches. On
   M6Tw076 at y+_1 ~ 10 the patch T will read ~224 K while the physical
   wall is 300 K. This is a *representation artifact of the tangent*, not
   a bug — but it leaks:
3. **Consumers of patch T/transport.** Patch mu, alpha, psi, rho are
   recomputed from the drifted T. Hence the wallTransport/Tw plumbing
   above for the wall function itself; audit post-processing that reads
   wall T or patch transport (yPlus FO, Nusselt-style FOs, paraview wall
   temperatures) and document that the *physical* wall temperature lives
   in the BC dict / diagnostics output, not in the T patch values.
4. **wallHeatFlux FO.** Verify it evaluates `alphaEff*snGrad(he)` (it then
   reports the imposed q_w correctly — one of the main wins). If any
   variant computes `kappa*(T_face - T_1)/Delta` from values instead of
   snGrad, it will report garbage under this BC; check the v2412 source.
5. **HiSA residual vs Jacobian.** The explicit residual
   (`fvc::laplacian(alphaEff, e)`) uses `snGrad()` of the he patch field =
   prescribed gradient — correct automatically. The viscous Jacobian
   (`jacobianViscous/laplacian.C`) must treat gradient-type boundaries
   consistently; adiabatic (zeroGradient) walls already work in HiSA, so
   the code path exists — confirm the new BC takes the same branch. An
   inconsistent Jacobian only slows convergence (it never changes the
   converged answer), so diagnose convergence stalls here.
6. **Iteration lag.** g_T is computed from last iteration's fields
   (q_w from cell states, alphaEff from patch state). This is the same lag
   the alphat BC already has; no extra treatment needed, but do not
   "freshen" alphaEff inside updateCoeffs with hand-computed transport —
   use the registered fields so the inversion matches the assembly.
7. **Startup.** First iteration has no converged q_w; initialise gradient
   to 0 (adiabatic start) rather than a large tangent from garbage fields.

## Steps

1. Extend `Coeffs` (+`writeEntries`) with `Tw`, `wallTransport`; extend
   `calcPatch` wall-state evaluation as specified; keep behaviour
   bit-identical when `fromPatch` + legacy alphat BC are used (regression
   guard).
2. Implement `LiWallHeatFluxTemperatureFvPatchScalarField` in
   `LiWallHeatFluxTemperature/`; move the diagnostics writer into a shared
   helper so alphat- and T-BC variants both use it; register with
   `makePatchTypeField`, add to `Make/files`, `wmake libso`.
3. Unit check (extend `tests/equations` or a tiny new integration case):
   with a manufactured qTarget on one face, assert
   `alphaEff*snGrad(he) == qTargetOF` to round-off after one updateCoeffs,
   for both signs of (T_1 - T_w).
4. Integration smoke tests (`tests/integration`): isothermal + adiabatic
   variants using the new BC (T: `LiWallHeatFluxTemperature`, alphat:
   `fixedValue 0`), serial + decomposed.
5. Update `examples/boundaryConditions.foam`, `README.md` (new BC snippet,
   the "patch T is not the wall temperature under this BC" warning), and
   `DESIGN_AND_DEVIATIONS.md` (energy-coupling section).

## Acceptance criteria (validation protocol)

- **Healthy-regime regression** (strongly cooled wall, e.g. the Gao
  Ma=3 / Tw=300 K flat plate the user has under `~/Xue_et_al`): new
  coupling reproduces the legacy alphat coupling's cf and qw within ~1 %
  (both are valid there; option 1 must not regress it).
- **M6Tw076, y+_1 ~ 10** (`~/Xue_et_al/flat_plate/M6Tw076/WF/Li/`),
  compared against the wall-resolved SA reference (TM/SA) and DNS:
  - diagnostics: `h1 - hw` ~ -34 kJ/kg mid-plate (T_1 -> ~266 K), smooth
    along x, no face-to-face flip-flopping;
  - model qw -> ~2000-2600 W/m^2 (DNS: 2142; expect the ~20-25 %
    near-adiabatic amplification floor, FINDINGS Sec. 6);
  - `wallHeatFlux` FO equals -data.qw to round-off;
  - Cf within a few % of wall-resolved SA (currently +15-20 %);
  - Ch(x) smooth (currently sign-reversing dips between stations).
- **Adiabatic mode**: identical results to a plain zeroGradient T wall.
- All existing `tests/equations` pass unchanged.

## Out of scope / don'ts

- Do not modify HiSA or OpenFOAM under `/opt`.
- Do not change the momentum (nut) or SA (nuTilda) coupling logic.
- Do not remove or change defaults of the existing
  `compressible::alphatLiWallFunction` (back-compat; it remains correct
  for strongly cooled walls) — but do add the envelope guard from FINDINGS
  Sec. 5 item 3 to it.
- Do not implement signed alphat here (that is option 4, documented
  separately in FINDINGS).
