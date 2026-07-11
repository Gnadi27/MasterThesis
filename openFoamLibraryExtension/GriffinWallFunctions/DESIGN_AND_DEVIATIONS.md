# Griffin wall function: design decisions and deviations from the papers

Reference papers (PDFs in `../GriffinWallFunction/`):

- **[WM]** Griffin, Fu, Moin, "Near-wall model for compressible turbulent
  boundary layers based on an inverse velocity transformation", CTR ARB
  2022 — the wall model itself (Eqs. 2.1-2.8).
- **[VT]** Griffin, Fu, Moin, PNAS 118 (2021) — total-stress-based
  velocity transformation (Eqs. 1-7).
- **[BL]** Griffin, Fu, Moin, PRF 6 (2021) — boundary-layer-thickness
  local reconstruction; Appendix D = compressible (stagnation-enthalpy)
  form.
- **[CE]** Griffin, Fu, Moin, CTR ARB 2021 (grid-point requirements) —
  background only; note its T-U relation variant (sPr = 0.8, r = 0.9) was
  NOT used: [WM] Eq. 2.8 with s = 1.14, r = Pr^(1/3) is authoritative.

## What is transcribed faithfully

- Inverse transformation [WM] Eq. 2.5 with the constant-stress-layer
  form of [VT] Eq. 5; the property-gradient terms are evaluated
  analytically via dT/dU of Eq. 2.8 (chain rule), not by finite
  differences.
- Van Driest-damped mixing length [WM] Eq. 2.2 in semi-local units and
  the algebraic incompressible shear S_t+ (no ODE, as the paper notes).
- Calibration table kappa/A+ at Re_tau = {182, 543, 1000, 1990, 5190}:
  {0.400, 0.408, 0.400, 0.391, 0.391} / {18.2, 17.4, 17.0, 16.5, 16.5},
  linear interpolation, constant extrapolation, indexed by
  Re_tau* = sqrt(rho_e tau_w) delta / mu_e ([WM] Sec. 2).
- T-U relation [WM] Eq. 2.8 exactly, including the matching-point
  correction term (U/U_m)^2 (T_m - T_w); s = 1.14, T_r = T_e +
  Pr^(1/3) U_e^2/(2 Cp). Endpoints are exact: T(0) = T_w, T(U_m) = T_m.
- Solution procedure [WM] Sec. 2.3: shooting on tau_w with per-node
  sub-iteration for the T(U)-consistent properties.
- q_w = s tau_w Cp (T_w - T_r)/U_e ([WM] Sec. 2.2). Adiabatic:
  q_w = 0 <=> T_w = T_r.
- BL edge: [BL] Appendix D. h_0 = h + |U|^2/2, h_0,ref = profile max,
  U_I^2 = 2(h_0,ref - h) - V^2, edge at the first wall-outward crossing
  of (U/U_I)^2 = threshold^2 (squared form; threshold 0.99), searching
  from the wall outward as the paper recommends. Valid across shocks.

### Sign convention (important)

Griffin's q_w (via B_q = q_w/(Cp T_w rho_w u_tau); cooled walls have
-B_q > 0) is **positive for heat entering the fluid**. That coincides
with OpenFOAM's outward-normal boundary-flux sign
(alphaEff*snGrad(he) > 0 = heating the fluid), so `qTargetOF = data.qw`
with **no sign flip** — unlike the Li family (q_w positive into the
wall, target `-q_w`). Never copy sign mappings between the sibling
libraries; the cooled-wall unit test pins this down.

## Deviations and assumptions (flagged)

1. **Matching location = first-cell centre**, not the paper's WMLES
   y_m = 0.3 delta. This is the RANS wall-function idiom (the discrete
   BC lives at the first cell); the paper only requires the matching
   data to lie in the inner layer, which a wall-function-spaced first
   cell satisfies. A one-time warning is emitted if y_1 > 0.3 delta.
2. **RANS turbulence coupling is not specified by the paper** (in WMLES
   only tau_w, q_w feed back). Adopted from the validated Li pattern:
   - nut BC inverts the discrete gradient so the delivered wall shear
     equals tau_w exactly (patch-field mu/rho used deliberately — they
     are drifted, and cancel);
   - nuTilda = 0 at the wall + first-cell pin via `matrix.setValues` to
     the value whose SA nut equals the Griffin eddy viscosity at the
     first-cell centre, mut(y_1) = (tau_w - mu dU/dy)/(dU/dy)
     (dimensional constant-stress form of [VT] Eq. 6), inverted through
     fv1.
3. **Energy coupling is Neumann-only** (fixedGradient T BC paired with
   gradientEnergy by `basicThermo::heBoundaryTypes()`; the alphat patch
   is a plain fixedValue). This matches the paper's direct-flux coupling
   semantics and deliberately abandons the OpenFOAM alphat-inversion
   idiom, which is structurally unable to represent T_1 < T_w
   (near-adiabatic walls at coarse y+); see
   `../LiWallFunctions/FINDINGS_M6Tw076_energy_coupling.md`.
4. **BL edge via single-processor cell walks.** The walk (greedy
   max-inward-advancement over `cellCells()`) cannot cross processor
   boundaries; affected faces fall back to the outermost sample (warned)
   — decompose wall-tangentially or use `edgeMode prescribed`. Two
   further caveats inherent to the method:
   - with near-constant stagnation enthalpy along the profile (adiabatic
     outer flow) U_I tracks the profile's own maximum, so a walk that
     stops inside the boundary layer produces an *underestimated* delta
     with a "successful" crossing — size `maxEdgeCells` generously;
   - on coarse outer meshes the crossing interpolation produces a
     face-to-face sawtooth in U_e/delta of a few percent (visible in the
     smoke cases).
5. **Constant-stress-layer form tau+ = 1** of the transformation ([VT]
   Eq. 5) — the papers' own practical recommendation ("no discernible
   difference" vs. the exact total-stress form in their data).
6. **Calorically-perfect-gas envelope.** The T-U relation and the
   isobaric ideal-gas density rho(T) = rho_1 T_1/T (thin-BL constant
   pressure) assume constant Cp and R; Cp is evaluated from the thermo
   at the first-cell state. Strongly varying-Cp thermo (janaf at high
   temperature ranges) is outside the faithful envelope. The edge
   reconstruction itself uses h from the thermo directly (general).
7. **Viscosity from an explicit dictionary law** (Sutherland / power
   law / constant), not from the thermo: the inner integration needs
   mu(T) at arbitrary T, and the patch transport fields are drifted
   under the Neumann coupling. Keep it consistent with
   `thermophysicalProperties` yourself.
8. **Model constants dict-overridable** (s, Pr, edgeThreshold) with
   paper defaults; the kappa/A+ calibration table is fixed in code (it
   is part of the model).
9. **Numerics of the inner solve** (not specified in the paper beyond
   "shooting"): log-spaced grid (48 intervals) from a viscous seed point
   (y+ <= 0.05, U = tau_w y/mu_w) to y_1, RK2 midpoint steps, <= 5
   fixed-point property sub-iterations per node, bracketed bisection on
   tau_w (60 iterations, relative tolerance 1e-6). Unit-tested: the
   incompressible limit recovers a manufactured tau_w to 0.2% and is
   grid-converged at 48 nodes to 0.5%.

## Wall alphat conditioning (important)

alphaEff appears once in the gradient back-solve and once in the flux
assembly, so it cancels from the delivered flux — it is a free
conditioning parameter. But it scales the face-state excursion

    TFace - T_1  ~  qTargetOF / (alphaEff * deltaCoeffs * cp),

and with `alphat = 0` (laminar alphaEff, ~1e-5 kg/m/s at cold walls /
low pressure) the required TFace can leave the physical range, bind the
`TFaceMin` clamp, and silently collapse the delivered flux (the
cancellation argument is void once clamped). Prescribe a constant wall
alphat, e.g. `fixedValue uniform 0.05`:

- HiSA Mach-3 smoke case, alphat 0: delivered ~24% of target, TFace
  pinned at the 20 K floor;
- same case, alphat 0.05: qError mean 1.7 W/m^2 on ~9.5 kW/m^2 targets
  (0.02%), TFace 119-177 K, no clamping.

Full analysis (including why the drifted face state cannot corrupt the
interior solution, and the list of contaminated wall-patch fields):
`../LiWallFunctions/FINDINGS_LiNeumann_TFace_limit`, addendum.

## Consumers of wall-patch fields

Same contamination rules as the Li Neumann family: state-type wall-patch
fields (T, rho, mu, psi, nut/alphat splits) are bookkeeping under this
coupling; flux-type outputs (`wallHeatFlux`, `wallShearStress`) are
exact by construction. Reconstruct wall states from dictionary `Tw` +
patch p, or read the diagnostics (`uTau`, `tauw`, `yPlus`, `qw`, `Tr`,
`delta`, `Ue`, `Te`, ...). The `yPlus()` of the nut BC is computed from
the physical wall state.

## HiSA integration notes (verified in the smoke test)

- Energy residual `fvc::laplacian(alphaEff, eCalc)` with the fresh
  `eCalc = rhoE/rho - |U|^2/2`: the prescribed gradient carries through
  because fixedGradient face values are internal + g/deltaCoeffs, so
  eCalc's boundary secant equals the prescribed gradient.
- Viscous Jacobian (`jacobianViscous/laplacian.C`) linearizes via
  `T.gradientInternalCoeffs()` = 0 for fixedGradient — consistent (same
  branch as adiabatic walls).
- The he-gradient mirroring in the T BC makes the `wallHeatFlux`
  function object report the imposed flux exactly (verified to all
  written digits through the `compact` energy laplacian).

## Cost

The inner solve is a shooting loop around a 48-node ODE integration per
wall face, re-evaluated by each of the three BCs every outer iteration
(the edge reconstruction is cached and shared; `edgeUpdateInterval`
throttles it). Measured impact on the smoke cases is negligible; on
large wall counts reduce `integrationNodes`/`maxIter` or raise
`edgeUpdateInterval` first.

## Validation status

Only engineering smoke tests so far (`tests/`): unit tests of every
model equation (incompressible-limit tau_w recovery, T-U endpoints,
adiabatic T_w = T_r, q_w sign and formula, edge finder analytics) and
integration tests (rhoSimpleFoam isothermal/adiabatic/prescribed/
parallel + HiSA Mach-3 flat plate with exact flux delivery). **No
physics validation against DNS/wall-resolved references has been
performed yet** (deliberately deferred); the natural campaign is the
Gao Ma=3 and M6Tw076 flat plates mirroring the Li acceptance protocol,
plus an a-priori test feeding DNS profiles through `solveInner`.
