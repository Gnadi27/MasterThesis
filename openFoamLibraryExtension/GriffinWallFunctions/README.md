# GriffinWallFunctions

RANS (Spalart-Allmaras) wall-function library for OpenFOAM v2412 + HiSA
implementing, as faithfully as possible, the near-wall model of

> K. P. Griffin, L. Fu, P. Moin, "Near-wall model for compressible
> turbulent boundary layers based on an inverse velocity transformation",
> CTR Annual Research Briefs 2022 ("the paper"),

built on the total-stress-based velocity transformation (Griffin, Fu,
Moin, PNAS 118, 2021) and the boundary-layer-edge local-reconstruction
method (Griffin, Fu, Moin, PRF 6, 2021; compressible form of Appendix D).
The reference PDFs live in `../GriffinWallFunction/`. Deviations from the
papers are catalogued in `DESIGN_AND_DEVIATIONS.md` — read it before
trusting results.

## What the model does

Per wall face and outer iteration, the shared evaluator
(`GriffinWallFunction::calcPatch`) solves the paper's inner problem:

1. Incompressible constant-stress law of the wall with Van Driest
   damping, `S_t+ = 1/(1 + kappa y* (1 - exp(-y*/A+))^2)`, with
   Reynolds-number-dependent `kappa(Re_tau*)`, `A+(Re_tau*)` linearly
   interpolated from the paper's channel-DNS calibration table
   (constant extrapolation outside Re_tau* in [182, 5190]).
2. Inverse total-stress transformation (paper Eq. 2.5) integrated in
   semi-local coordinates from the wall to the matching point (the
   first-cell centre), with the algebraic temperature-velocity relation
   (paper Eq. 2.8, s = 1.14, r = Pr^(1/3)) supplying T(U), rho(T), mu(T)
   at every node — no temperature ODE.
3. Shooting (bracketed bisection) on tau_w such that U(y_1) matches the
   first-cell velocity.
4. Wall heat flux algebraically: `q_w = s tau_w Cp (T_w - T_r)/U_e`,
   positive for heat entering the fluid (Griffin's B_q convention =
   OpenFOAM's outward-normal boundary-flux sign). Adiabatic mode gives
   q_w = 0 and T_w = T_r in closed form.

Edge quantities (U_e, T_e, delta) come from the PRF-2021
local-reconstruction method: stagnation enthalpy along a wall-normal cell
walk, `U_I^2 = 2(h_0,ref - h) - V^2`, edge at the first wall-outward
crossing of U/U_I = 0.99. Cached per patch, recomputed every
`edgeUpdateInterval` time indices, shared by all three boundary
conditions. A `prescribed` mode overrides with dictionary values.

## Boundary conditions

All three must carry **identical** wall-function coefficients (they each
re-evaluate `calcPatch`; diverging coefficients silently decouple the
momentum, turbulence, and energy closures):

| field | type | role |
|---|---|---|
| `nut` | `nutGriffinWallFunction` | delivers tau_w exactly through the discrete stress `(mu + rho nut) snGrad(U)` |
| `nuTilda` | `nuTildaGriffinWallFunction` | wall value 0; pins the first-cell nuTilda to match the model eddy viscosity (matrix manipulation) |
| `T` | `GriffinWallHeatFluxTemperature` | fixedGradient (Neumann) BC imposing q_w directly; pairs with `gradientEnergy` on he |
| `alphat` | `fixedValue` (a conditioning constant, see below) | not a wall function — energy coupling is Neumann-only |

Example (isothermal wall; full snippets in
`examples/boundaryConditions.foam`):

```
wall
{
    type            GriffinWallHeatFluxTemperature;
    thermalMode     isothermal;     // auto | isothermal | adiabatic
    Tw              300;            // physical wall temperature [K]
    viscosity
    {
        mode        sutherland;     // sutherland | powerLaw | constant
        As          1.458e-06;
        Ts          110.4;
    }
    gradientRelaxation 0.3;
    TFaceMin        20;
    diagnostics     true;
    value           uniform 300;
}
```

## Two warnings you must not skip

1. **The wall patch T is NOT the wall temperature.** Under the Neumann
   coupling the patch T is the tangent-extrapolated face temperature and
   drifts far below the physical `Tw` (that is the point: it makes the
   T_1 < T_w regime representable). Everything evaluated from wall-patch
   thermo fields (patch T, rho, mu, psi, alphat, nut splits) is a
   bookkeeping state. The physical wall state lives in the BC dictionary
   (`Tw`, `viscosity`) and in the diagnostics output. Flux-type outputs
   (`wallHeatFlux`, `wallShearStress`) are correct by construction; the
   `wallHeatFlux` function object reports exactly the imposed q_w. See
   `../LiWallFunctions/FINDINGS_LiNeumann_TFace_limit` for the full
   analysis.
2. **Set wall `alphat` to a conditioning constant, not 0.** alphaEff
   cancels identically from the delivered flux, but it scales the face
   temperature excursion `TFace - T_1 ~ q / (alphaEff deltaCoeffs cp)`.
   With laminar-only alphaEff (alphat 0) at low density/coarse spacing
   the required TFace becomes unphysically cold and hits the `TFaceMin`
   clamp, and the delivered flux silently collapses below target.
   `fixedValue uniform 0.05` keeps TFace within a few K of T_1 (verified
   in the smoke tests: qError drops from ~75% to ~0.02%).

## Dictionary reference (shared `Coeffs`)

| entry | default | meaning |
|---|---|---|
| `s` | 1.14 | Reynolds-analogy factor of Eq. 2.8 |
| `Pr` | 0.72 | molecular Prandtl number (r = Pr^(1/3)) |
| `cv1` | 7.1 | SA fv1 constant for the nuTilda inversion |
| `maxIter` | 60 | shooting bisection iterations |
| `tolerance` | 1e-6 | relative shooting/sub-iteration tolerance |
| `integrationNodes` | 48 | log-spaced ODE intervals wall -> first cell |
| `subIterations` | 5 | per-node property sub-iterations |
| `thermalMode` | auto | auto / isothermal / adiabatic |
| `strict` | false | fatal error instead of laminar fallback |
| `Tw` | -1 | physical wall temperature (required for isothermal) |
| `viscosity` | air Sutherland | mu(T) law: sutherland (As, Ts), powerLaw (muRef, TRef, n), constant (mu) |
| `edgeMode` | reconstruct | reconstruct / prescribed (then Ue, Te, delta required) |
| `maxEdgeCells` | 200 | cell-walk length limit |
| `maxEdgeDistance` | GREAT | cell-walk distance limit |
| `edgeUpdateInterval` | 1 | recompute edge every N time indices |
| `edgeThreshold` | 0.99 | U/U_I defining the edge |

T-BC-only entries: `gradientRelaxation` (0,1], `TFaceMin`/`TFaceMax`
(face-temperature clamp; a *binding* clamp means the flux is NOT being
delivered — fix conditioning/mesh instead of widening silently),
`diagnostics*` switches (per-face `.dat` under
`postProcessing/GriffinWallFunctionDiagnostics/<patch>/`).

## Status codes (diagnostics column 5)

0 OK; 1 near-zero first-cell velocity (trivial state); 2 invalid edge
data (Ue <= U_1, bad Te, or delta <= y_1) — laminar fallback; 3 shooting
did not converge — laminar fallback. Edge status (column 6): 0 crossing
found, 1 no crossing (outermost sample used — enlarge maxEdgeCells /
check decomposition), 2 prescribed.

## Parallel limitation

The edge cell walk cannot cross processor boundaries. Decompose
wall-tangentially (e.g. `simple` with n = (Nx 1 Nz) over a wall in the
y-normal); otherwise affected faces fall back to the outermost reached
sample (warned once) — or use `edgeMode prescribed`.

## Build & test

```
. /opt/of-v2412/etc/bashrc
wmake libso                      # -> $FOAM_USER_LIBBIN/libGriffinWallFunctions.so
tests/equations/Allrun           # unit tests (31 assertions)
tests/integration/Allrun         # rhoSimpleFoam + HiSA smoke tests
```

Load in a case with `libs ("libGriffinWallFunctions.so");` in
controlDict. Works unchanged with HiSA (`hysa`): the smoke test verifies
exact flux delivery through HiSA's `laplacian(alphaEff,e) ... compact`
scheme and correct `wallHeatFlux` reporting.
