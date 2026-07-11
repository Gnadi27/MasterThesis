# Wang Wall Functions - Design And Deviations Log

## Scope And Packaging

OpenFOAM v2412 user library implementing the compressible Analytical Wall
Function **AWF^Cmp** of

> X. Wang, J. Chen, T. Craft, H. Iacovides, "Extension of Analytical Wall
> Functions to Supersonic and Hypersonic Flows", *Flow, Turbulence and
> Combustion* 111:1-34 (2023),

with the wall heat flux imposed through a Neumann (prescribed-gradient)
temperature boundary condition, following the delivery mechanism of the
`GaoNeumann*` family in `../GaoWallFunctions`
(`DESIGN_GaoNeumann_energy_coupling.md` there documents the coupling
mechanics; they are identical here and only summarised below).

Equation labels: **WFxx** = paper equation numbers; **Txx.yy** = Wang's
PhD thesis (Manchester, 2018); **Cx.y** = Craft, Gerasimov, Iacovides,
Launder, *Int. J. Heat Fluid Flow* 23 (2002) 148-160.

The implemented model is the paper's final AWF^Cmp = thesis
**hyper-CMAWF**:

- MAWF: convection terms scaled quadratically with wall distance
  (WF17-WF20),
- CMAWF: wall-normal viscous dissipation of the *total*-energy equation
  retained in the analytical temperature equation (WF22-WF23) and, via
  the fvOption, in the discretized near-wall energy equation (WF24),
- hyperbolic molecular-viscosity variation across the zero-turbulent-
  viscosity sublayer with fluid properties evaluated at the sublayer
  edge (WF25-WF28),
- full appendix solution WF30-WF45, plus the thesis Appendix E thin-cell
  branch for y*_n < y*_v.

Packaging:

- Source directory: `$WM_PROJECT_USER_DIR/src/WangWallFunction`
- Library: `$FOAM_USER_LIBBIN/libWangWallFunctions.so`
- Build: `source /opt/of-v2412/etc/bashrc && wmake libso`
- No OpenFOAM core library or HiSA source is edited.
- Verification artifacts: `tests/TestWangWallFunction` (wmake app) and
  `tests/reference_wang_awf.py` (independent Python reference; run the
  app, pipe its output to the script).

Boundary condition / model type names:

| field    | type |
|----------|------|
| `T`      | `WangNeumannWallHeatFluxTemperature` |
| `nut`    | `nutWangNeumannWallFunction` |
| `epsilon`| `epsilonWangNeumannWallFunction` (faithful k-epsilon variant) |
| `omega`  | `omegaWangNeumannWallFunction` (**k-omega translation, deviation**) |
| `alphat` | standard `compressible::alphatWallFunction` (no Wang alphat class; see coupling notes) |
| `k`      | standard `kqRWallFunction` (zero-gradient) |
| fvOption | `WangNeumannDissipationSource` (CMAWF-D cell term) |

Shared evaluator: namespace `WangNeumannWallFunction`
(`WangNeumannWallFunction/`), `Coeffs` + `PatchData` + `calcPatch`. All
conditions and the fvOption on a patch must carry identical coefficients;
each evaluates `calcPatch` independently (same redundancy accepted in the
Gao library).

## Model summary as implemented

Non-dimensional wall distance (WF26): `y* = y rho_v sqrt(kP) / mu_v`,
with `kP` the near-wall-cell turbulence energy and the `v` state evaluated
at the sublayer-edge temperature `T_v` (see fixed-point below). Sublayer
edge `y*_v = 10.8`; turbulent viscosity `mu_t = mu_v alpha (y* - y*_v)`
beyond it, zero below (WF27, alpha = c_mu c_l = 0.2295); hyperbolic
molecular viscosity below it (WF28):
`mu = mu_v / (1 + b_mu (y* - y*_v))`, `b_mu = (mu_w - mu_v)/(mu_w y*_v)`.

Simplified near-wall equations solved analytically over the cell (with
first integrals `mu dU1/dy* = D y*^3/3 + C y* + A1` and
`mu_v Y dU2/dy* = D y*^3/3 + C y* + A2`, `Y = 1 + alpha(y*-y*_v)`):

- momentum (WF17): RHS `D y*^2 + C`,
- total energy (WF22): RHS `Dth y*^2 - (1/Cp) d/dy*(mu U dU/dy*)`.

Outputs per wall face:

- wall shear stress from `A1` (WF32/E.2; sign convention below),
- wall heat flux `q_w = -(rho_v cp sqrt(kP)/mu_v) A_th1` with `A_th1`
  from WF44 (isothermal) or wall temperature from WF42 with `A_th1 = 0`
  (adiabatic),
- cell-averaged k production `Pk` (WF13) and dissipation `epsBar` (WF16),
  nodal `eps_P = kP^{3/2}/(c_l y_P)` (C4.9),
- cell-averaged wall-normal viscous dissipation `dissAvg` for the
  CMAWF-D term (WF24), evaluated exactly as `Un tau(y_n)/y_n` (see
  below).

The monolithic printed expressions N (WF35), B2 (WF34), BT (WF41), A_th1
(WF44), T_w (WF42) are not transcribed literally; they are assembled from
the profile/continuity conditions with `D1 = D2`, `C1 = C2`,
`Dth1 = Dth2` (see "Single nodal coefficients" below), which reproduces
the printed forms exactly in that case. This was verified term-by-term by
hand against the paper appendix, including the checks that the
implemented U1/U2 differentiate back to the exact first integrals with
the hyperbolic mu(y*).

## Constants and dictionary entries

All BCs and the fvOption accept the same entries (defaults = paper
values):

```foam
yStarV              10.8;        // WF: sublayer edge (see note 12)
cl                  2.55;        // c_l; alpha = Cmu*cl, y*_d = 2*cl
Cmu                 0.09;
Prt                 0.9;
nIntegrationPoints  40;          // Simpson intervals (numerics)
TvIterations        3;           // Tv/adiabatic-Tw fixed-point iterations
TvRelaxation        0.5;         // under-relaxation of that fixed point
viscosityVariation  hyperbolic;  // hyperbolic (paper) | constant (b_mu=0)
thermalMode         auto;        // auto | isothermal | adiabatic
Tw                  -1;          // physical wall T, required isothermal
wallTransport
{
    mode    sutherland;          // sutherland | constant | fromPatch
    As      1.458e-06;
    Ts      110.4;
    // Pr   0.7;                 // paper value; omitted -> modified Eucken
    // mode constant: mu <val>; kappa <val>;
}
```

`viscosityVariation constant` sets `b_mu = 0` and recovers the
constant-property CMAWF limit at runtime. The paper used constant
Pr = 0.7 for air; supply `wallTransport.Pr` to reproduce that, otherwise
kappa is estimated by modified Eucken (`mu Cv (1.32 + 1.77 R/Cv)`) as in
the Gao library.

The T condition additionally accepts the Neumann-coupling guards and
diagnostics (identical semantics to GaoNeumann):
`gradientRelaxation` (default 1), `TFaceMin`/`TFaceMax` (1/GREAT),
`diagnostics`, `diagnosticLog`, `diagnosticWrite`,
`diagnosticWriteInterval`, `diagnosticWritePrecision`. Face-wise output
goes to `postProcessing/WangNeumannWallFunctionDiagnostics/<patch>/
<time>.dat` (per-processor subdirectories in parallel).

## The Neumann energy coupling (same as GaoNeumann)

Each outer iteration the T condition computes q_w from `calcPatch` and
places the face enthalpy so the discrete wall diffusion flux equals the
target exactly:

```text
qTargetOF   = +q_w                       // conventions coincide, see below
he_f,needed = he_1 + qTargetOF/(alphaEff*deltaCoeffs)
T_f,needed  = THE(he_f,needed, p_f)      // thermo's own inverse
gradient()  = (T_f,needed - T_1)*deltaCoeffs
```

`alphaEff` cancels between inversion and flux; there is no sign
constraint, so heat flux into the wall with `T_1 < T_w` is representable
(the reason this family exists; see
`../LiWallFunctions/FINDINGS_M6Tw076_energy_coupling.md`). The delivered
he-gradient is mirrored onto the paired `gradientEnergy` he patch so
`wallHeatFlux` reports the imposed flux. Adiabatic mode prescribes a zero
gradient (identical to `zeroGradient`) while the closure uses the WF42
recovery wall temperature. Under this coupling the **patch T value is the
tangent-extrapolated face temperature, not the physical wall
temperature**; hence the dictionary `Tw` and `wallTransport` requirements
(refused otherwise), exactly as in the Gao family.

## Sign conventions

- **q_w**: the paper's q_w (WF43, `-lambda dT/dy` at the wall) is
  positive when the wall heats the fluid. OpenFOAM's discrete wall
  diffusion flux `alphaEff*snGrad(he)` (also what the `wallHeatFlux`
  function object reports) is positive in exactly the same situation,
  so the imposed target is **`+q_w` — no sign flip**. This is the
  opposite mapping to the Gao library, whose `-q_w` reflects the
  Nichols/Gao convention (q_w positive from the fluid into the wall);
  copying Gao's mapping here inverted the flux and was caught by the
  first restart smoke test (Wang wallHeatFlux mirrored Gao's converged
  values with the wrong sign). Note the paper itself plots Q with the
  opposite (into-the-wall positive) sign in its results sections.
- **tau_w**: the papers print `tau_wall = -(rho_v sqrt(kP)/mu_v) A1`
  (T8.12, C4.5) — the reaction force *on the fluid*, negative for
  attached flow (cf. Craft Fig. 6, tau_wall drawn in -x under the CV).
  OpenFOAM's wall-stress path reproduces the conventional
  `mu dU/dy|_0 = +(rho_v sqrt(kP)/mu_v) A1`, positive for attached flow.
  `PatchData::tauw` stores the **conventional** sign; this is a pure
  sign-bookkeeping translation, not a model change. (Craft's own Eq. 2.1
  log-law tau_wall is positive for positive U_P, so the two conventions
  coexist within the papers as well.)

## Deviations and mapping choices

Numbered for reference; 1-4 are faithful-intent FV mappings, 5-9 are
genuine deviations or extensions, 10-13 are typo adjudications between
paper and thesis, 14+ are numerics.

1. **Unstructured-FV frame.** The paper works in structured wall-aligned
   (x, y). Here, per wall face: `nHat` = inward wall normal, `tHat` =
   unit projection of the near-wall cell velocity onto the wall plane.
   Nodal terms:
   - momentum convection (non-conservative nodal form of T8.105):
     `conv = rho_P (U_t d(u.tHat)/dx_t + U_n d(u.tHat)/dx_n)` from a
     local Gauss cell gradient of U;
   - pressure gradient: `C ~ tHat . grad(p)` (local Gauss);
   - energy convection: `Dth ~ div((rho E + p) u)` (local Gauss cell
     divergence), which equals the paper's
     `d(rho U E + U p)/dx + d(rho V E + V p)/dy` for 2-D flows and is its
     coordinate-free generalisation. `E = e + |u|^2/2`, with
     `e = he - p/rho` if the thermo energy variable is an enthalpy
     (detected from the field name "h..." vs "e...").
   The local Gauss sums use linear face interpolation with the mesh's
   own weights; wall-face values enter with their (possibly drifted)
   patch values. Faces on processor boundaries use the evaluated coupled
   patch values.

2. **North-face reconstruction (Un, Tn, y_n).** The paper's `U_n`, `T_n`,
   `y_n` live on the face of the near-wall cell opposite the wall (Craft
   Fig. 4). On a general mesh this face is identified as the cell face
   whose outward normal is most aligned with `nHat` (alignment > 0.5
   required); `y_n` is the wall-normal distance of its centre, `Un/Tn`
   are linearly interpolated face values (`Un` relative to the wall
   velocity, projected on `tHat`). Fallback for degenerate cells:
   `y_n = 2 y_P`, `Un = |u_P - u_wall|` along tHat, `Tn = T_P`.

3. **Single nodal coefficients.** The paper evaluates the convection
   terms once from nodal values (WF18: "D represents the convection
   terms evaluated from the nodal values") and scales both regions by
   `1/y*_n^2`; hence `D1 = D2 = D`, `C1 = C2 = C`, `Dth1 = Dth2 = Dth`
   here. The `(D1-D2)`, `(C1-C2)`, `(Dth1/3 - Dth2/3)` groups in
   WF33/WF35/WF40/WF42/WF44 then vanish identically.

4. **CMAWF-D dissAvg is evaluated exactly, not by quadrature.** WF24's
   cell average telescopes:
   `(1/y_n) [ (mu+mu_t) U dU/dy ]_0^{y_n} = Un tau(y_n)/y_n` since
   `U(0) = 0` and the first integral of the momentum equation gives
   `(mu+mu_t) dU/dy` in closed form at `y_n`. The paper integrates the
   expanded two-term form numerically; the telescoped form is the same
   number to quadrature accuracy.

5. **CMAWF-D delivered as an fvOption correction difference.** Wang
   modified the solver to *replace* the discretized near-wall viscous
   dissipation by dissAvg. In HiSA (and any unmodified OpenFOAM solver)
   that term is embedded in the face-based `div(sigma . u)` flux
   assembly and cannot be removed, so
   `fv::WangNeumannDissipationSource` adds the **difference**
   `S = dissAvg - dissLocal` [W/m^3] in near-wall cells, with
   `dissLocal = mu_eff Un (d(u.tHat)/dn) / y_n` a discrete estimate of
   the same wall-normal dissipation from the resolved cell gradient
   (wall-face contribution exactly zero because U = 0 there). **This is
   the roughest paper-to-OpenFOAM mapping in the library**: dissLocal is
   a cell-centred reconstruction, not the solver's literal face-flux
   contribution. Omit the fvOption to run plain CMAWF (dissipation only
   inside the analytical temperature profile). HiSA picks the source up
   through `fvOptions(rho, thermo.he())` in its rhoE residual
   (`modules/hysa/src/modules/hisa/residualsUpdate.H`, `#else` branch
   for v2412).

6. **epsilon variant: cell-averaged sink through the G field.** The
   paper solves the k equation in the near-wall cell with cell-averaged
   `Pk` and `epsBar` and prescribes nodal
   `eps_P = kP^{3/2}/(c_l y_P)`. OpenFOAM's kEpsilon sinks the *nodal*
   epsilon implicitly (`fvm::Sp(rho eps/k, k)` structure), so
   `epsilonWangNeumannWallFunction` (derived from `epsilonWallFunction`
   for its master-patch/corner machinery and matrix constraint) sets
   `epsilon0 = eps_P` and `G = Pk/rho + eps_P - epsBar` (kinematic,
   OpenFOAM G convention). At steady state the net k source is
   `Pk/rho - epsBar` as in the paper; the transient/implicit
   linearisation differs from Wang's in-house discretization. No Craft
   F_epsilon laminarisation factor (C3.14-C3.16) is applied — Wang omits
   it.

7. **omega variant is an untested translation (deviation).** The papers
   never couple the AWF to k-omega. `omegaWangNeumannWallFunction` sets
   `omega_P = eps_P/(Cmu kP)` and the same G combination; with the SST
   k-sink `beta* omega k`, `beta* = Cmu = 0.09`, the implicit near-wall
   sink equals eps_P, so the same steady-state identity holds. No
   Nichols-style inner/outer blending. Flagged clearly: results with
   this variant are not covered by any validation in the papers.

8. **k wall treatment.** `kqRWallFunction` (zero-gradient), consistent
   with Craft's neglect of molecular diffusion of k to the wall. No
   custom k class.

9. **Tv fixed point.** The papers state that `T_v` comes from the
   thermal AWF and `rho_v`, `mu_v` from `T_v` (WF25-WF26 discussion) but
   not the update scheme. Implemented: `TvIterations` (default 3)
   fixed-point passes per calcPatch call, under-relaxed by
   `TvRelaxation` (default 0.5), initialised from the cell temperature
   each call (no persistent state); the adiabatic wall temperature (WF42
   with A_th1 = 0) participates in the same loop. `T_v` is evaluated as
   `T1(y*_v)` (continuity value). The stored outputs of the final pass
   use the pre-update Tv/Tw of that pass; with the relaxed fixed point
   the difference is a higher-order residual and vanishes as the outer
   iterations converge.

10. **Region-1 coefficient scaling: paper over thesis.** Thesis
    T8.78/T8.90/T8.105/T8.117 print `1/y*_v^2` scalings for the
    region-1 (sublayer) coefficients while the paper (WF18, WF23) scales
    both regions by `1/y*_n^2`. The `1/y*_v^2` variant makes the assumed
    RHS `D y*^2` discontinuous at `y*_v` for equal nodal convection
    values; the peer-reviewed paper's single `1/y*_n^2` is used.
    Adjudicated as a thesis typo.

11. **Property state in D, C, Dth: thesis over paper.** Paper WF18/WF23
    still print `mu_w, rho_w` (they predate the variable-property
    Sect. 3.5), while the thesis hyper-CMAWF (T8.105 ff.) uses
    `mu_v, rho_v` consistently with the WF26 definition of y*. The `v`
    state is used, matching the final variable-property model.

12. **y*_v = 10.8 (paper/Craft) vs 10.1 (thesis Sect. 6.2.3).**
    Dictionary-adjustable; default 10.8. Adjudicated as a thesis typo
    (Craft 2002 and the paper both say 10.8 from pipe-flow numerical
    experiments).

13. **Sign of WF13 (Pk).** The paper prints a leading minus; production
    of k must be positive (`Pk = + (1/y_n) int mu_t (dU/dy)^2 dy`, cf.
    C4.6 and thesis T8.13 without the minus). Implemented positive.
    Similarly, thesis Appendix E carries `C_th` terms absent from the
    paper's final closed system (the pressure-work convection is inside
    Dth per WF23); the paper's form is implemented and the thesis E.3/E.4
    are used with those terms folded into Dth (i.e. dropped as separate
    coefficients).

14. **Thin-cell branch (y*_n < y*_v), thesis Appendix E.** `mu_t = 0`
    across the whole cell, `b_mu = (mu_w - mu_n)/(mu_w y*_n)` (E.1; the
    `v` state relaxes to the `n` state through the Tv fixed point),
    momentum constant from E.2, thermal from E.3 (adiabatic) / E.4
    (isothermal) with the numerically integrated dissipation integral
    over [0, y*_n]; `Pk = 0`; `epsBar = 2 nu_v kP / y_d^2` when
    `y*_n <= y*_d` (C4.8), where `y_d` is the physical position of
    `y*_d = 2 c_l = 5.1`.

15. **nut inversion (as Gao).** `nut = max(0, (tau_w - mu_patch
    |snGrad(U)|)/(rho_patch |snGrad(U)|))` with **patch-field** mu and
    rho — the values OpenFOAM's stress operator actually uses, which
    differ from the physical wall state under the Neumann coupling
    because the patch T drifts. Negative tau_w (reversed near-wall flow
    relative to the discrete gradient) and zero-gradient faces fall back
    to `nut = 0`; the wall then carries at least the laminar stress.
    `yPlus()` on the nut condition reports the AWF y* (WF26), not a
    log-law y+.

16. **Numerical guards (not physics).** `kP` floored at SMALL (negative
    powers of k would overflow from an uninitialised k = 0 — initialise
    k > 0); `ln` arguments floored at SMALL; `y*_n` floored at SMALL;
    temperatures kept positive; the A1/A_th1 denominators are provably
    positive (`b_mu < 1/y*_v` implies `denominator > y*_v/2 > 0`) but
    are guarded anyway; Simpson quadrature with a fixed even interval
    count (`nIntegrationPoints`, default 40) for the Pk integral and the
    two thermal dissipation integrals I1, I2 (the integrands are smooth
    closed-form expressions built from the analytical profiles and first
    integrals — no nested numerical differentiation).

17. **Evaluation cost/consistency.** Every consumer (T, nut,
    epsilon/omega, fvOption) calls `calcPatch` independently per outer
    iteration, so the closure is evaluated up to 4x per patch per
    iteration with identical inputs except fields updated in between
    (same accepted redundancy and same caveat as the Gao library: all
    conditions must carry identical coefficient entries; there is no
    cross-check).

18. **fromPatch transport mode.** Only valid for adiabatic walls (where
    the floating patch T is physical). Even there, `mu_v` is then the
    patch (wall-temperature) viscosity rather than the T_v viscosity —
    an additional approximation; prefer `sutherland` everywhere.

## Wall alphat conditioning (important)

The T condition inverts the target flux through
`alphaEff*deltaCoeffs` to place the face enthalpy, and the solver
delivers `alphaEff*snGrad(he)` with the same `alphaEff` — the value of
the wall `alphat` therefore **cancels exactly** from the delivered heat
flux. It does NOT cancel from the *face temperature*: the extrapolated
face temperature is `T_1 + qTarget/(alphaEff*deltaCoeffs*cp-ish)`, so a
small `alphaEff` (e.g. `alphat = 0` with a large flux on a coarse wall
cell) demands an absurd face temperature, the `TFaceMin/Max` clamps
bind, the delivered flux undershoots, and — because the drifted face
temperature feeds the solver's face-state (rho, mu) evaluations — the
coupled system can destabilise.

Quantified on the Mach-3 flat-plate smoke case (y* ~ 1100, q_w ~
-0.6 MW/m2): with `alphat = 0` the required face temperature is
~-2400 K, the clamp pins TFace at its floor, and only ~0.5% of the flux
is delivered; with `compressible::alphatWallFunction` (alphat =
rho nut/Prt evaluated at the drifted patch state) ~60-85%; with a
prescribed `fixedValue alphat = 0.05` the face temperature is physical
(270-380 K), no clamp binds, and delivery is exact (qError ~ 3 W/m2 on
a 572 kW/m2 target). Recommendation: give the wall `alphat` a fixed
value large enough that `|q_w|/(alphaEff*deltaCoeffs)` is a moderate
enthalpy increment; this has no effect on the imposed physics. (This
deviates from the Gao/Li Neumann write-ups, which pinned alphat to 0 —
their regimes evidently kept the required face temperatures
representable.)

## Solver coupling notes (HiSA)

- Verified in the local HiSA fork
  (`~/OpenFOAM/mgnadlinger-v2412/modules/hysa`):
  `src/modules/hisa/residualsUpdate.H` builds
  `rhoESource = (fvOptions(rho, thermo.he()) & thermo.he()) + ...`, so
  the dissipation fvOption enters the rhoE residual. The v2412 build
  takes the `#else` (`fv::options`) branch.
- The T/nut/turbulence conditions act through standard patch-field
  evaluation and need no solver hooks; this matches how the GaoNeumann
  family is used with HiSA.
- A prescribed-gradient T condition enters the energy matrix as an
  explicit flux with no diagonal feedback: crude startup fields may need
  `gradientRelaxation`, the `TFace` clamps, or the solver's
  pseudo-transient damping (same guidance as GaoNeumann; the gradient
  starts from zero = adiabatic when absent from the initial field file).

## Example boundary conditions

```foam
// system/controlDict (or fvOptions-loading solver dict)
libs ("libWangWallFunctions.so");

// 0/T
wall
{
    type        WangNeumannWallHeatFluxTemperature;
    thermalMode isothermal;          // or adiabatic
    Tw          300;
    wallTransport
    {
        mode    sutherland;
        As      1.458e-06;
        Ts      110.4;
        Pr      0.7;                 // paper value for air
    }
    // gradientRelaxation 0.3;       // startup aids (optional)
    // TFaceMin           50;
    // TFaceMax           5000;
    diagnostics false;
    value       uniform 300;         // initial face value only
}

// 0/alphat -- wall alphat is pure conditioning here: alphaEff cancels
// between the T-gradient inversion and the delivered flux, so ANY value
// gives the same wall heat flux. A LARGER alphat keeps the extrapolated
// face temperature physical (see "Wall alphat conditioning" below).
wall { type fixedValue; value uniform 0.05; }
// or: wall { type compressible::alphatWallFunction; Prt 0.9; value uniform 0; }

// 0/nut
wall
{
    type        nutWangNeumannWallFunction;
    thermalMode isothermal;
    Tw          300;
    wallTransport {mode sutherland; As 1.458e-06; Ts 110.4; Pr 0.7;}
    value       uniform 0;
}

// 0/k
wall { type kqRWallFunction; value uniform 1e-8; }

// 0/epsilon  (kEpsilon runs)
wall
{
    type        epsilonWangNeumannWallFunction;
    thermalMode isothermal;
    Tw          300;
    wallTransport {mode sutherland; As 1.458e-06; Ts 110.4; Pr 0.7;}
    value       uniform 1;
}

// 0/omega  (kOmegaSST runs; translation variant, see deviation 7)
wall
{
    type        omegaWangNeumannWallFunction;
    thermalMode isothermal;
    Tw          300;
    wallTransport {mode sutherland; As 1.458e-06; Ts 110.4; Pr 0.7;}
    value       uniform 1;
}

// constant/fvOptions  (CMAWF-D term; omit for plain CMAWF)
WangDissipation
{
    type        WangNeumannDissipationSource;
    patches     (wall);
    thermalMode isothermal;
    Tw          300;
    wallTransport {mode sutherland; As 1.458e-06; Ts 110.4; Pr 0.7;}
}
```

## Status

- `wmake libso` clean on the first pass; library at
  `$FOAM_USER_LIBBIN/libWangWallFunctions.so`.
- Hand verification of the appendix algebra done (U1/U2 differentiate to
  the exact first integrals with the hyperbolic mu(y*); N/B2/BT/A_th1/T_w
  reduce to the printed WF forms with D1=D2, C1=C2, Dth1=Dth2; profile
  continuity at y*_v holds by construction).
- Numerical cross-check passed (`tests/`): `TestWangWallFunction`
  (includes the library translation unit and drives the per-face solver
  on 8 synthetic cases — heating/cooling/adiabatic, adverse/favourable
  gradients, thick and thin cells) against
  `tests/reference_wang_awf.py`, an independent Python transcription of
  the printed paper equations. 120 values compared, 0 mismatches at
  1e-9 relative tolerance; the reference's own internal checks
  (U/stress/T continuity at y*_v, U2(y*_n)=Un, T2(y*_n)=Tn, exact
  reduction to the Craft 2002 forced-convection constants for
  b_mu = D = C = Dth = 0, and the q_w <-> (Tn - Tw) inversion round
  trip) all pass.
- HiSA smoke tests (Mach-3 flat plate, `~/Gao_TurbZPGFP/1e-3` mesh and
  freestream, restarted from the converged Gao solution; hysa, serial):
  - **kOmegaSST + full Wang BC set** (T/nut/omega + kqRWallFunction +
    fixed wall alphat 0.05): stable over 500 continued iterations, no
    NaN; tau_w 910-1620 Pa; imposed q_w mean -572 kW/m2 with qError
    mean ~3 W/m2 (exact delivery); `wallHeatFlux` reports the imposed
    flux (he-gradient mirroring works). The converged Gao closure on
    the same case gives mean -591 kW/m2 — 3% apart.
  - **kEpsilon + epsilonWangNeumannWallFunction** (epsilon field built
    as Cmu k omega from the converged fields): stable 100 iterations,
    qError mean ~93 W/m2, wall quantities consistent with the SST
    variant.
  - **WangNeumannDissipationSource**: stable when enabled on a healthy
    coupled state (qError 0.06% after 100 further iterations). Enabling
    it on a badly conditioned state (TFace clamp pinned, see the alphat
    section) drove a divergence — enable the fvOption only once the
    Neumann T coupling delivers its target (qError ~ 0 in the
    diagnostics).
  - Restart round-trip of the written T gradient exercised by the
    100 -> 600 continuation runs.
  - Two implementation bugs were found and fixed only by these tests:
    the q_w sign mapping (see "Sign conventions") and out-of-bounds
    boundary-field access on empty patches in the local Gauss sums.
