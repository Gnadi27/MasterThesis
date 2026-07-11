# Findings: Li Wall Function on M6Tw076 — Energy-Coupling Breakdown at Coarse y+

Status notes from the paper-consistency review and the M6Tw076 (Zhang/Duan/Choudhari
DNS database, Ma 5.86, Tw = 300 K, Tw/Tr = 0.76, -B_q = 0.02) validation
campaign. Case family: `~/Xue_et_al/flat_plate/M6Tw076/`, solver HiSA,
SA turbulence model, perfect gas (hConst/Sutherland).

## 1. Summary

- The implementation transcribes Li et al. Eqs. (32)-(42) faithfully. The
  implemented unified law (Eq. 42) reproduces the paper's explicit log-layer
  law (Eq. 30) at large y+ and u+ = y+ in the sublayer. No transcription
  errors found.
- A-priori, against DNS: velocity law matches below y+ = 10; the
  enthalpy-velocity relation (Eq. 12/40) matches the DNS T-u curve within
  ~2 %. Fed with the *true* DNS cell states, the implementation's
  Eq. (40)->(41)->(33) chain returns qw within ~25 % of DNS. **The laws are
  fine.**
- A-posteriori, at y+_1 ~ 9-10, both `thermalFluxTarget qHalf` (old default
  via `auto`) and `qw` produce large errors (Cf +15-20 %, Ch +50 % to +200 %,
  spatially erratic). Root cause is **not** the wall-function model and
  **not** FVM per se, but the OpenFOAM coupling idiom: imposing the wall heat
  flux by inverting a non-negative turbulent diffusivity `alphat`.

## 2. The core issue: alphat inversion vs. sublayer temperature peak

**This is the headline finding.**

> Gao et al Ma=3_Tw=300K flat plate couldn't expose this, because up until edge (288K), temperature within the BL > Tw
> In real flight scenario might become the case, when hor surface (easily 400 - 600K or more) and highly compressed BL (=> peak and drop-after within the BL)

### Physics of the near-adiabatic case

For Tw/Tr = 0.76 the DNS temperature profile peaks *inside the viscous
sublayer*:

| location   | DNS value                       |
|------------|---------------------------------|
| T peak     | 300.4 K at y+ = 2.6             |
| T(y+ 5.2)  | 294.8 K                         |
| T(y+ 10.7) | 266.3 K  (**34 K below Tw!**)   |
| T(y+ 31.8) | 202.4 K                         |
| qw         | +2142 W/m^2 (*into* the wall)   |

Heat is generated around the peak and flows *both* ways: down into the wall
(qw > 0) and up into the boundary layer. Consequently the physically correct
coarse-mesh state for y+_1 > ~3 has **T_1 < T_w together with heat flowing
into the wall**. The paper's own Eq. (34) anticipates this: the total heat
flux q(y) = qw - r*u*tau_w/Pr_t changes sign inside the sublayer.

### Why the OpenFOAM idiom cannot represent this state

OpenFOAM-style thermal wall functions deliver the wall flux as

```
q_OF = alphaEff * snGrad(he),   alphaEff = CpByCpv*(alpha_lam + alphat),  alphat >= 0
```

i.e. the flux must be expressible as a **positive diffusivity times the
resolved cell-to-wall gradient**. A positive diffusivity can only move heat
from hot to cold. With T_1 < T_w the discrete gradient points *out of* the
wall while the physical flux points *into* it: the required flux is outside
the representable set, for any target `qw`.

What the coupled system then does (confirmed by diagnostics in both runs):

1. The inversion clamps `alphat` to 0; the laminar leak *heats* the first
   cell; T_1 can never settle below T_w.
2. The solution is pushed to the constraint boundary **T_1 = T_w**:
   `h1 - hw` pinned within ~100 J/kg on 538/570 faces, `snGrad(he) ~ 0`,
   `alphat` exploding to 1e10-1e11 (spatially alternating with O(1) values).
3. Eq. (40), fed the degenerate input h1 = hw plus a cold second cell,
   infers a spurious h_rg - h_w ~ 150 kJ/kg (true value with physical
   states: ~41 kJ/kg) -> model qw of 5000-10000 W/m^2 vs. DNS 2142.
4. With `thermalFluxTarget qw` the full spurious amount is extracted
   (positive feedback, worse than `qHalf`, station-to-station oscillations
   in Ch including sign-reversing dips). The corrupted beta*, Gamma* also
   inflate u_tau -> Cf +15-20 % in both runs.

> expressing the flux as alphat*snGrad(he) turns "T1 > T_w" into a solvability requirement at equilibrium: a positive diffusivity can only carry heat into the wall if the discrete gradient points that way, so a steady state with the target flux delivered must have T1 > T_w.

### Mental model: tangent, secant, and the self-consistency loop

Two different slopes live at the wall face:

- the **tangent**: the true wall gradient, dT/dy|_w = qw/lambda_w. This is
  what the wall function *computes* from the first/second-cell states via
  the closures (Eq. 40 -> 41 -> 33).
- the **secant**: the discrete gradient OpenFOAM's boundary discretization
  sees, snGrad(he) = (he(T_1) - he(T_w))/Delta — the straight chord from
  the wall value to the current cell-centre value. On a coarse mesh this
  differs from the tangent *by design*; the mismatch is not an error.

**alphat is pure bookkeeping.** The energy equation can only transmit wall
heat in the form `flux = alphaEff * snGrad(he)`, so the BC solves that
single equation for the single free coefficient:

```
alphat = -qw / (CpByCpv * snGrad(he)) - alpha_laminar
```

Nothing is modeled in this step; alphat carries zero independent
information. It is not a local turbulence quantity — it is the conversion
factor between the flux the closure already decided and the format the
discretization demands (fixed total = unit price x quantity; the quantity
(secant) is dictated by the current fields, so the unit price is backed
out). Large values (1e2-1e6 x laminar) are normal and harmless; alphat
becomes a *symptom* only when the secant it compensates for is degenerate
(T_1 pinned at T_w). Corollary: nonlinearity of the true T(y) between wall
and cell centre does **not** corrupt the delivered flux — the "wrong"
secant slope is expected and absorbed exactly. (This differs from
classical Jayatilleke-type alphat wall functions, where alphat comes from
an assumed T+ law and profile-shape error maps directly into flux error.)

**There is no T_1 setpoint.** The wall function never targets a first-cell
temperature; the only thing it ever enforces is flux-through-secant. Per
outer iteration:

1. current T_1, T_2  ->  closure infers qw           (the tangent)
2. qw / current secant  ->  alphat                   (bookkeeping)
3. cell-1 energy balance  ->  new T_1
   (imposed wall flux + top-face turbulent flux via the prescribed
   mut_1/nuTilda_1 + advection + dissipation)
4. repeat.

The converged state is a **self-consistency fixed point**: T_1 settles at
the value where the qw inferred *from* T_1 equals the heat the surrounding
balance actually transports to the wall face. With exact closures and an
exact cell-1 balance the fixed point lies on the true profile and qw is
exact. Errors in the other balance terms (mut_1 prescription,
under-resolved dissipation in a coarse cell) shift T_1, and hence the
inferred qw, boundedly — the ordinary wall-function closure error present
in every code, including the paper's FD solver.

This picture also explains the failure signature of the pinned state: a
target-tracking scheme would fail gracefully (slightly wrong T_1);
self-consistency fails *structurally*. When no valid fixed point exists —
required flux and secant sign-incompatible while alphat >= 0 — the system
does not settle near a slightly-wrong state but slides to the constraint
boundary T_1 = T_w, where the closure inputs (h_1 = h_w) become garbage.

### Attribution

| layer                          | verdict |
|--------------------------------|---------|
| Li et al. laws (Eqs. 12-42)    | Sound here. Eq. (40) with h_rg > h_w describes the sublayer T-peak; correct inputs give sensible qw. |
| Paper's coupling (Fig. 16)     | Sound. Direct (Neumann-style) flux replacement in the discretized energy equation; no positivity/sign constraint. |
| FVM in general                 | Not the issue. A boundary flux is the natural FV boundary condition and can be imposed directly. |
| **OpenFOAM alphat inversion**  | **The failure point.** Any code (FD/FV/FE) delivering the flux as (positive diffusivity) x (resolved gradient) fails identically once y+_1 exceeds the T-peak location. |

### Regime boundary

The general envelope condition is **T(y+_1) > T_w** (secant sign compatible
with a cooled-wall flux) — *not* "first cell below the T-peak" per se.
Every cooled wall has a temperature peak; what matters is where the profile
falls back below T_w again:

- For a strongly cooled wall the profile stays above T_w far beyond the
  peak (generic example: Tw = 101 K, peak 148 K at z+ ~ 20, recrosses T_w
  only at z+ ~ 250 -> alphat coupling well-posed for y+_1 up to ~200).
- For near-adiabatic M6Tw076 (Theta = 0.76) the profile recrosses T_w
  already at y+ ~ 4-5 (peak 300.4 K at y+ = 2.6), so only y+_1 <~ 3 is
  inside the envelope.

Pre-run check: from the T-u relation the peak sits at
s* = (1 + A - Tw/Te)/(2A) in u/u_e; more directly, find where the relation
returns T = T_w and keep u+_1 (hence y+_1) left of it.

Observed in the y+ sweep for M6Tw076 (all vs. wall-resolved SA):

- y+_1 ~ 0.5, 0.9, 1.8: Li WF collapses onto resolved SA. Works.
- y+_1 ~ 6: degraded (T_1 near/above the T_w-recrossing).
- y+_1 ~ 9.5: broken (pinned state, both flux targets).

For strongly cooled walls (Theta <~ 0.3, the paper's entire test envelope,
Theta <= 0.2) T_1 > T_w holds at any sane y+_1 and the alphat inversion is
well-posed. The pathology is specific to walls near recovery temperature at
coarse y+_1. (Walls *above* recovery temperature would flip both signs
consistently — flux out of the wall, T_1 < T_w — and be well-posed again,
but are currently blocked by the separate rg <= 0 guard.)

## 3. Secondary finding: qHalf vs qw flux target

- The paper's q_1/2 is the correct total heat flux **at y_1/2**, the FD
  half-node face. An OpenFOAM boundary face sits at y = 0, where Eq. (34)
  gives exactly q_w. So `qw` is both the physical *and* the discretely
  consistent target in a cell-centred FV code; `qHalf` applies the right
  number at the wrong location.
- For near-adiabatic cases the difference is O(1): here
  (qw - qHalf)/qw = (r/Pr_t)*u+_1/2*u_tau^2/(B_q*cp*Tw) ~ 0.34*u+_1/2,
  i.e. qHalf is sign-reversed already at y+_1 ~ 8. For the paper's
  Theta ~ 0.1 cases the offset is a small fraction of qw, which is why the
  paper never noticed.
- Everything that reads the boundary flux (`wallHeatFlux` FO, Ch, CHT)
  measures the *applied target*; with `qHalf` the reported quantity is not
  the wall heat flux.
- **But**: selecting `qw` does not rescue the near-adiabatic case, because
  the delivery mechanism (Sec. 2) is the binding constraint. `qw` should be
  the default once the coupling is fixed.

## 4. Workaround validated for M6Tw076: adiabatic WF mode with isothermal T-BC

Because -B_q = 0.02 (wall within ~4 K of its local recovery state), run
`thermalMode adiabatic;` on the Li BCs while keeping `T` as
`fixedValue 300` (the WF needs the true Tw for rho_w, mu_w, h_w):

- Skips the Eq. (40) inversion -> the pinned-h1 attractor cannot form.
- alphat = 0 -> T_1 is finally free to float below T_w to its physical
  value (~266 K at y+ ~ 10).
- Residual errors: neglected wall heat extraction and a wrong-signed
  laminar leak (~1.7 kW/m^2), each ~2 % of the BL enthalpy flux; momentum
  side loses the qw share of the recovery factor (Gamma* low ~25-30 %) ->
  expected u_tau error ~1-2 % (vs. 15-20 % broken isothermal coupling).
- **Do not read `wallHeatFlux`** in this mode (it shows the meaningless
  laminar leak). Obtain qw by post-processing converged cell-1/cell-2
  states through Eq. (40)->(41)->(33). Intrinsic accuracy floor ~20 %
  at Theta = 0.76 (qw is a small residual of large terms; relation errors
  amplify by ~ rg/rg_heatflux ~ 4x).
- Case-specific trick, justified by B_q ~ 0. Not a general recipe.

## 5. Proper fixes (library work items)

1. **Rework the energy coupling to impose flux directly** (restores the
   paper's coupling semantics). Four delivery routes — Neumann BC,
   deferred-correction source, solver-level substitution, signed alphat —
   are compared in "Option survey" below; recommended: Neumann BC
   (option 1), fallback: deferred-correction source (option 2).
2. **Change the `auto` flux-target default to `qw`** (keep rejecting `qHalf`
   for adiabatic walls; document `qHalf` as legacy/diagnostic).
3. **Add a degeneracy guard**: currently the pinned state reports status 0
   everywhere. Flag faces where |h1 - hw| is tiny while the flux target is
   significant, or where target flux and snGrad(he) have incompatible signs
   ("outside coupling envelope"), and warn loudly.
4. `thermalMode auto` misclassifies fixed-flux T BCs (e.g.
   externalWallHeatFluxTemperature) as adiabatic (`!fixesValue()` test);
   restrict or warn.
5. Documented earlier (DESIGN_AND_DEVIATIONS review): omega* = 0.72 is
   coupled to the paper's enthalpy datum; rg <= 0 (heated walls) silently
   falls back; recovery factor r evaluated with wall Pr.

### Option survey: imposing the wall heat flux without a diffusivity

The paper itself never mentions a diffusivity: Fig. 16 overwrites tau_1/2
and q_1/2 *directly in the discretized viscous-flux routine* of their FD
solver. The alphat route is purely this port's packaging choice ("wall
function = boundary field on nut/alphat" is the OpenFOAM idiom). The
delivery options:

| # | option | mechanism | flux exact? | free of the sign constraint? | `wallHeatFlux` reports qw? | matrix risk | effort / portability |
|---|--------|-----------|-------------|------------------------------|----------------------------|-------------|----------------------|
| 1 | Neumann (prescribed-gradient) BC on T/he | set snGrad = qTarget/alphaEff each iteration (`externalWallHeatFluxTemperature` pattern); `fvc::laplacian` and `wallHeatFlux` read the *prescribed* gradient, not the secant | yes | yes | yes | none (standard machinery) | new T-BC class + plumbing: WF must take Tw (and wall transport props) from its dict because the patch T *value* drifts (see guide) |
| 2 | Deferred-correction source (fvOptions) | benign wall BC + volumetric source of the flux deficit in cell 1; variant: zero out alphaEff at the face and source the *full* qw | yes at convergence | yes | **no** (part of the heat enters as a source; report qw via WF diagnostics) | none | low; no new BC classes; source lags one iteration |
| 3 | Solver-level flux substitution | overwrite the wall-face viscous energy flux in HiSA `residualsUpdate.H` — the literal Fig. 16 analog | yes | yes | yes | none | HiSA fork; dies on solver updates; least portable |
| 4 | Signed alphat | drop the >= 0 clamp in the existing inversion | yes while snGrad != 0 | yes, but singular at snGrad = 0 (needs cap + relaxation) | yes | anti-diffusive boundary row (see below) | smallest code diff |

Recommendation: **option 1 long-term** (exact, sign-free, correct
reporting, standard fvPatchField machinery; cost is plumbing, not
numerics — see `IMPLEMENTATION_GUIDE_OPTION1_neumann_energy_BC.md`),
option 2 as the low-effort fallback, option 4 as the quick experiment,
option 3 only if a HiSA fork is maintained anyway.

Why momentum never had this problem: tau_w and the wall-tangential
velocity gradient always point the same way (the wall always retards the
flow), so the nut inversion never meets a sign conflict. It is uniquely
the energy equation — with its interior heat-generation peak — where the
diffusivity idiom can become unrepresentable.

### The signed-alphat fix in detail (option 4 above)

**The idea.** The structural failure is purely the alphat >= 0 constraint:
with T_1 < T_w the secant points out of the wall while the flux must point
into it. Allowing *signed* alphat makes the representable set "any flux for
any nonzero secant" — a heat flux into the wall coexists with a
positive-outward gradient. The physical fixed point then exists, T_1 is
free to settle below T_w, and the Eq. (40) inputs stay healthy.

**How mild it actually is.** At the healthy M6Tw076 fixed point
(T_1 ~ 266 K, T_w = 300 K, y_1 = 0.5 mm): snGrad(he) ~ c_v * 34 K / 5e-4 m
~ 4.9e7 J/kg/m, so delivering qw = 2142 W/m^2 needs
alphaEff ~ -4.4e-5 kg/(m s) — about **1.7x laminar alpha in magnitude,
just negative**. The bookkeeping coefficient is tame at convergence; it is
not a monstrous anti-diffusivity.

**Repercussion 1 — the singularity moves, it does not vanish.** No finite
alphat delivers a nonzero flux through snGrad = 0. Previously T_1 = T_w
was an unreachable-from-below attractor; with signed alphat it becomes a
point the transient must *pass through* (initial fields usually have
T_1 > T_w). Requires a |alphat| magnitude cap; near the crossing the
delivered flux transiently undershoots (harmless here: the top-face
exchange with the colder cell 2 is what drags T_1 down; extraction resumes
below T_w). Cases sitting *at* the envelope edge (true T_1 ~ T_w, y+_1 at
the T_w-recrossing) will hover near the singularity and may flicker
face-to-face — cap plus under-relaxation of the alphat update needed.

**Repercussion 2 — anti-diffusion in the energy matrix.** A negative
boundary conductance contributes with the wrong sign to the cell-1
diagonal/Jacobian (HiSA's viscous Jacobian uses boundary alphaEff
directly). With alphaEff frozen within an outer iteration, the linearized
feedback inverts: "cell warms -> wall extracts less" — locally
destabilizing. Mitigations: the magnitude is only O(laminar) at the fixed
point, and pseudo-transient continuation adds diagonal damping. Expected
outcome: converges, possibly with slower energy residuals / transient wall
wiggle. If the implicit solve objects, fall back to the
deferred-correction source (benign alphat + explicit flux-deficit source
in cell 1): identical converged cell balance, no matrix pathology, but the
`wallHeatFlux` FO then under-reports (part of the heat enters as a
source) and qw must be reported from the WF diagnostics instead.

**Repercussion 3 — ecosystem assumptions.** alphat >= 0 is an unwritten
invariant in OpenFOAM: CHT couplings, Nusselt post-processing, or other
BCs doing their own alphaEff arithmetic may silently misbehave on a
negative patch value. Non-issue in a standalone HiSA aero case, but the
feature must be gated behind an explicit flag (e.g. `allowSignedAlphat`),
never default.

**Recipe.** Flag-gated signed alphat + |alphat| cap + relaxed updates;
validate on M6Tw076 at y+_1 ~ 10 against the wall-resolved SA reference.
Success signatures in the diagnostics: alphat smoothly negative and
O(1e-5) along the plate; h_1 - h_w ~ -34 kJ/kg; model qw dropping from
~8000 to the ~2000 W/m^2 range. All three together confirm the structural
fix behaves as the theory says.

## 6. Method-level (not implementation) caveats that remain

- Near-adiabatic walls are outside the paper's validated envelope
  (Theta <= 0.2). Nothing in the derivation breaks, but accuracy is
  unproven there.
- At Theta -> 1, qw is a small residual: relative qw errors are amplified
  ~4x (for M6Tw076) from enthalpy-relation errors regardless of how the
  flux is imposed. A ~2 % T-u relation error -> ~10-30 % qw error. This is
  an intrinsic accuracy floor, not a well-posedness problem.
- Expectation confirmed by the y+ sweep: with y+_1 <= ~2-3 (below the
  T-peak) the WF run reproduces wall-resolved SA; "as good as the
  underlying turbulence model" holds for Cf, and for qw within the
  amplification floor above.
