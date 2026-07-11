# WangWallFunctions — faithful OpenFOAM v2412 library for the Wang et al. AWF^Cmp with Neumann heat-flux coupling

## Context

Implement, as faithfully as possible, the compressible Analytical Wall Function of
**Wang, Chen, Craft & Iacovides (2023), "Extension of Analytical Wall Functions to
Supersonic and Hypersonic Flows", Flow Turb. Combust. 111:1-34** (the "paper",
`Wang et al.pdf`), which extends Craft et al. (2002) (`Craft et al.pdf`) and is fully
derived in Wang's PhD thesis (`Wang Thesis.pdf`). The final model is the paper's
**AWF^Cmp** = thesis **hyper-CMAWF**: MAWF (convection scaled y*²) + total-energy
dissipation term + hyperbolic molecular-viscosity variation in the sublayer.

Heat flux is imposed through a **Neumann (fixedGradient) T boundary condition**,
cloning the delivery mechanism of `../GaoWallFunctions` (`GaoNeumannWallHeatFluxTemperature`),
because the alphat-inversion coupling cannot represent flux into the wall when T_1 < T_w.

User decisions (already made):
- **Turbulence models: both** — faithful standard k-ε variant AND a k-ω SST translation (flagged deviation).
- **Solver: HiSA** (user's `hysa` fork at `~/OpenFOAM/mgnadlinger-v2412/modules/hysa`). Verified:
  `src/modules/hisa/residualsUpdate.H:35-37` adds `fvOptions(rho, thermo.he())` into the ρE residual → fvOption energy sources work.
- **CMAWF-D near-wall-cell dissipation replacement: fvOption** correction source.
- **Viscosity variation: hyperbolic** (paper form); `bMu = 0` limit recovers constant-μ CMAWF (runtime switch).

## Source-equation map (for implementation-time re-reading)

| Content | Where |
|---|---|
| Model overview, Eqs. 1–29 (SWF, AWF, MAWF, CMAWF, μ(T)) | paper PDF pages 6–13 |
| **Appendix a/b: full AWF^Cmp formulas, Eqs. 30–45** | paper PDF pages 29–32 |
| Simple AWF + P̄k/ε̄ (Eqs. 6.10–6.52) | thesis PDF pages 79–88 |
| MAWF (Eqs. 8.1–8.39) | thesis PDF pages 119–127 |
| CMAWF: dissipation in energy eq, N1..N7, Eq. 8.64 | thesis PDF pages 128–136 |
| hyper-CMAWF (Eqs. 8.104–8.130) — matches paper appendix | thesis PDF pages 147–152 |
| Thin-cell fallback y*n < y*v: Appendix A (AWF), C (CMAWF), E (hyper) | thesis PDF pages 250–258 |
| Craft k/ε near-wall treatment (Eqs. 4.6–4.9, App. A/B) | Craft PDF pages 7–12 |

## Constants (dictionary-adjustable, defaults faithful to the paper)

`yStarV = 10.8` (paper/Craft; thesis 6.2.3 prints 10.1 — treated as typo),
`alpha = cmu*cl = 0.2295`, `cl = 2.55`, `yStarD = 2*cl = 5.1`, `Cmu = 0.09`,
`Prt = 0.9`, molecular `Pr` from wallTransport (paper used 0.7 const — allow override),
`kappa=0.4178`, `E=9.793` (only for optional SWF-style diagnostics; not used by the model).

## Library layout (`src/WangWallFunction`, lib `libWangWallFunctions.so`)

```
Make/files, Make/options            (clone Gao Make/options; add -I to hisa lnInclude NOT needed — BC-level only)
WangNeumannWallFunction/            shared evaluator (namespace, Coeffs + PatchData + calcPatch)
WangNeumannWallHeatFluxTemperature/ T BC (fixedGradient)  — the Neumann coupling
nutWangNeumannWallFunction/         nut BC — imposes tau_w through OF's stress path
epsilonWangNeumannWallFunction/     k-ε variant: fixes ε_P, injects P̄k, ε̄ via G
omegaWangNeumannWallFunction/       k-ω SST translation (flagged deviation)
WangNeumannDissipationSource/       fv::option — CMAWF-D cell-averaged dissipation correction
DESIGN_AND_DEVIATIONS.md            document every mapping choice (Gao-style)
```
alphat at the wall: plain `fixedValue 0` (as in GaoNeumann design). k at the wall: standard
`kqRWallFunction` (zero-gradient; Craft: molecular diffusion of k to wall neglected, k~y²).

## Shared evaluator `WangNeumannWallFunction::calcPatch`

### Inputs / wall state (Gao-Neumann plumbing, reuse pattern from `GaoNeumannWallFunction.C`)
- `Coeffs` from dict: model constants above + `thermalMode auto|isothermal|adiabatic`,
  `Tw` (physical wall T, required isothermal), `wallTransport {mode sutherland|constant; As; Ts; [Pr]}`
  (patch T drifts under the Neumann coupling → patch props unusable, same errors/refusals as Gao),
  `nIntegrationPoints` (default 40), `TvRelax` (default 0.5), `bMuZero` switch (constant-μ CMAWF).
- Per patch: kP = k[celli] (turbModel field via `turbModel.k()`), y_P = `turbModel.y()`,
  ρ from `turbModel.rho()`, p_wall from patch p, μw = wallTransport(Tw), cp = thermo.Cp(p_w,Tw),
  ρ_w-style EoS calls via `thermo.rhoEoS` (as Gao).

### North-face state (Un, Tn, y_n) — Craft Fig. 4
For each wall face: find the cell face of `celli` most anti-parallel to the wall face normal
(the "north" face); `Un` = wall-parallel velocity magnitude and `Tn` = T interpolated to that
face (linear owner/neighbour weights); `y_n` = wall-normal distance of that face centre.
Fallback for degenerate cells: `y_n = 2*y_P`, `Un = |U_P|·(y_n/y_P interpolation) → just U_P`, `Tn = T_P`.
**Mapping choice to document** (paper: structured wall-aligned grid).

### Sublayer-edge state Tv, ρv, μv (paper Eqs. 25–28)
Fixed-point loop (3 iterations, under-relaxed by `TvRelax`, initialised from T_P; store previous
Tv per face in the T BC for warm starts): evaluate thermal profile → Tv = T1(y*v) → μv = wallTransport(Tv),
ρv = thermo.rhoEoS(p_w, Tv) → y* = y·√kP·ρv/μv (Eq. 26), μt = μv·α(y*−y*v) (Eq. 27),
b_μ = (μw − μv)/(μw·y*v) (Eq. 28; 0 if `bMuZero`).

### Nodal convection/pressure coefficients (wall-aligned frame; D1=D2=D, C1=C2=C, Dth1=Dth2=Dth)
t̂ = unit tangential projection of U_P onto the wall plane; n̂ = inward wall normal.
Using cached volFields (computed once per patch update): `gradU = fvc::grad(U)`, `gradp = fvc::grad(p)`,
`divPhiE = fvc::div(((rho*E + p)*U))` with `E = e + ½|U|²`, `e` from thermo
(HiSA: he is sensible internal energy; if enthalpy-based thermo, e = he − p/ρ).
- Momentum (paper Eq. 18, non-conservative nodal form of thesis 8.105):
  `conv = ρ_P·(U_t·(t̂·gradU·t̂) + U_n·(n̂·gradU·t̂))`,
  `D = μv²/(ρv²·kP·y*n²)·conv`, `C = μv²/(ρv²·kP)·(t̂·gradp)`.
- Energy (paper Eq. 23): `Dth = μv²/(Cp·ρv²·kP·y*n²)·divPhiE[celli]`.
Document: paper scales BOTH regions by 1/y*n²; thesis Eqs. 8.50/8.78/8.105 print 1/y*v² for the
region-1 coefficients — inconsistent with RHS continuity at y*v and with the peer-reviewed paper;
follow the paper (suspected thesis typo), note in DESIGN doc.

### Momentum solution (paper Eqs. 30–36; thesis 8.104–8.113)
With Y = 1+α(y*−y*v), Yn = Y(y*n):
- `A1 = (μv·Un − N) / (lnYn/α + y*v − 0.5·bμ·y*v²)`  (Eq. 32)
- N per Eq. 35 (transcribe exactly; with D1=D2, C1=C2 the (D1−D2), (C1−C2) terms vanish but keep them in code for clarity)
- `A2 = (D1−D2)y*v³/3 + (C1−C2)y*v + A1`, B2 per Eq. 34.
- **τ_w = −(ρv·√kP/μv)·A1** (thesis Eq. 8.12; μ∂U1/∂y*|₀ = A1).
- Profiles: U1(y*) Eq. 30, U2(y*) Eq. 31; first-integral identities used everywhere:
  `μ·U1' = D1y*³/3 + C1y* + A1` (sublayer), `μvY·U2' = D2y*³/3 + C2y* + A2` (outer).

### k-equation quantities
- `P̄k = (1/y_n)·(ρv√kP/μv)·∫_{y*v}^{y*n} μv·α(y*−y*v)·(∂U2/∂y*)² dy*` (paper Eq. 13/thesis 8.13),
  Simpson with `nIntegrationPoints`; zero if y*n ≤ y*v.
- `ε̄ = (1/y_n)·[2·kP^{3/2}/y*d + (kP^{3/2}/cl)·ln(y*n/y*d)]` (paper Eq. 16); if y*n ≤ y*d:
  `ε̄ = 2·νv·kP/y_d²` with `y_d = y*d·μv/(ρv√kP)` (Craft Eq. 4.8). No Craft Fε factor (Wang omits it).
- `ε_P = kP^{3/2}/(cl·y_P)` (Craft Eq. 4.9, Fε=1).

### Thermal solution (paper Eqs. 37–45; thesis 8.114–8.130)
αt = α·Pr/Prt, YT = 1+αt(y*−y*v), YTn = YT(y*n).
Dissipation integrals (Simpson, exploiting first-integral identities; per-face 1-D quadrature):
- `I1(b) = (1/Cp)·∫₀^b (1+bμ(s−y*v))·U1(s)·(D1s³/3 + C1s + A1) ds`
- `I2(b) = (1/Cp)·∫_{y*v}^b (1/YT(s))·U2(s)·(D2s³/3 + C2s + A2) ds`
- `BT = Dth1·y*v⁴/12 − (Dth1/60)·bμ·y*v⁵ − I1(y*v) − (Dth2/3)(1/αt³)·[αt²/3·y*v³ − αt(1−αt·y*v)/2·y*v² + (1−αt·y*v)²·y*v]` (paper Eq. 41)
- `Ath2 = (Dth1/3 − Dth2/3)y*v³ + Ath1` (Eq. 40), `Bth2 = Ath1·y*v − ½bμ·Ath1·y*v² + BT + (μv/Pr)·Tw` (Eq. 39).
- **Isothermal** (Tw given): Ath1 from paper Eq. 44:
  `Ath1 = { (μv/Pr)(Tn−Tw) − (Dth2/3)(1/αt³)[αt²/3·y*n³ − αt(1−αt y*v)/2·y*n² + (1−αt y*v)²y*n − (1−αt y*v)³(1/αt)lnYTn] − (lnYTn/αt)(Dth1/3 − Dth2/3)y*v³ − BT + I2(y*n) } / (lnYTn/αt + y*v − ½bμ·y*v²)`
  then **q_w = −(ρv·cp·√kP/μv)·Ath1** (Eq. 43). Paper sign: q_w positive from wall into fluid;
  OpenFOAM target flux is `−q_w` (Gao convention, keep identical).
- **Adiabatic / prescribed q_w**: `Ath1 = −q_w·μv/(cp·ρv·√kP)` (=0 adiabatic); wall temperature
  from paper Eq. 42 (thesis 8.128) — needed for the wall-state properties loop.
- Tv for the property loop: evaluate T2(y*v) = Pr/μv·{Ath2·... } — cheaper: T1(y*v) via Eq. 37
  with I1(y*v) (equal by continuity).

### Thin-cell fallback y*n < y*v (thesis Appendix E)
- `bμ = (μw−μn)/(μw·y*n)` with μn at Tn (E.1); μt = 0 everywhere in the cell.
- `A1 = [μn·Un − (D1y*n⁴/12 + C1y*n²/2 − bμ·y*n³·(D1y*n²/60 + C1y*n/6))]/(y*n − ½bμ·y*n²)` (E.2 rearranged); τ_w as above; **P̄k = 0**.
- Thermal: E.3 (Tw for prescribed flux) / E.4 (Ath1 isothermal), with the I1-type integral over [0, y*n].
- ε̄ thin-cell branch as above.

### PatchData outputs
`tauw, qw, Tw, Tv, yStar, yStarN, A1, Ath1, Pk̄, ε̄, εP, dissAvg, muv, rhov, muw, cp, converged/iterations` (+ diagnostics fields).

## Boundary conditions

### `WangNeumannWallHeatFluxTemperature` (T, derives fixedGradient)
Direct clone of `GaoNeumannWallHeatFluxTemperatureFvPatchScalarField` mechanics
(`GaoNeumannWallHeatFluxTemperature/...C:364-486`): compute q_w via calcPatch each outer
iteration; `heFaceNeeded = he_1 + (−q_w)/(alphaEff·deltaCoeffs)`; `TFace = thermo.THE(...)`
clamped `[TFaceMin, TFaceMax]`; `gradient() = (TFace − T_1)·deltaCoeffs` under `gradientRelaxation`;
mirror delivered he-gradient onto the paired `gradientEnergy` he patch; adiabatic → zero gradient;
same startup guards and face-wise diagnostics (dir `postProcessing/WangNeumannWallFunctionDiagnostics`).
Also stores the per-face `Tv` field (read/written for restart warm-start) shared by the property loop.
`adiabaticIntent()` hook so the other BCs resolve `thermalMode auto` (as Gao's `isAdiabatic`).

### `nutWangNeumannWallFunction` (nut)
calcPatch → τ_w; invert through OF's actual stress path (Gao pattern, `GaoNeumannWallFunction.C:712-726`):
`nut = max(0, (τ_w − μ_patch·|snGrad(U)|)/(ρ_patch·|snGrad(U)|))`, patch-field μ/ρ deliberately
(they drift; documented). `yPlus()` reporting via y*.

### `epsilonWangNeumannWallFunction` (ε; faithful variant)
Derive from `epsilonWallFunction` to reuse master-patch/corner machinery. In near-wall cells:
- fix `ε_P = kP^{3/2}/(cl·y_P)` (matrix constraint as the base class does),
- set `G[celli] = P̄k + ε_sink − ε̄` where `ε_sink = ε_P` — so the k-equation's net source is
  `P̄k − ε̄` while its implicit sink uses the fixed nodal ε_P. Exact at steady state; the implicit
  Sp linearisation difference is documented (deviation: OF's kEqn structure vs Wang's in-house code).

### `omegaWangNeumannWallFunction` (ω; **translation, flagged deviation** — not in the papers)
Derive from `omegaWallFunction`. `ω_P = ε_P/(Cμ·kP) = √kP/(Cμ·cl·y_P)`;
`G[celli] = P̄k + β*·ω_P·kP − ε̄` (β* = Cμ = 0.09, SST k-sink). No Nichols blending.

## `fv::WangNeumannDissipationSource` (fvOption, CMAWF-D term)

Paper Eq. 24 / thesis Eq. 8.64 — but note the exact simplification (document it):
`avg = (1/V_cell)·∫ ∂/∂y((μ+μt)U ∂U/∂y) dV = (1/y_n)·[(μ+μt)U ∂U/∂y]₀^{y_n} = Un·τ(y_n)/y_n`
since U(0)=0, with `τ(y_n) = (ρv√kP/μv)·(D2y*n³/3 + C2y*n + A2)` (or region-1 form if y*n<y*v).
The solver's own version of this term is embedded in HiSA's face-based `∇·(σ·u)` flux and cannot
be removed, so the option adds the **difference**:
`S[celli] = avg − localEstimate`, `localEstimate = n̂·∇(φ)` with volField
`φ = (μ+μt)·U_t·(n̂·∇U_t)` built from resolved fields (wall-normal direction from `wallDist`).
Registered on the energy field (`he`); HiSA picks it up via `fvOptions(rho, thermo.he())`
(`modules/hysa/src/modules/hisa/residualsUpdate.H:36`, `#else` branch active for v2412).
Own dict: `patches (...)` + same Coeffs entries (must match the BCs); `enabled` switch.
**This is the roughest mapping in the whole library — clearly marked in DESIGN doc; default off?
No: default ON (CMAWF-D is part of AWF^Cmp), but trivially disabled by omitting the fvOption.**

## Deviations / mapping choices to flag in DESIGN_AND_DEVIATIONS.md

1. Unstructured-FV frame: t̂/n̂ per face, nodal gradients via `fvc::grad`/`fvc::div`; energy convection
   uses full divergence of (ρE+p)u (equals the paper's 2-D x+y form for 2-D flows).
2. D1=D2, C1=C2, Dth1=Dth2 (single nodal values; paper Eq. 18/20/23). Thesis y*v² scaling of
   region-1 coefficients treated as typo (breaks RHS continuity at y*v; paper uses y*n²).
3. North-face (Un, Tn, y_n) reconstruction on general meshes + fallback.
4. Tv fixed-point (3 iters, relaxed, warm-started) — papers don't specify the update scheme.
5. ε̄ injected through G-trick (`G = P̄k + ε_P − ε̄`); implicit-sink linearisation differs from
   Wang's in-house discretisation. No Craft Fε laminarisation factor (Wang omits it).
6. ω variant entirely a translation (ω = ε/(Cμk)); untested by Wang.
7. fvOption dissipation source adds a *correction difference*, not a literal replacement, and its
   local-estimate is a cell-centred reconstruction of the wall-normal dissipation only.
8. y*v = 10.8 (paper) vs 10.1 (thesis 6.2.3) — dictionary-adjustable, default 10.8.
9. Neumann coupling inherits all Gao-Neumann caveats: drifted patch T, dictionary Tw + wallTransport
   required for isothermal walls, alphat pinned to 0, startup clamps/relaxation, sign convention −q_w.
10. Numerical guards (not physics): Y, YT, denominators, ln arguments bounded by SMALL; y* floors;
    exp/ln clamps; Simpson point count fixed.

## Build

`Make/files` (all six units) → `LIB = $(FOAM_USER_LIBBIN)/libWangWallFunctions`;
`Make/options` = Gao's includes/libs + `-I$(LIB_SRC)/finiteVolume/lnInclude` etc. (same list;
fvOptions needs `-lfvOptions`? — in v2412 `fv::option` lives in finiteVolume; verify at build).
Build: `source /opt/of-v2412/etc/bashrc && wmake libso`.

## Verification

1. `wmake libso` clean build.
2. **Algebra cross-check**: small `TestWangWallFunction` app (pattern: `LiWallFunctions/TestLiWallFunction`)
   + independent Python/numpy implementation of paper Eqs. 30–45 in the scratchpad; feed identical
   synthetic inputs (several y*n above/below y*v, heating & cooling, adverse/favourable dP/dx);
   assert machine-precision agreement of A1, N, τw, BT, Ath1, q_w, Tw, P̄k, ε̄ and Simpson convergence.
3. **Analytic limit checks** in the test app:
   - continuity: U1(y*v)=U2(y*v), μU1'(y*v)=μvY U2'(y*v), T1(y*v)=T2(y*v), U2(y*n)=Un, T2(y*n)=Tn;
   - bμ=0, D=C=Dth=0 reduces to Craft (2002) forced-convection formulas (App. A.5/A.11/A.12);
   - q_w from Eq. 43 with Ath1 from Eq. 44 reproduces Tn−Tw consistency (invert round-trip).
4. **HiSA smoke test**: coarse flat-plate case (adapt user's existing Gao/Li validation setup, e.g.
   `Gao_TurbZPGFP`), kOmegaSST + ω variant + T/nut/alphat BCs + fvOption: runs 100 iters, no NaNs,
   diagnostics show `qError ≈ 0` where clamps inactive, τ_w/q_w magnitudes plausible vs the Gao WF run.
5. k-ε variant smoke test with `kEpsilon` on the same mesh.
6. Restart/write-read round trip of the T BC (gradient + Tv fields).
