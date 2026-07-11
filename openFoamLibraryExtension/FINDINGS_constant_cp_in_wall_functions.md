# Findings: the constant-cp assumption in the Crocco–Busemann-type wall functions

Date: 2026-07-09. Scope: all four wall-function families under `src/`
(Gao, Wang, Griffin, Li). Question investigated: are Gao, Wang, and
Griffin constant-cp (calorically perfect) models, where exactly does cp
enter, what does the choice of cp value do, and is a path-averaged cp a
sensible improvement?

## TL;DR

- **Gao, Wang, and Griffin are all constant-cp models — in the papers
  and in these implementations.** The assumption is structural (it is
  what makes their closed-form solutions exist), not a property-
  evaluation detail. Li is the only family formulated in enthalpy and
  therefore the only one valid for variable cp.
- The implementations sample the single cp at **different states**:
  Gao/Wang/Li at the wall state (p_face, T_w), Griffin at the
  first-cell state. Inside the calorically perfect envelope these
  coincide; outside it they diverge — but that divergence is a symptom,
  not the disease.
- The physical picture: in a cold-wall hypersonic layer the temperature
  peaks *inside* the layer; a physical T peak is also a physical cp
  peak (vibrational excitation), which absorbs dissipation and lowers
  the real peak. **No single cp — cp(T_w), cp(T_1), or a mean — sees
  that interior feedback.**
- Of the single-value picks, the **path-mean
  c̄p = (h(p_w,T_1) − h(p_w,T_w))/(T_1 − T_w)** is the best: exact
  integral mean (h is the antiderivative of cp, two `thermo.he()`
  calls, no quadrature), endpoint-exact by construction, reduces the
  first-order sampling error (~Δcp/cp) to a second-order interior-shape
  error. It preserves every closed-form expression in all three
  families.
- Current stance: Li is the principled variable-cp answer but is **not
  yet generally trusted here**, and the other families may still
  predict better through other aspects of their modelling. c̄p is the
  pragmatic hedge for hot cases while keeping Gao/Wang/Griffin in play.

## Where each implementation evaluates cp

All four call `basicThermo::Cp(p, T, patchi)` — the mixture cp at a
caller-supplied (p, T) — with wall-face pressure
`p.boundaryField()[patchi]` (thin-BL constant pressure).

| Family | Call site | T used for cp | Iterated? |
|---|---|---|---|
| Gao (Dirichlet) | `GaoWallFunctions/GaoWallFunction/GaoWallFunction.C:381,400` | wall T (patch value) | 3-pass fixed point |
| Gao (Neumann) | `GaoWallFunctions/GaoNeumannWallFunction/GaoNeumannWallFunction.C:562,589` | wall T (dict `Tw`) | 3-pass fixed point |
| Wang | `WangWallFunction/WangNeumannWallFunction/WangNeumannWallFunction.C:1095` | wall T (dict `Tw`) | inside Tv loop |
| Griffin | `GriffinWallFunctions/GriffinWallFunction/GriffinWallFunction.C:776` | **first-cell T** | once |
| Li | `LiWallFunctions/LiWallFunction/LiWallFunction.C:514`, `LiNeumannWallFunction/LiNeumannWallFunction.C:612` | wall T | once |

Notes:

- Gao adiabatic mode is self-consistent: cp(p_w, T_w) alternates with
  the recovery update `Tw = Tcell + r·u²/(2cp)` (r = Pr^{1/3} built
  from that cp) over 3 iterations. The Dirichlet variant mixes
  wall-state cp with patch-field conductivity in Pr; the Neumann
  variant rebuilds mu/kappa at T_w from the `wallTransport` law each
  pass (incl. `thermo.Cv` for the modified-Eucken kappa).
- Wang deliberately evaluates rho_v, mu_v at the sublayer-edge T_v but
  cp at T_w (`cpw // cp at (p_wall, Tw)`); the variable-property
  treatment covers mu and rho only.
- Griffin's code comments it explicitly: "Calorically-perfect-gas form
  of the model: Cp evaluated at the matching (first-cell) state."
  Deliberate, since the patch T drifts under the Neumann coupling.
- Li also calls `thermo.he(pWall, Tw, patchi)` — its model works in the
  thermo's real h(p,T), which is exactly the variable-cp-general part.

## How cp enters the derivation (Gao/White chain, maps to `GaoWallFunction.C`)

cp enters at two different levels: once as a **coefficient** in the
final formulas (harmless, tunable) and once as an **enabling
assumption** in the derivation (structural, not tunable).

### Step 1: energy balance in the wall layer (no cp assumption yet)

Constant-stress layer (Couette approximation), both equations
integrate once in y:

```
mu_eff du/dy = tau_w                                (constant shear)
lambda_eff dT/dy + mu_eff u du/dy = -q_w            (constant total-energy flux)
```

Heat conducted toward the wall plus kinetic energy transported by
shear work is the same at every height — exact within the thin-layer
approximation.

### Step 2: ODE in u — the load-bearing cp step

Divide energy by momentum to eliminate y, use lambda = cp·mu/Pr:

```
cp · dT/du = -Pr_eff · (q_w/tau_w) - Pr_eff · u
```

**If cp is constant**, the left side integrates to cp·(T − T_w):

```
T(u) = T_w - (Pr·q_w)/(cp·tau_w) · u - Pr·u²/(2·cp)
```

— the Crocco–Busemann relation: **T is a downward parabola in u**.
The u² term is viscous heating, the linear term is wall heat flux
tilting the parabola. In White/Gao nondimensional form:

```
T/T_w = 1 + beta·u+ - Gamma·(u+)²
```

with `Gamma = r*sqr(uTau)/(2*cp*Tw)` and beta built from q_w and
k_w (= cp·mu_w/Pr) — the symbols in the code. Physical meaning of cp
here: **how much temperature rise the gas suffers per unit of
dissipated kinetic energy.** Big cp → flat parabola; small cp → hot
layer.

### Step 3: temperature feeds the velocity law via density

Ideal gas at constant wall-layer pressure: rho(u) = rho_w·T_w/T(u).
The compressible law of the wall is a Van Driest-type transformation:

```
int sqrt(rho/rho_w) du  =  int du / sqrt(1 + beta·u+ - Gamma·(u+)²)
```

**The integral of 1/sqrt(quadratic) is an arcsine.** That is the only
reason Gao Eq. 3–4 (White's law) has the closed
`asin((2*Gamma*u+ - beta)/Q)` form — the very expression the
implementation clamps to [−1, 1] in its numerical guards. The printed
wall law exists because T(u) is quadratic, which exists because cp
came out of the integral in Step 2.

### What changing the *value* of cp does (coefficient level)

- `Gamma ∝ 1/cp`; recovery temperature `T_aw = T_cell + r·u²/(2cp)`
  drops as cp grows. At M6 flight (u ≈ 1800 m/s): u²/(2cp) ≈ 1600 K
  with cold cp = 1005 vs ≈ 1430 K with vibrationally excited
  cp ≈ 1130 — ~10% in the driving temperature.
- Cold-wall chain for too-small cp: parabola too steep → mid-layer T
  over-predicted → density under-, viscosity over-predicted → solved
  u_tau shifts and (T_aw − T_w) too large → **q_w over-predicted**,
  roughly by the cp error (10–15% class at M6–8 flight).
- Sampling cp at wall vs first-cell vs freestream only moves you
  within this band; it changes beta, Gamma, r numerically but nothing
  structural.

### Why no "better" cp value fixes it (structural level)

With cp(T), Step 2 still integrates — but the left side becomes
∫cp dT = h(T) − h(T_w). The parabola moves to **enthalpy**:

```
h(u) = h_w - Pr·(q_w/tau_w)·u - Pr·u²/2
```

Then: T(u) = h⁻¹(parabola) is no longer quadratic; rho(u) is no longer
1/(quadratic); the Van Driest integral no longer yields an arcsine —
no closed form at all. cp's constancy is the algebraic property
(∫cp dT = cp·ΔT) that lets T stand in for h, keeps the parabola in T,
keeps density rational in u, and makes the arcsine exist. Break it and
the printed formulas are not "less accurate" — they are no longer
solutions of the model equations.

**Li's move is exactly the fix**: keep the parabola in h (Zhang's
generalized recovery factor handles the cold wall), invert through the
thermo's real h(p,T) at the end. Wang's `(1/Cp)` pulled outside the
WF22 integral and Griffin's quadratic-in-U relation (Eq. 2.8, a
Crocco–Busemann → Walz → Zhang descendant) are the **same Step-2 move
in different costumes — all three inherit the same envelope for the
same reason.** Retrofittability differs: Griffin (numerical inner
solve) could swap in h(U) + h⁻¹ most easily; Gao would lose the
arcsine; Wang is unretrofittable in spirit (analytic integrability is
the AWF's defining feature).

## Paper positions (verified in the PDFs)

- **Wang et al. 2023**: explicit — "working medium is assumed as a
  calorically perfect gas … Cp = 1004.06 J/(kgK)"; variable-property
  sophistication spent on mu(T), rho(T) only.
- **Griffin et al. (WM 2022)**: Cp a plain constant in Eq. 2.3/2.8,
  "ideal gas equation of state closes the system"; s = 1.14 and
  Pr^{1/3} calibrated against calorically perfect DNS.
- **Gao et al.**: temperature law "deduced from the Crocco-Busemann
  energy integration equation" — constant cp by construction.
- **Li et al. 2024**: notes the classical T–u relation "is derived
  based on the calorically perfect-gas assumption" and generalizes to
  an enthalpy–velocity relation; one-temperature 5-species air,
  derived for a noncatalytic wall, thermal nonequilibrium neglected.

## Validity envelope vs flight condition (~220 K ambient)

- Layer peak T < ~600–800 K (M ≲ 4–5 flight, or cryogenic-freestream
  tunnel/DNS conditions incl. the Xue M6Tw025/M6Tw076 validation
  cases): cp constant to ~1% — **all four families in envelope**; the
  M6 DNS validations are legitimate despite "Mach 6".
- ~800–2000 K (M ~ 5–8 flight): vibrational excitation, cp(T) rises
  ~15% by 1500 K, no chemistry — thermally perfect; needs an
  enthalpy-based wall law (Li, non-reacting limit) on a janaf thermo.
- > ~2000–2500 K (M ≳ 8–10, reentry): dissociation; finite-rate
  chemistry, wall catalycity first-order for q_w — beyond stock
  OpenFOAM/HiSA; Li's full HRTBL form in a reacting solver.

## The path-mean cp option (c̄p)

Analogous to path-averaging mu(T), k(T) over [T_w, T_1] (as seen in
another code's `sampleWallTransportPath`; for a constant-flux sublayer
the T-path mean of k is the *exact* effective conductivity). For cp no
quadrature is needed — h is the exact antiderivative:

```
c̄p = (h(p_w, T_1) − h(p_w, T_w)) / (T_1 − T_w)     // two thermo.he() calls
```

- **Endpoint-exact by construction**: c̄p·(T_1 − T_w) = h_1 − h_w, so
  the T–u relation becomes exact at the wall and the matching cell;
  the residual is the interior shape, second-order in Δcp/cp
  (≈ a percent or two on T where the sampling ambiguity was 10–15%).
- Slots into each family without touching the closed forms: Gao —
  replace cp in Gamma/beta/recovery (adiabatic path [T_cell, T_aw] is
  implicit but the existing 3-pass fixed point covers it); Wang — a
  different constant outside the WF22 integral, algebra unchanged
  (path [T_w, T_n]); Griffin — c̄p over [T_w, T_m] in Eq. 2.8 +
  recovery (mild off-calibration of s, Pr^{1/3}).
- **What it cannot fix**: the interior T(u) is still a parabola in T;
  the density weighting ∫sqrt(rho/rho_w)du integrates through the
  interior; the in-layer temperature peak (which can exceed both T_w
  and T_1, and is where cp deviates most) is not bracketed by
  [T_w, T_1]. The physical feedback — T peak → cp peak → real peak
  lowered — is only captured by the enthalpy formulation.
- **Validation caveat**: inside the calorically perfect envelope
  (cryogenic DNS cases) all cp picks coincide — the retrofit cannot be
  validated against the current reference cases; a genuinely hot case
  is required.

## Decision / next step

Li is the principled answer for variable cp but is not yet generally
trusted here, and Gao/Wang/Griffin may retain advantages from other
modelling aspects. Before implementing anything: run one hot
flight-condition case three ways — cp(T_w), cp(T_1), c̄p — via a
dictionary-overridable cp and compare the q_w spread against the other
error sources (Pr_t, transition, turbulence model). If the spread is a
few percent, skip the retrofit and reach for Li only where the
envelope demands it.
