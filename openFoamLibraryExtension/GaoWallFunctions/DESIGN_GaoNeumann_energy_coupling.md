# GaoNeumann Wall Functions - Neumann Energy Coupling Design Notes

This documents the `GaoNeumann*` boundary-condition family: a copy of the
Gao et al. wall function (see `DESIGN_AND_DEVIATIONS.md`, which remains
authoritative for the closure itself) whose energy coupling imposes the
wall heat flux directly through a prescribed-gradient (Neumann)
temperature boundary condition instead of packaging it into a
non-negative `alphat`. It is the Gao analogue of the `LiNeumann*` family
in `../LiWallFunctions` (see `DESIGN_LiNeumann_energy_coupling.md` there
for the full delivery-mechanism analysis, verified against the v2412 and
HiSA sources); the coupling mechanics are identical and are only
summarised here. The original `*GaoWallFunction` conditions are unchanged
and remain the right choice for strongly cooled walls.

## Why

`FINDINGS_M6Tw076_energy_coupling.md` (in `../LiWallFunctions`): the
OpenFOAM idiom `flux = alphaEff * snGrad(he)` with `alphat >= 0` cannot
represent a heat flux into the wall when the first cell is colder than
the wall (near-adiabatic walls at coarse y+, where the temperature peak
sits inside the sublayer). The coupled system then pins `T_1 = T_w` and
corrupts the closure inputs. The Gao closure consumes exactly the same
kind of inputs — `beta` is inverted from the Crocco-Busemann temperature
law (Gao Eq. 6) using `T_cell` and `T_w`, and
`q_w = beta*rho_w*T_w*k_w*u_tau/(1.2*mu_w)` (Gao Eq. 16) — so the alphat
version fails in the same regime the same way.

## Classes

| field    | boundary condition | notes |
|----------|--------------------|-------|
| `T`      | `GaoNeumannWallHeatFluxTemperature` | derives from `fixedGradient`; computes q_w each outer iteration and prescribes the wall-normal T gradient |
| `alphat` | plain `fixedValue 0`                | no Gao alphat class in this variant |
| `nut`    | `nutGaoNeumannWallFunction`         | identical momentum logic to `nutGaoWallFunction` |
| `omega`  | `omegaGaoNeumannWallFunction`       | identical Nichols Eqs. 19-21 logic |
| `k`      | `kGaoNeumannWallFunction`           | identical Nichols Eq. 22 logic |

Shared evaluator: namespace `GaoNeumannWallFunction`
(`GaoNeumannWallFunction/`), a copy of `GaoWallFunction` extended with
the wall-state plumbing below. All conditions on a patch must carry
identical coefficients, including `Tw` and `wallTransport`.

## The coupling (same as LiNeumann)

Each outer iteration the T condition places the face enthalpy so that the
discrete wall diffusion flux equals the target exactly:

```text
qTargetOF   = -q_w                       // paper sign -> outward patch normal
he_f,needed = he_1 + qTargetOF/(alphaEff*deltaCoeffs)
T_f,needed  = THE(he_f,needed, p_f)      // thermo's own inverse
gradient()  = (T_f,needed - T_1)*deltaCoeffs
```

`alphaEff` (the assembly's `CpByCpv*(alpha + alphat)`, with the wall
alphat fixed at 0) appears once in the inversion and once in the flux, so
its evaluation at the drifted patch state cancels; there is no sign
constraint, and `T_1 < T_w` with heat flowing into the wall is
representable. The delivered energy gradient is mirrored onto the paired
`gradientEnergy` he patch so the `wallHeatFlux` function object reports
the imposed q_w. `thermalFluxTarget` defaults to `qw`; `qHalf` remains a
diagnostic option (rejected for adiabatic walls). Adiabatic mode
prescribes a zero gradient — identical to `zeroGradient` — and keeps the
Gao Eq. 7 recovery-temperature evaluation of the closure.

## Wall state: Tw and wallTransport

Under this coupling the patch T value is the tangent-extrapolated face
temperature, NOT the physical wall temperature, and everything evaluated
at patch T (patch mu, kappa, alpha, rho, yPlus-style post-processing)
sees the drifted value. Therefore:

- `Tw` (dictionary, required for isothermal walls) supplies the physical
  wall temperature used throughout the closure (`beta`, `Gamma`, q_w,
  `rho_w` via `rhoEoS`, `cp_w`).
- `wallTransport` supplies `mu_w`/`k_w`: `sutherland` (`As`, `Ts`;
  `kappa` from `Pr` if given, else modified Eucken
  `mu*Cv*(1.32 + 1.77*R/Cv)`) or `constant` (`mu`, `kappa`). The default
  `fromPatch` is refused for isothermal walls under this coupling and
  remains valid for adiabatic mode, where the floating patch T is
  physical. The Nichols `omega_i` and the Eq. 8 `mu_t` also use the
  physical `mu_w`.
- The momentum inversion to `nut` deliberately keeps the *patch-field*
  `mu` and `rho` (the values OpenFOAM's stress operator uses), while
  `u_tau`, `yPlus`, `Pr_w`, and q_w use the physical wall state. In the
  legacy conditions patch mu equals `mu_w` because the patch T is pinned
  at `T_w`; here they differ.
- `thermalMode auto` recognises the new T condition (its
  `adiabaticIntent()`), so nut/omega/k do not silently misclassify the
  wall as adiabatic.

## Startup guards

Same knobs as LiNeumann, same reasoning: `gradient` starts at zero when
absent (adiabatic start), `TFaceMin`/`TFaceMax` (defaults 1/GREAT) clamp
the extrapolated face temperature (a violent transient can otherwise
demand a negative face temperature and crash the thermo), and
`gradientRelaxation` (default 1 = off) under-relaxes the gradient update.
While a clamp binds, the delivered flux undershoots the target — visible
as `qError` in the diagnostics. A prescribed gradient enters the energy
matrix as an explicit flux with no diagonal feedback, so crude startup
fields need the relaxation/clamp aids (or the solver's pseudo-transient
damping); the alphat coupling is self-limiting in this respect and does
not.

## Diagnostics

Placed on the T condition (same switches as the legacy alphat condition),
written to `postProcessing/GaoNeumannWallFunctionDiagnostics/<patch>/
<time>.dat` with columns:

```text
face x y z converged iterations residual uTau tauw yPlus uPlus Gamma beta
Tw cp rhoW qw qHalf uHalf nutShear mutEq8 omega k gradT TFace alphaEff
heGradient qTargetOF qOpenFOAM qError
```

`qOpenFOAM = alphaEff*heGradient` is the flux the assembly will deliver;
`qError` is nonzero only where the clamp or relaxation intervened.

## Example boundary conditions

```foam
// system/controlDict
libs ("libGaoWallFunctions.so");

// 0/T
wall
{
    type        GaoNeumannWallHeatFluxTemperature;
    thermalMode isothermal;         // or adiabatic (== zeroGradient)
    Tw          300;                // physical wall temperature [K]
    wallTransport
    {
        mode    sutherland;         // or constant (mu/kappa entries)
        As      1.458e-06;
        Ts      110.4;
        // Pr   0.72;               // optional; omitted -> modified Eucken
    }
    thermalFluxTarget qw;           // qw (default) | qHalf (diagnostic)
    // Startup/stability aids (optional):
    // gradientRelaxation 0.3;
    // TFaceMin           100;
    // TFaceMax           5000;
    diagnostics false;
    value       uniform 300;        // initial face value only
}

// 0/alphat
wall { type fixedValue; value uniform 0; }

// 0/nut
wall
{
    type        nutGaoNeumannWallFunction;
    thermalMode isothermal;
    Tw          300;
    wallTransport {mode sutherland; As 1.458e-06; Ts 110.4;}
    value       uniform 0;
}

// 0/omega
wall
{
    type        omegaGaoNeumannWallFunction;
    thermalMode isothermal;
    Tw          300;
    wallTransport {mode sutherland; As 1.458e-06; Ts 110.4;}
    value       uniform 1;
}

// 0/k
wall
{
    type        kGaoNeumannWallFunction;
    thermalMode isothermal;
    Tw          300;
    wallTransport {mode sutherland; As 1.458e-06; Ts 110.4;}
    value       uniform 1e-12;
}
```

## Status

- `wmake libso` clean; the original `*GaoWallFunction` classes are
  byte-identical (only `Make/files`/`Make/options` gained entries).
- Smoke tests (rhoSimpleFoam, kOmegaSST, coarse flat-channel toy):
  - isothermal GaoNeumann: runs; on faces where the `TFaceMin` clamp is
    inactive the delivered flux matches the target to round-off;
  - adiabatic GaoNeumann: `gradT` identically zero (== zeroGradient);
  - legacy Gao alphat case on the same setup: unchanged behaviour,
    `qError = 0` at its (well-posed, strongly cooled) state.
- No validation against the Gao flat-plate cases (`WF_yPlus_47_91`
  family) was run here; healthy-regime regression against the legacy
  alphat coupling and a near-adiabatic coarse-y+ case are the intended
  acceptance tests, mirroring the LiNeumann protocol.
