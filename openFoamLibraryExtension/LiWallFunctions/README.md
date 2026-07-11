# Li Frozen-Composition Wall Functions

This directory contains an OpenFOAM v2412 user library implementing the Li et
al. improved wall function for a single-temperature, frozen-composition ideal
gas. Constant-property and temperature-dependent thermodynamics are supported.
The reacting-species part of the paper is deliberately not implemented.

The implementation uses the reference-point form in Li et al. Eqs. (40)-(42).
For each wall face, the first off-wall cell is the boundary-face owner and the
second point is the adjacent cell with the largest displacement in the inward
wall-normal direction. Processor-neighbour cells are included.

## Build

```bash
source /opt/of-v2412/etc/bashrc
cd "$WM_PROJECT_USER_DIR/src/LiWallFunctions"
wmake libso
```

This writes:

```text
$FOAM_USER_LIBBIN/libLiWallFunctions.so
```

Load the library from `system/controlDict`:

```foam
libs ("libLiWallFunctions.so");
```

## Boundary conditions

The v1 implementation is Spalart-Allmaras-only:

- `nutLiWallFunction` on `0/nut` imposes the Li wall shear through the
  OpenFOAM viscous-stress operator.
- `compressible::alphatLiWallFunction` on `0/alphat` imposes the Li energy
  flux through the OpenFOAM energy-diffusion operator.
- `nuTildaLiWallFunction` on `0/nuTilda` keeps `nuTilda=0` at the wall and
  defines the Li Eq. (39) first-cell matrix constraint. Stock OpenFOAM v2412
  SA does not invoke patch `manipulateMatrix()`, so the constraint is inactive
  unless `SpalartAllmarasBase::correct()` is changed to call
  `boundaryManipulate()` before solving.
  The local WSL `/opt/of-v2412` build contains this change as of 2026-07-11;
  other installations, including cluster builds, must be patched separately.

Complete snippets are in [examples/boundaryConditions.foam](examples/boundaryConditions.foam).

Shared defaults are:

```foam
kappa          0.41;
B              6.0;
Prt            0.9;
omegaH         0.72;
cv1            7.1;
enthalpyOffset 0;
maxIter        80;
tolerance      1e-8;
thermalMode    auto;
strict         false;
```

All three fields must use the same shared coefficients.

`thermalMode auto` treats a fixed-value temperature patch as isothermal and a
non-fixed-value temperature patch as adiabatic. `isothermal` and `adiabatic`
can be selected explicitly.

For `alphat`, `thermalFluxTarget auto` selects Li's numerical half-cell flux
`qHalf` for isothermal walls and the physical flux `qw=0` for adiabatic walls.
Explicit `qw` and `qHalf` modes are also available, but `qHalf` is rejected for
an adiabatic wall because it would impose a nonzero flux at OpenFOAM's physical
boundary face.

## LiNeumann variant: prescribed-gradient energy coupling

The library also contains a second, independently named copy of the wall
function whose energy coupling imposes the Li heat flux directly through a
Neumann (prescribed-gradient) temperature boundary condition instead of
inverting a non-negative `alphat`. Use it for walls near recovery
temperature at coarse `y+` (first cell colder than the wall), where the
`alphat` inversion is structurally unrepresentable — see
`FINDINGS_M6Tw076_energy_coupling.md` and
`DESIGN_LiNeumann_energy_coupling.md`.

- `LiNeumannWallHeatFluxTemperature` on `0/T` computes `qw` each outer
  iteration and prescribes the wall-normal temperature gradient.
- `0/alphat` becomes a plain `fixedValue 0`.
- `nutLiNeumannWallFunction` and `nuTildaLiNeumannWallFunction` are
  unchanged copies of the momentum and SA couplings using the same shared
  evaluator.

Additional entries beyond the shared defaults (all four fields must carry
identical coefficients):

```foam
Tw              300;      // physical wall temperature, required (isothermal)
wallTransport             // wall mu/kappa at (p_wall, Tw), required (isothermal)
{
    mode        sutherland;   // or constant (mu/kappa)
    As          1.458e-06;
    Ts          110.4;
    // Pr       0.72;         // optional; omitted -> modified Eucken
}
thermalFluxTarget auto;   // auto = qw (qHalf remains a legacy/diagnostic target)
// gradientRelaxation 0.3;   // optional startup aids
// TFaceMin           100;
// TFaceMax           5000;
```

**Warning:** under this coupling the patch `T` value is the
tangent-extrapolated face temperature, *not* the physical wall
temperature (e.g. ~224 K against a 300 K wall on M6Tw076 at `y+ ~ 10`).
The physical wall temperature lives in the `Tw` entry and the diagnostics
output; audit any post-processing that reads wall-patch `T` or transport
properties. `wallHeatFlux` reports the imposed flux correctly.

Diagnostics (same switches as the alphat condition, placed on `T`) are
written to `postProcessing/LiNeumannWallFunctionDiagnostics/<patch>/`.
Adiabatic mode (`thermalMode adiabatic`) reduces identically to
`zeroGradient` and needs neither `Tw` nor `wallTransport`. Complete
snippets are in [examples/boundaryConditions.foam](examples/boundaryConditions.foam).

## Mesh and flow requirements

The method requires a distinct second off-wall cell with tangential speed
larger than the first-cell speed. It is intended for smooth, attached,
approximately zero-pressure-gradient equilibrium boundary layers. It is not a
stagnation, separation, roughness, or shock-interaction model.

If a face has zero tangential velocity or no valid reference state, the default
is a laminar fallback for that face. Set `strict true` to stop instead. The
diagnostic `status` values are:

| Status | Meaning |
| ---: | --- |
| 0 | Li wall function evaluated successfully |
| 1 | Zero tangential velocity; zero-turbulence limit |
| 2 | Missing/non-monotone second-cell reference |
| 3 | Degenerate Eq. (40) velocity ratio |
| 4 | Invalid enthalpy or generalized recovery factor |
| 5 | Eq. (42) nonlinear solve did not converge |

## Diagnostics

Set the following on `alphat`:

```foam
diagnostics              true;
diagnosticLog            true;
diagnosticWrite          true;
diagnosticWriteInterval  1;
diagnosticWritePrecision 12;
```

Face data are written under:

```text
postProcessing/LiWallFunctionDiagnostics/<patch>/<time>.dat
```

Parallel runs add `processor<N>` below the diagnostics directory. The output
contains the selected reference, solve status, all principal Eq. (33)-(42)
variables, the SA first-cell values, and the requested and realised OpenFOAM
energy fluxes.

## Equation tests

```bash
source /opt/of-v2412/etc/bashrc
cd tests/equations
wmake
TestLiWallFunction
```

The test checks the Spalding limit of Eq. (42), Eq. (38) by finite difference,
the Eq. (39) SA inversion, the Eq. (40) reference-point inversion, and the
calorically perfect enthalpy limit.

The compressible integration smoke tests exercise isothermal and adiabatic
walls, the local `polynomialTransport`/JANAF internal-energy thermo, and a
two-rank decomposed run:

```bash
source /opt/of-v2412/etc/bashrc
cd tests/integration
bash Allrun
```

The generated cases and diagnostics are written under `/tmp`.

See [DESIGN_AND_DEVIATIONS.md](DESIGN_AND_DEVIATIONS.md) before using the
library for production calculations.
